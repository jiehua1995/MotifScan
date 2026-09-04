//! Species-aware fuzzy scanning for noisy long-read FASTA/FASTQ data.
//!
//! The scanner deliberately separates two questions:
//! 1. locus identity: long-window approximate matching with sequencing-error tolerance;
//! 2. species identity: voting only at substitution SNPs extracted from a mel/sim pair.
//!
//! Known differences between the paired references (substitution SNPs and pair-specific
//! indels) do not count as ordinary locus-matching errors. Pair-specific indels are not
//! used as species votes because long-read indel errors are comparatively difficult to
//! model robustly.

use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;

use aho_corasick::AhoCorasick;
use anyhow::{anyhow, bail, Context, Result};
use csv::{ReaderBuilder, Trim, Writer};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use crate::cli::SpeciesArgs;
use crate::io::{open_record_reader, Record};
use crate::motif::load_motif_file;

const CHUNK_SIZE: usize = 512;
const DEFAULT_NO_QUAL: u8 = 60;

#[derive(Debug, Clone)]
struct DiagnosticSnp {
    mel_pos: usize,
    sim_pos: usize,
    mel_base: u8,
    sim_base: u8,
}

#[derive(Debug, Clone)]
struct Locus {
    name: String,
    mel_name: String,
    sim_name: String,
    mel_seq: Vec<u8>,
    sim_seq: Vec<u8>,
    snps: Vec<DiagnosticSnp>,
    indel_columns: usize,
    /// Positions present in the mel reference but aligned to a gap in sim.
    /// They are known pair differences and are excluded from shared-site identity.
    mel_only_positions: HashSet<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum ScanStrand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone)]
struct AnchorMeta {
    locus_index: usize,
    strand: ScanStrand,
    ref_pos: usize,
}

/// Aho-Corasick patterns are unique strings. One pattern may intentionally point to
/// several loci/positions, which avoids losing candidates when shared anchors are equal.
struct AnchorIndex {
    ac: AhoCorasick,
    meta: Vec<Vec<AnchorMeta>>,
}

#[derive(Debug, Clone)]
struct AlignmentResult {
    edit_distance: usize,
    ref_to_target: Vec<Option<usize>>,
}

#[derive(Debug, Clone)]
struct ReadHit {
    locus_index: usize,
    strand: ScanStrand,
    shared_identity: f64,
    aligned_ref_bases: usize,
    edit_distance: usize,
    mel_support: u32,
    sim_support: u32,
    other_support: u32,
    informative_sites: u32,
    class: SpeciesClass,
    snp_observations: Vec<SnpObservation>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SpeciesClass {
    Mel,
    Sim,
    Ambiguous,
    Conflict,
    NoSnp,
}

#[derive(Debug, Clone)]
struct SnpObservation {
    snp_index: usize,
    category: SnpCategory,
    quality: u8,
}

#[derive(Debug, Clone, Copy)]
enum SnpCategory {
    Mel,
    Sim,
    Other,
    LowQuality,
}

#[derive(Debug, Clone, Default)]
struct SnpCount {
    mel: u64,
    sim: u64,
    other: u64,
    lowq: u64,
    q_sum: u64,
    q_n: u64,
}

#[derive(Debug, Clone, Default)]
struct LocusStats {
    reads_with_hit: u64,
    total_hits: u64,
    forward_hits: u64,
    revcomp_hits: u64,
    mel_reads: u64,
    sim_reads: u64,
    ambiguous_reads: u64,
    conflict_reads: u64,
    no_snp_reads: u64,
    shared_identity_sum: f64,
    informative_sum: u64,
    snp_counts: Vec<SnpCount>,
}

#[derive(Debug, Clone)]
struct InputSpec {
    sample: String,
    path: PathBuf,
}

pub fn run_species(args: &SpeciesArgs) -> Result<()> {
    args.validate()?;

    let motifs_path = args.motifs.as_ref().unwrap();
    let pairs_path = args.pairs.as_ref().unwrap();
    let output_path = args.output.as_ref().unwrap();

    let motif_map: HashMap<String, Vec<u8>> = load_motif_file(motifs_path)?
        .into_iter()
        .map(|m| (m.name, m.sequence.trim().to_ascii_uppercase().into_bytes()))
        .collect();
    let loci = load_pairs(pairs_path, &motif_map)?;
    if loci.is_empty() {
        bail!("pair file did not contain any locus pairs");
    }

    let anchors = build_anchor_index(&loci, args.anchor_k, args.anchors_per_locus)?;
    let inputs = collect_inputs(args)?;

    let mut summary_writer = Writer::from_path(output_path)
        .with_context(|| format!("failed to create {}", output_path.display()))?;
    write_summary_header(&mut summary_writer)?;

    let snp_path = args.snp_output.clone().unwrap_or_else(|| {
        let mut p = output_path.to_path_buf();
        let stem = output_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("species_scan");
        p.set_file_name(format!("{stem}.snps.csv"));
        p
    });
    let mut snp_writer = Writer::from_path(&snp_path)
        .with_context(|| format!("failed to create {}", snp_path.display()))?;
    write_snp_header(&mut snp_writer)?;

    if let Some(path) = &args.pair_qc_output {
        write_pair_qc(path, &loci)?;
    }

    for (idx, input) in inputs.iter().enumerate() {
        eprintln!(
            "[{}/{}] sample={} input={}",
            idx + 1,
            inputs.len(),
            input.sample,
            input.path.display()
        );
        let stats = scan_one_file(input, &loci, &anchors, args)?;
        write_sample_summary(&mut summary_writer, input, &loci, &stats)?;
        write_sample_snps(&mut snp_writer, input, &loci, &stats)?;
        summary_writer.flush()?;
        snp_writer.flush()?;
    }

    eprintln!("summary: {}", output_path.display());
    eprintln!("SNP QC : {}", snp_path.display());
    Ok(())
}

fn collect_inputs(args: &SpeciesArgs) -> Result<Vec<InputSpec>> {
    let mut out = Vec::new();
    if let Some(input) = &args.input {
        out.push(InputSpec {
            sample: args
                .sample
                .clone()
                .unwrap_or_else(|| infer_sample_name(input)),
            path: input.clone(),
        });
    }
    if let Some(list_path) = &args.input_list {
        let reader = BufReader::new(
            File::open(list_path)
                .with_context(|| format!("failed to open {}", list_path.display()))?,
        );
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = trimmed.split('\t').collect();
            let (sample, path) = if fields.len() >= 2 {
                (fields[0].trim().to_string(), PathBuf::from(fields[1].trim()))
            } else {
                let p = PathBuf::from(trimmed);
                (infer_sample_name(&p), p)
            };
            if sample.is_empty() {
                bail!("empty sample name at {}:{}", list_path.display(), line_no + 1);
            }
            out.push(InputSpec { sample, path });
        }
    }
    if out.is_empty() {
        bail!("one of --input or --input-list is required");
    }
    for x in &out {
        if !x.path.is_file() {
            bail!("input file not found: {}", x.path.display());
        }
    }
    Ok(out)
}

fn infer_sample_name(path: &Path) -> String {
    let name = path.file_name().and_then(|s| s.to_str()).unwrap_or("sample");
    for suffix in [
        ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fastq", ".fq", ".fasta", ".fa",
    ] {
        if name.to_ascii_lowercase().ends_with(suffix) {
            return name[..name.len() - suffix.len()].to_string();
        }
    }
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("sample")
        .to_string()
}

fn load_pairs(path: &Path, motifs: &HashMap<String, Vec<u8>>) -> Result<Vec<Locus>> {
    let mut reader = ReaderBuilder::new()
        .has_headers(false)
        .comment(Some(b'#'))
        .trim(Trim::All)
        .from_path(path)
        .with_context(|| format!("failed to open pair file {}", path.display()))?;

    let mut loci = Vec::new();
    for (row_idx, rec) in reader.records().enumerate() {
        let rec = rec?;
        if rec.is_empty() {
            continue;
        }
        if rec.len() != 3 {
            bail!(
                "{}:{} expected 3 columns: locus,mel,sim",
                path.display(),
                row_idx + 1
            );
        }
        let name = rec.get(0).unwrap().trim();
        let mel_name = rec.get(1).unwrap().trim();
        let sim_name = rec.get(2).unwrap().trim();
        if row_idx == 0
            && name.eq_ignore_ascii_case("locus")
            && mel_name.eq_ignore_ascii_case("mel")
            && sim_name.eq_ignore_ascii_case("sim")
        {
            continue;
        }
        let mel_seq = motifs
            .get(mel_name)
            .ok_or_else(|| anyhow!("pair '{}' references missing motif '{}'", name, mel_name))?
            .clone();
        let sim_seq = motifs
            .get(sim_name)
            .ok_or_else(|| anyhow!("pair '{}' references missing motif '{}'", name, sim_name))?
            .clone();
        validate_dna(mel_name, &mel_seq)?;
        validate_dna(sim_name, &sim_seq)?;
        let (snps, indel_columns, mel_only_positions) = extract_pair_differences(&mel_seq, &sim_seq);
        if snps.is_empty() {
            bail!("pair '{}' has no diagnostic substitution SNPs", name);
        }
        loci.push(Locus {
            name: name.to_string(),
            mel_name: mel_name.to_string(),
            sim_name: sim_name.to_string(),
            mel_seq,
            sim_seq,
            snps,
            indel_columns,
            mel_only_positions,
        });
    }
    Ok(loci)
}

fn validate_dna(name: &str, seq: &[u8]) -> Result<()> {
    if seq.is_empty() {
        bail!("motif '{}' has empty sequence", name);
    }
    if let Some(&bad) = seq
        .iter()
        .find(|&&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
    {
        bail!("motif '{}' contains unsupported base '{}'", name, bad as char);
    }
    Ok(())
}

/// Return substitution SNPs, number of pairwise indel columns, and mel coordinates that
/// are absent from the sim reference. The latter are excluded from locus identity.
fn extract_pair_differences(
    mel: &[u8],
    sim: &[u8],
) -> (Vec<DiagnosticSnp>, usize, HashSet<usize>) {
    let (a, b) = global_align(mel, sim);
    let mut mel_pos = 0usize;
    let mut sim_pos = 0usize;
    let mut snps = Vec::new();
    let mut indel_columns = 0usize;
    let mut mel_only_positions = HashSet::new();

    for (&x, &y) in a.iter().zip(b.iter()) {
        match (x, y) {
            (Some(mx), Some(sx)) => {
                if mx != sx {
                    snps.push(DiagnosticSnp {
                        mel_pos,
                        sim_pos,
                        mel_base: mx,
                        sim_base: sx,
                    });
                }
                mel_pos += 1;
                sim_pos += 1;
            }
            (Some(_), None) => {
                mel_only_positions.insert(mel_pos);
                mel_pos += 1;
                indel_columns += 1;
            }
            (None, Some(_)) => {
                sim_pos += 1;
                indel_columns += 1;
            }
            (None, None) => unreachable!(),
        }
    }
    (snps, indel_columns, mel_only_positions)
}

/// Simple edit-cost global alignment. References are short diagnostic windows, so an
/// O(m*n) matrix is tiny and this runs only once per pair at startup.
fn global_align(a: &[u8], b: &[u8]) -> (Vec<Option<u8>>, Vec<Option<u8>>) {
    let mut dp = vec![vec![0usize; b.len() + 1]; a.len() + 1];
    for i in 0..=a.len() {
        dp[i][0] = i;
    }
    for j in 0..=b.len() {
        dp[0][j] = j;
    }
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let sub = dp[i - 1][j - 1] + usize::from(a[i - 1] != b[j - 1]);
            dp[i][j] = sub.min(dp[i - 1][j] + 1).min(dp[i][j - 1] + 1);
        }
    }

    let (mut i, mut j) = (a.len(), b.len());
    let mut aa = Vec::new();
    let mut bb = Vec::new();
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            if dp[i][j] == dp[i - 1][j - 1] + cost {
                aa.push(Some(a[i - 1]));
                bb.push(Some(b[j - 1]));
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i][j] == dp[i - 1][j] + 1 {
            aa.push(Some(a[i - 1]));
            bb.push(None);
            i -= 1;
        } else {
            aa.push(None);
            bb.push(Some(b[j - 1]));
            j -= 1;
        }
    }
    aa.reverse();
    bb.reverse();
    (aa, bb)
}

fn build_anchor_index(loci: &[Locus], k: usize, per_locus: usize) -> Result<AnchorIndex> {
    if k < 5 {
        bail!("--anchor-k must be at least 5");
    }
    if per_locus == 0 {
        bail!("--anchors-per-locus must be greater than 0");
    }

    let mut patterns = Vec::<String>::new();
    let mut meta = Vec::<Vec<AnchorMeta>>::new();
    let mut pattern_ids = HashMap::<String, usize>::new();

    let mut register = |pattern: String, anchor_meta: AnchorMeta| {
        if let Some(&id) = pattern_ids.get(&pattern) {
            meta[id].push(anchor_meta);
        } else {
            let id = patterns.len();
            pattern_ids.insert(pattern.clone(), id);
            patterns.push(pattern);
            meta.push(vec![anchor_meta]);
        }
    };

    for (locus_index, locus) in loci.iter().enumerate() {
        if locus.mel_seq.len() < k {
            bail!("locus '{}' is shorter than --anchor-k", locus.name);
        }
        let snp_positions: HashSet<usize> = locus.snps.iter().map(|s| s.mel_pos).collect();
        let mut candidates = Vec::new();
        for pos in 0..=locus.mel_seq.len() - k {
            if (pos..pos + k).any(|p| {
                snp_positions.contains(&p) || locus.mel_only_positions.contains(&p)
            }) {
                continue;
            }
            let mer = &locus.mel_seq[pos..pos + k];
            if mer.contains(&b'N') || !contains_subslice(&locus.sim_seq, mer) {
                continue;
            }
            candidates.push((pos, mer.to_vec()));
        }
        if candidates.is_empty() {
            bail!(
                "locus '{}' has no shared {}-mer anchors; reduce --anchor-k",
                locus.name,
                k
            );
        }

        for (pos, mer) in evenly_spaced(&candidates, per_locus) {
            register(
                String::from_utf8(mer.clone()).unwrap(),
                AnchorMeta {
                    locus_index,
                    strand: ScanStrand::Forward,
                    ref_pos: pos,
                },
            );
            register(
                String::from_utf8(revcomp(&mer)).unwrap(),
                AnchorMeta {
                    locus_index,
                    strand: ScanStrand::Reverse,
                    ref_pos: locus.mel_seq.len() - (pos + k),
                },
            );
        }
    }

    let refs: Vec<&str> = patterns.iter().map(String::as_str).collect();
    Ok(AnchorIndex {
        ac: AhoCorasick::new(&refs),
        meta,
    })
}

fn evenly_spaced(candidates: &[(usize, Vec<u8>)], n: usize) -> Vec<(usize, Vec<u8>)> {
    if candidates.len() <= n {
        return candidates.to_vec();
    }
    if n == 1 {
        return vec![candidates[candidates.len() / 2].clone()];
    }
    (0..n)
        .map(|i| candidates[i * (candidates.len() - 1) / (n - 1)].clone())
        .collect()
}

fn contains_subslice(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|x| x == needle)
}

fn scan_one_file(
    input: &InputSpec,
    loci: &[Locus],
    anchors: &AnchorIndex,
    args: &SpeciesArgs,
) -> Result<Vec<LocusStats>> {
    let mut reader = open_record_reader(&input.path)?;
    let snapshot = reader.progress_snapshot();
    let bar = if args.progress {
        let pb = ProgressBar::new(snapshot.total_bytes);
        pb.set_style(ProgressStyle::with_template(
            "{spinner:.green} {msg}\n[{bar:40.cyan/blue}] {percent:>3}% | {bytes}/{total_bytes} | eta {eta_precise}",
        )?);
        Some(pb)
    } else {
        None
    };

    let started = Instant::now();
    let mut reads_processed = 0u64;
    let mut bases_processed = 0u64;
    let mut stats: Vec<LocusStats> = loci
        .iter()
        .map(|l| LocusStats {
            snp_counts: vec![SnpCount::default(); l.snps.len()],
            ..Default::default()
        })
        .collect();

    loop {
        let chunk = reader.next_chunk(CHUNK_SIZE)?;
        if chunk.is_empty() {
            break;
        }
        reads_processed += chunk.len() as u64;
        bases_processed += chunk.iter().map(|r| r.seq.len() as u64).sum::<u64>();

        let batch: Vec<Vec<ReadHit>> = chunk
            .par_iter()
            .map(|r| scan_record_species(r, loci, anchors, args))
            .collect();
        for hits in batch {
            for hit in hits {
                merge_hit(&mut stats[hit.locus_index], &hit);
            }
        }

        if let Some(pb) = &bar {
            let snap = reader.progress_snapshot();
            pb.set_position(snap.bytes_read.min(snap.total_bytes));
            let elapsed = started.elapsed().as_secs_f64().max(1e-9);
            let avg = if reads_processed > 0 {
                bases_processed as f64 / reads_processed as f64
            } else {
                0.0
            };
            pb.set_message(format!(
                "{} | reads {} | {:.0} reads/s | avg {:.0} bp",
                input.sample,
                reads_processed,
                reads_processed as f64 / elapsed,
                avg
            ));
        }
    }

    if let Some(pb) = &bar {
        pb.finish_and_clear();
    }
    eprintln!(
        "  done: reads={} bases={} elapsed={:.1}s",
        reads_processed,
        bases_processed,
        started.elapsed().as_secs_f64()
    );
    Ok(stats)
}

fn scan_record_species(
    record: &Record,
    loci: &[Locus],
    anchors: &AnchorIndex,
    args: &SpeciesArgs,
) -> Vec<ReadHit> {
    let mut candidates: HashMap<(usize, ScanStrand), Vec<isize>> = HashMap::new();
    for m in anchors.ac.find_iter(&record.seq) {
        for meta in &anchors.meta[m.pattern()] {
            let estimated_start = m.start() as isize - meta.ref_pos as isize;
            let starts = candidates.entry((meta.locus_index, meta.strand)).or_default();
            if !starts.iter().any(|&x| (x - estimated_start).abs() <= 3) {
                starts.push(estimated_start);
            }
        }
    }

    let mut best_by_locus: HashMap<usize, ReadHit> = HashMap::new();
    for ((locus_index, strand), starts) in candidates {
        for estimated_start in starts {
            if let Some(hit) = evaluate_candidate(
                locus_index,
                &loci[locus_index],
                strand,
                estimated_start,
                &record.seq,
                record.qual.as_deref(),
                args,
            ) {
                let replace = best_by_locus
                    .get(&locus_index)
                    .map(|old| hit_key(&hit) > hit_key(old))
                    .unwrap_or(true);
                if replace {
                    best_by_locus.insert(locus_index, hit);
                }
            }
        }
    }

    let mut hits: Vec<ReadHit> = best_by_locus.into_values().collect();
    if args.locus_mode == "all" {
        return hits;
    }
    hits.sort_by(|a, b| hit_key(b).partial_cmp(&hit_key(a)).unwrap());
    hits.truncate(1);
    hits
}

fn hit_key(hit: &ReadHit) -> (f64, usize, u32, Reverse<usize>) {
    (
        hit.shared_identity,
        hit.aligned_ref_bases,
        hit.informative_sites,
        Reverse(hit.edit_distance),
    )
}

fn evaluate_candidate(
    locus_index: usize,
    locus: &Locus,
    strand: ScanStrand,
    estimated_start: isize,
    read: &[u8],
    qual: Option<&[u8]>,
    args: &SpeciesArgs,
) -> Option<ReadHit> {
    let reference = match strand {
        ScanStrand::Forward => locus.mel_seq.clone(),
        ScanStrand::Reverse => revcomp(&locus.mel_seq),
    };
    let slack = args.alignment_slack as isize;
    let from = (estimated_start - slack).max(0) as usize;
    let to = (estimated_start + reference.len() as isize + slack)
        .max(0)
        .min(read.len() as isize) as usize;
    if to <= from {
        return None;
    }
    let target = &read[from..to];
    let alignment = semi_global_align(&reference, target);

    let mut snp_by_pos: HashMap<usize, (usize, u8, u8)> = HashMap::new();
    let mut excluded_pair_indels = HashSet::new();
    for (idx, snp) in locus.snps.iter().enumerate() {
        let (pos, mel_base, sim_base) = match strand {
            ScanStrand::Forward => (snp.mel_pos, snp.mel_base, snp.sim_base),
            ScanStrand::Reverse => (
                locus.mel_seq.len() - 1 - snp.mel_pos,
                complement(snp.mel_base),
                complement(snp.sim_base),
            ),
        };
        snp_by_pos.insert(pos, (idx, mel_base, sim_base));
    }
    for &pos in &locus.mel_only_positions {
        excluded_pair_indels.insert(match strand {
            ScanStrand::Forward => pos,
            ScanStrand::Reverse => locus.mel_seq.len() - 1 - pos,
        });
    }

    let mut shared_match = 0usize;
    let mut shared_error = 0usize;
    let mut aligned_ref_bases = 0usize;
    let mut mel_support = 0u32;
    let mut sim_support = 0u32;
    let mut other_support = 0u32;
    let mut observations = Vec::new();

    for ref_pos in 0..reference.len() {
        if excluded_pair_indels.contains(&ref_pos) {
            continue;
        }
        let target_pos = alignment.ref_to_target[ref_pos];
        if let Some(&(snp_idx, mel_base, sim_base)) = snp_by_pos.get(&ref_pos) {
            let Some(tp) = target_pos else { continue };
            aligned_ref_bases += 1;
            let read_pos = from + tp;
            let base = read[read_pos];
            // io.rs already decodes FASTQ Phred+33 to numeric Phred values.
            let q = qual
                .and_then(|qv| qv.get(read_pos).copied())
                .unwrap_or(DEFAULT_NO_QUAL);
            let category = if q < args.min_snp_baseq {
                SnpCategory::LowQuality
            } else if base == mel_base {
                mel_support += 1;
                SnpCategory::Mel
            } else if base == sim_base {
                sim_support += 1;
                SnpCategory::Sim
            } else {
                other_support += 1;
                SnpCategory::Other
            };
            observations.push(SnpObservation {
                snp_index: snp_idx,
                category,
                quality: q,
            });
            continue;
        }

        match target_pos {
            Some(tp) => {
                aligned_ref_bases += 1;
                if reference[ref_pos] == target[tp] {
                    shared_match += 1;
                } else {
                    shared_error += 1;
                }
            }
            None => shared_error += 1,
        }
    }

    if aligned_ref_bases < args.min_aligned_bases {
        return None;
    }
    let denominator = shared_match + shared_error;
    if denominator == 0 {
        return None;
    }
    let shared_identity = shared_match as f64 / denominator as f64;
    if shared_identity < args.min_shared_identity {
        return None;
    }

    let informative = mel_support + sim_support + other_support;
    let class = if informative < args.min_informative_snps as u32 {
        SpeciesClass::NoSnp
    } else {
        let mel_fraction = mel_support as f64 / informative as f64;
        let sim_fraction = sim_support as f64 / informative as f64;
        if mel_fraction >= args.species_fraction {
            SpeciesClass::Mel
        } else if sim_fraction >= args.species_fraction {
            SpeciesClass::Sim
        } else if mel_support > 0 && sim_support > 0 {
            SpeciesClass::Conflict
        } else {
            SpeciesClass::Ambiguous
        }
    };

    Some(ReadHit {
        locus_index,
        strand,
        shared_identity,
        aligned_ref_bases,
        edit_distance: alignment.edit_distance,
        mel_support,
        sim_support,
        other_support,
        informative_sites: informative,
        class,
        snp_observations: observations,
    })
}

/// Full reference vs best target substring (HW-style semi-global edit alignment).
/// Insertions in the read are represented by target positions not mapped to a reference
/// coordinate. They are intentionally not included in shared-site identity because they
/// may correspond either to long-read insertion error or to a known sim-only pair indel.
fn semi_global_align(reference: &[u8], target: &[u8]) -> AlignmentResult {
    let mut dp = vec![vec![0usize; target.len() + 1]; reference.len() + 1];
    for i in 0..=reference.len() {
        dp[i][0] = i;
    }
    for j in 0..=target.len() {
        dp[0][j] = 0;
    }
    for i in 1..=reference.len() {
        for j in 1..=target.len() {
            let sub = dp[i - 1][j - 1] + usize::from(reference[i - 1] != target[j - 1]);
            dp[i][j] = sub.min(dp[i - 1][j] + 1).min(dp[i][j - 1] + 1);
        }
    }

    let (mut j, &best) = dp[reference.len()]
        .iter()
        .enumerate()
        .min_by_key(|(_, v)| *v)
        .unwrap();
    let mut i = reference.len();
    let mut ref_to_target = vec![None; reference.len()];
    while i > 0 {
        if j > 0 {
            let cost = usize::from(reference[i - 1] != target[j - 1]);
            if dp[i][j] == dp[i - 1][j - 1] + cost {
                ref_to_target[i - 1] = Some(j - 1);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if dp[i][j] == dp[i - 1][j] + 1 {
            i -= 1;
        } else if j > 0 {
            j -= 1;
        } else {
            i -= 1;
        }
    }
    AlignmentResult {
        edit_distance: best,
        ref_to_target,
    }
}

fn merge_hit(stats: &mut LocusStats, hit: &ReadHit) {
    stats.reads_with_hit += 1;
    stats.total_hits += 1;
    match hit.strand {
        ScanStrand::Forward => stats.forward_hits += 1,
        ScanStrand::Reverse => stats.revcomp_hits += 1,
    }
    match hit.class {
        SpeciesClass::Mel => stats.mel_reads += 1,
        SpeciesClass::Sim => stats.sim_reads += 1,
        SpeciesClass::Ambiguous => stats.ambiguous_reads += 1,
        SpeciesClass::Conflict => stats.conflict_reads += 1,
        SpeciesClass::NoSnp => stats.no_snp_reads += 1,
    }
    stats.shared_identity_sum += hit.shared_identity;
    stats.informative_sum += hit.informative_sites as u64;
    for obs in &hit.snp_observations {
        let c = &mut stats.snp_counts[obs.snp_index];
        match obs.category {
            SnpCategory::Mel => c.mel += 1,
            SnpCategory::Sim => c.sim += 1,
            SnpCategory::Other => c.other += 1,
            SnpCategory::LowQuality => {
                c.lowq += 1;
                continue;
            }
        }
        c.q_sum += obs.quality as u64;
        c.q_n += 1;
    }
}

fn write_summary_header(w: &mut Writer<File>) -> Result<()> {
    w.write_record([
        "sample",
        "input_file",
        "locus",
        "mel_motif",
        "sim_motif",
        "mel_sequence",
        "sim_sequence",
        "mel_length",
        "sim_length",
        "diagnostic_snp_count",
        "pairwise_indel_columns",
        "diagnostic_snps",
        "reads_with_hit",
        "total_hits",
        "forward_hits",
        "revcomp_hits",
        "mel_reads",
        "sim_reads",
        "ambiguous_reads",
        "conflict_reads",
        "no_snp_reads",
        "species_assigned_reads",
        "mel_fraction_among_assigned",
        "sim_fraction_among_assigned",
        "mean_shared_identity",
        "mean_informative_snps_per_hit",
    ])?;
    Ok(())
}

fn write_snp_header(w: &mut Writer<File>) -> Result<()> {
    w.write_record([
        "sample",
        "input_file",
        "locus",
        "snp_index",
        "mel_pos_1based",
        "sim_pos_1based",
        "mel_base",
        "sim_base",
        "mel_count",
        "sim_count",
        "other_count",
        "lowq_count",
        "hq_depth",
        "mel_fraction",
        "sim_fraction",
        "other_fraction",
        "mean_hq_baseq",
    ])?;
    Ok(())
}

fn diagnostic_string(locus: &Locus) -> String {
    locus
        .snps
        .iter()
        .enumerate()
        .map(|(i, s)| {
            format!(
                "{}:mel{}{}>sim{}{}",
                i + 1,
                s.mel_pos + 1,
                s.mel_base as char,
                s.sim_pos + 1,
                s.sim_base as char
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn write_sample_summary(
    w: &mut Writer<File>,
    input: &InputSpec,
    loci: &[Locus],
    stats: &[LocusStats],
) -> Result<()> {
    for (locus, st) in loci.iter().zip(stats) {
        let assigned = st.mel_reads + st.sim_reads;
        let mel_fraction = if assigned > 0 {
            format!("{:.8}", st.mel_reads as f64 / assigned as f64)
        } else {
            String::new()
        };
        let sim_fraction = if assigned > 0 {
            format!("{:.8}", st.sim_reads as f64 / assigned as f64)
        } else {
            String::new()
        };
        let mean_identity = if st.reads_with_hit > 0 {
            format!("{:.6}", st.shared_identity_sum / st.reads_with_hit as f64)
        } else {
            String::new()
        };
        let mean_snps = if st.reads_with_hit > 0 {
            format!("{:.4}", st.informative_sum as f64 / st.reads_with_hit as f64)
        } else {
            String::new()
        };
        w.write_record([
            input.sample.clone(),
            input.path.display().to_string(),
            locus.name.clone(),
            locus.mel_name.clone(),
            locus.sim_name.clone(),
            String::from_utf8_lossy(&locus.mel_seq).to_string(),
            String::from_utf8_lossy(&locus.sim_seq).to_string(),
            locus.mel_seq.len().to_string(),
            locus.sim_seq.len().to_string(),
            locus.snps.len().to_string(),
            locus.indel_columns.to_string(),
            diagnostic_string(locus),
            st.reads_with_hit.to_string(),
            st.total_hits.to_string(),
            st.forward_hits.to_string(),
            st.revcomp_hits.to_string(),
            st.mel_reads.to_string(),
            st.sim_reads.to_string(),
            st.ambiguous_reads.to_string(),
            st.conflict_reads.to_string(),
            st.no_snp_reads.to_string(),
            assigned.to_string(),
            mel_fraction,
            sim_fraction,
            mean_identity,
            mean_snps,
        ])?;
    }
    Ok(())
}

fn write_sample_snps(
    w: &mut Writer<File>,
    input: &InputSpec,
    loci: &[Locus],
    stats: &[LocusStats],
) -> Result<()> {
    for (locus, st) in loci.iter().zip(stats) {
        for (idx, snp) in locus.snps.iter().enumerate() {
            let c = &st.snp_counts[idx];
            let depth = c.mel + c.sim + c.other;
            let fraction = |x: u64| {
                if depth > 0 {
                    format!("{:.8}", x as f64 / depth as f64)
                } else {
                    String::new()
                }
            };
            let mean_q = if c.q_n > 0 {
                format!("{:.3}", c.q_sum as f64 / c.q_n as f64)
            } else {
                String::new()
            };
            w.write_record([
                input.sample.clone(),
                input.path.display().to_string(),
                locus.name.clone(),
                (idx + 1).to_string(),
                (snp.mel_pos + 1).to_string(),
                (snp.sim_pos + 1).to_string(),
                (snp.mel_base as char).to_string(),
                (snp.sim_base as char).to_string(),
                c.mel.to_string(),
                c.sim.to_string(),
                c.other.to_string(),
                c.lowq.to_string(),
                depth.to_string(),
                fraction(c.mel),
                fraction(c.sim),
                fraction(c.other),
                mean_q,
            ])?;
        }
    }
    Ok(())
}

fn write_pair_qc(path: &Path, loci: &[Locus]) -> Result<()> {
    let mut w = Writer::from_path(path)?;
    w.write_record([
        "locus",
        "mel_motif",
        "sim_motif",
        "mel_length",
        "sim_length",
        "diagnostic_snp_count",
        "pairwise_indel_columns",
        "diagnostic_snps",
    ])?;
    for locus in loci {
        w.write_record([
            locus.name.clone(),
            locus.mel_name.clone(),
            locus.sim_name.clone(),
            locus.mel_seq.len().to_string(),
            locus.sim_seq.len().to_string(),
            locus.snps.len().to_string(),
            locus.indel_columns.to_string(),
            diagnostic_string(locus),
        ])?;
    }
    w.flush()?;
    Ok(())
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy_locus() -> Locus {
        let mel = b"ACGTACGTACGTACGTACGTACGT".to_vec();
        let mut sim = mel.clone();
        sim[5] = b'T';
        sim[12] = b'G';
        sim[19] = b'A';
        let (snps, indel_columns, mel_only_positions) = extract_pair_differences(&mel, &sim);
        Locus {
            name: "toy".into(),
            mel_name: "mel".into(),
            sim_name: "sim".into(),
            mel_seq: mel,
            sim_seq: sim,
            snps,
            indel_columns,
            mel_only_positions,
        }
    }

    fn toy_args() -> SpeciesArgs {
        SpeciesArgs {
            input: None,
            input_list: None,
            sample: None,
            motifs: None,
            pairs: None,
            output: None,
            snp_output: None,
            pair_qc_output: None,
            threads: 1,
            progress: false,
            revcomp: true,
            anchor_k: 5,
            anchors_per_locus: 4,
            alignment_slack: 8,
            min_shared_identity: 0.80,
            min_aligned_bases: 15,
            min_snp_baseq: 10,
            min_informative_snps: 2,
            species_fraction: 0.75,
            locus_mode: "best".into(),
            help: false,
        }
    }

    #[test]
    fn extracts_expected_substitution_snps() {
        let locus = toy_locus();
        assert_eq!(locus.snps.len(), 3);
        assert_eq!(locus.indel_columns, 0);
        assert!(locus.mel_only_positions.is_empty());
        assert_eq!(locus.snps[0].mel_pos, 5);
    }

    #[test]
    fn semi_global_tolerates_shared_substitution_and_indel() {
        let a = semi_global_align(b"ACGTACGTACGT", b"TTACGTTCGTTACGTGG");
        assert!(a.edit_distance <= 3);
        assert_eq!(a.ref_to_target.len(), 12);
    }

    #[test]
    fn mel_and_sim_votes_are_separated() {
        let locus = toy_locus();
        let args = toy_args();
        let mel = evaluate_candidate(
            0,
            &locus,
            ScanStrand::Forward,
            0,
            &locus.mel_seq,
            None,
            &args,
        )
        .unwrap();
        assert_eq!(mel.class, SpeciesClass::Mel);
        assert_eq!(mel.mel_support, 3);
        let sim = evaluate_candidate(
            0,
            &locus,
            ScanStrand::Forward,
            0,
            &locus.sim_seq,
            None,
            &args,
        )
        .unwrap();
        assert_eq!(sim.class, SpeciesClass::Sim);
        assert_eq!(sim.sim_support, 3);
    }

    #[test]
    fn shared_errors_do_not_change_species_call() {
        let locus = toy_locus();
        let args = toy_args();
        let mut read = locus.mel_seq.clone();
        read[2] = b'A';
        read.insert(9, b'T');
        let hit = evaluate_candidate(0, &locus, ScanStrand::Forward, 0, &read, None, &args)
            .unwrap();
        assert_eq!(hit.class, SpeciesClass::Mel);
    }

    #[test]
    fn reverse_complement_is_classified() {
        let locus = toy_locus();
        let args = toy_args();
        let read = revcomp(&locus.sim_seq);
        let hit = evaluate_candidate(0, &locus, ScanStrand::Reverse, 0, &read, None, &args)
            .unwrap();
        assert_eq!(hit.class, SpeciesClass::Sim);
    }

    #[test]
    fn known_pair_indels_do_not_penalize_sim_identity() {
        let mel = b"AACCGGTTAACCGGTTAACCGGTT".to_vec();
        let mut sim = mel.clone();
        sim[4] = b'T';
        sim[15] = b'A';
        sim.drain(19..22);
        let (snps, indel_columns, mel_only_positions) = extract_pair_differences(&mel, &sim);
        let locus = Locus {
            name: "indel".into(),
            mel_name: "mel".into(),
            sim_name: "sim".into(),
            mel_seq: mel,
            sim_seq: sim.clone(),
            snps,
            indel_columns,
            mel_only_positions,
        };
        let mut args = toy_args();
        args.min_shared_identity = 0.95;
        args.min_aligned_bases = 10;
        let hit = evaluate_candidate(0, &locus, ScanStrand::Forward, 0, &sim, None, &args)
            .unwrap();
        assert_eq!(hit.class, SpeciesClass::Sim);
        assert!((hit.shared_identity - 1.0).abs() < 1e-9);
    }
}
