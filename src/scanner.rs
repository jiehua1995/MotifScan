use std::time::Duration;

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use memchr::memchr_iter;
use rayon::prelude::*;

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use crate::cli::CountArgs;
use crate::io::{open_record_reader, ProgressSnapshot, Record, RecordReader};
use crate::motif::{
    canonical_base_mask, compile_motifs, load_motif_file, load_single_motif, CompiledMotif,
    Strand,
};
use crate::output::{
    create_writer, write_count_summary, write_read_hit_headers, write_read_hit_rows, CountRow,
    ReadHitRow, TableWriter,
};

const DEFAULT_CHUNK_SIZE: usize = 512;

#[derive(Debug, Clone, Copy)]
struct ScanOptions {
    use_iupac: bool,
}

#[derive(Debug, Clone)]
struct RecordResult {
    motif_hits: Vec<MotifHitSummary>,
    read_hits: Vec<ReadHitRow>,
}

#[derive(Debug, Clone)]
struct MotifHitSummary {
    motif_index: usize,
    total_hits: u64,
    forward_hits: u64,
    revcomp_hits: u64,
    read_has_hit: bool,
}

#[derive(Debug, Clone, Default)]
struct PatternScanResult {
    hit_count: u64,
    positions: Vec<usize>,
}

pub fn run_count(args: &CountArgs) -> Result<()> {
    args.validate()?;
    let raw_motifs = if let Some(sequence) = &args.motif {
        load_single_motif(&args.motif_name, sequence)?
    } else {
        load_motif_file(args.motifs.as_ref().unwrap())?
    };
    let motifs = compile_motifs(&raw_motifs, args.revcomp, args.iupac)?;
    let scan_options = ScanOptions {
        use_iupac: args.iupac,
    };

    let mut reader = open_record_reader(&args.input)?;
    let mut rows = initialize_rows(&motifs);
    let mut hit_writer = maybe_open_hit_writer(args.report_read_hits.as_deref())?;
    let mut progress =
        ScanProgress::new(&reader, args.progress, "count", &args.input, motifs.len())?;

    scan_records(
        &mut reader,
        &motifs,
        &mut progress,
        scan_options,
        hit_writer.as_mut(),
        &mut rows,
    )?;

    if let Some(writer) = &mut hit_writer {
        writer.flush()?;
    }
    progress.finish();
    write_count_summary(&args.output, &rows)
}

fn maybe_open_hit_writer(path: Option<&std::path::Path>) -> Result<Option<TableWriter>> {
    match path {
        Some(path) => {
            let mut writer = create_writer(path)?;
            write_read_hit_headers(&mut writer)?;
            Ok(Some(writer))
        }
        None => Ok(None),
    }
}

fn initialize_rows(motifs: &[CompiledMotif]) -> Vec<CountRow> {
    motifs
        .iter()
        .map(|motif| CountRow {
            motif: motif.name.clone(),
            sequence: motif.sequence.clone(),
            length: motif.len(),
            reads_with_hit: 0,
            total_hits: 0,
            forward_hits: 0,
            revcomp_hits: 0,
        })
        .collect()
}

fn scan_records(
    reader: &mut RecordReader,
    motifs: &[CompiledMotif],
    progress: &mut ScanProgress,
    options: ScanOptions,
    mut hit_writer: Option<&mut TableWriter>,
    rows: &mut [CountRow],
) -> Result<()> {
    loop {
        let chunk = reader.next_chunk(DEFAULT_CHUNK_SIZE)?;
        if chunk.is_empty() {
            break;
        }
        let chunk_reads = chunk.len() as u64;
        let chunk_bases = chunk.iter().map(|record| record.seq.len() as u64).sum();

        let emit_read_hits = hit_writer.is_some();
        let record_results: Vec<RecordResult> = chunk
            .into_par_iter()
            .map(|record| scan_record(&record, motifs, options, emit_read_hits))
            .collect();

        let mut chunk_read_hits = Vec::new();
        for record_result in record_results {
            merge_record_result(&record_result, rows);
            if emit_read_hits {
                chunk_read_hits.extend(record_result.read_hits);
            }
        }

        if let Some(writer) = hit_writer.as_deref_mut() {
            write_read_hit_rows(writer, &chunk_read_hits)?;
        }

        progress.update(chunk_reads, chunk_bases, reader.progress_snapshot());
    }

    Ok(())
}

enum ScanProgress {
    Enabled(ProgressState),
    Disabled,
}

struct ProgressState {
    bar: ProgressBar,
    reads_processed: u64,
    bases_processed: u64,
    input_name: String,
    mode: &'static str,
    motif_count: usize,
}

impl ScanProgress {
    fn new(
        reader: &RecordReader,
        enabled: bool,
        mode: &'static str,
        input_path: &std::path::Path,
        motif_count: usize,
    ) -> Result<Self> {
        if !enabled {
            return Ok(Self::Disabled);
        }

        let snapshot = reader.progress_snapshot();
        let bar = ProgressBar::new(snapshot.total_bytes);
        let style = ProgressStyle::with_template(
            "{spinner:.green} {msg} [{bar:36.cyan/blue}] {percent:>3}% | {bytes}/{total_bytes} | elapsed {elapsed_precise} | eta {eta_precise}",
        )?;
        bar.set_style(style.progress_chars("##-"));
        bar.enable_steady_tick(Duration::from_millis(120));
        bar.set_position(snapshot.bytes_read);
        let input_name = input_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("input")
            .to_string();
        let mut state = ProgressState {
            bar,
            reads_processed: 0,
            bases_processed: 0,
            input_name,
            mode,
            motif_count,
        };
        state.refresh_message();
        Ok(Self::Enabled(state))
    }

    fn update(&mut self, chunk_reads: u64, chunk_bases: u64, snapshot: ProgressSnapshot) {
        if let Self::Enabled(state) = self {
            state.reads_processed += chunk_reads;
            state.bases_processed += chunk_bases;
            state
                .bar
                .set_position(snapshot.bytes_read.min(snapshot.total_bytes));
            state.refresh_message();
        }
    }

    fn finish(&self) {
        if let Self::Enabled(state) = self {
            state.bar.finish_and_clear();
        }
    }
}

impl ProgressState {
    fn refresh_message(&mut self) {
        let elapsed = self.bar.elapsed().as_secs_f64();
        let reads_per_sec = if elapsed > 0.0 {
            self.reads_processed as f64 / elapsed
        } else {
            0.0
        };
        let avg_read_len = if self.reads_processed > 0 {
            self.bases_processed as f64 / self.reads_processed as f64
        } else {
            0.0
        };
        self.bar.set_message(format!(
            "{} {} | motifs {} | reads {} | avg_len {:.1} bp | {:.1} reads/s",
            self.mode,
            self.input_name,
            self.motif_count,
            self.reads_processed,
            avg_read_len,
            reads_per_sec
        ));
    }
}

fn merge_record_result(record_result: &RecordResult, rows: &mut [CountRow]) {
    for motif_hit in &record_result.motif_hits {
        let row = &mut rows[motif_hit.motif_index];
        row.total_hits += motif_hit.total_hits;
        row.forward_hits += motif_hit.forward_hits;
        row.revcomp_hits += motif_hit.revcomp_hits;
        if motif_hit.read_has_hit {
            row.reads_with_hit += 1;
        }
    }
}

fn scan_record(
    record: &Record,
    motifs: &[CompiledMotif],
    options: ScanOptions,
    emit_read_hits: bool,
) -> RecordResult {
    let read_masks = options.use_iupac.then(|| {
        record
            .seq
            .iter()
            .copied()
            .map(canonical_base_mask)
            .collect::<Vec<_>>()
    });
    let mut motif_hits = Vec::with_capacity(motifs.len());
    let mut read_hits = Vec::new();

    for (motif_index, motif) in motifs.iter().enumerate() {
        let forward_scan = scan_pattern(
            record,
            read_masks.as_deref(),
            &motif.forward,
            options,
            emit_read_hits,
        );
        if forward_scan.hit_count > 0 {
            if emit_read_hits {
                append_read_hits(
                    &mut read_hits,
                    record,
                    motif,
                    Strand::Forward,
                    &forward_scan.positions,
                    motif.len(),
                );
            }
        }

        let reverse_scan = motif
            .reverse
            .as_ref()
            .map(|pattern| {
                scan_pattern(
                    record,
                    read_masks.as_deref(),
                    pattern,
                    options,
                    emit_read_hits,
                )
            })
            .unwrap_or_default();
        if reverse_scan.hit_count > 0 {
            if emit_read_hits {
                append_read_hits(
                    &mut read_hits,
                    record,
                    motif,
                    Strand::Reverse,
                    &reverse_scan.positions,
                    motif.len(),
                );
            }
        }

        let forward_hits = forward_scan.hit_count;
        let revcomp_hits = reverse_scan.hit_count;
        let total_hits = forward_hits + revcomp_hits;
        let read_has_hit = total_hits > 0;

        motif_hits.push(MotifHitSummary {
            motif_index,
            total_hits,
            forward_hits,
            revcomp_hits,
            read_has_hit,
        });
    }

    RecordResult {
        motif_hits,
        read_hits,
    }
}

fn append_read_hits(
    sink: &mut Vec<ReadHitRow>,
    record: &Record,
    motif: &CompiledMotif,
    strand: Strand,
    positions: &[usize],
    motif_length: usize,
) {
    for position in positions {
        let window = &record.seq[*position..(*position + motif_length)];
        sink.push(ReadHitRow {
            read_id: record.id.clone(),
            motif: motif.name.clone(),
            strand,
            position: *position,
            matched_sequence: String::from_utf8_lossy(window).into_owned(),
        });
    }
}

fn scan_pattern(
    record: &Record,
    read_masks: Option<&[u8]>,
    pattern: &crate::motif::Pattern,
    options: ScanOptions,
    collect_positions: bool,
) -> PatternScanResult {
    if pattern.sequence.len() > record.seq.len() {
        return PatternScanResult::default();
    }

    let mut result = PatternScanResult {
        hit_count: 0,
        positions: if collect_positions {
            Vec::with_capacity(4)
        } else {
            Vec::new()
        },
    };

    if options.use_iupac {
        if let Some(masks) = read_masks {
            let last_start = masks.len() - pattern.masks.len();
            'outer: for start in 0..=last_start {
                for (offset, motif_mask) in pattern.masks.iter().enumerate() {
                    if motif_mask & masks[start + offset] == 0 {
                        continue 'outer;
                    }
                }
                result.hit_count += 1;
                if collect_positions {
                    result.positions.push(start);
                }
            }
        }
    } else {
        for position in exact_positions_iter(&record.seq, &pattern.sequence) {
            result.hit_count += 1;
            if collect_positions {
                result.positions.push(position);
            }
        }
    }

    result
}

#[cfg(test)]
fn exact_positions(sequence: &[u8], pattern: &[u8]) -> Vec<usize> {
    exact_positions_iter(sequence, pattern).collect()
}

fn exact_positions_iter<'a>(
    sequence: &'a [u8],
    pattern: &'a [u8],
) -> impl Iterator<Item = usize> + 'a {
    let pattern_len = pattern.len();
    let second_base = pattern.get(1).copied();
    let last_base = pattern.last().copied().unwrap_or(pattern[0]);

    memchr_iter(pattern[0], sequence).filter(move |&position| {
        if position + pattern_len > sequence.len() {
            return false;
        }
        if let Some(second_base) = second_base {
            if sequence[position + 1] != second_base {
                return false;
            }
        }
        if sequence[position + pattern_len - 1] != last_base {
            return false;
        }
        exact_match_window(&sequence[position..position + pattern_len], pattern)
    })
}

#[cfg(test)]
fn iupac_positions(read_masks: &[u8], motif_masks: &[u8]) -> Vec<usize> {
    if motif_masks.len() > read_masks.len() {
        return Vec::new();
    }

    let last_start = read_masks.len() - motif_masks.len();
    let mut hits = Vec::new();
    'outer: for start in 0..=last_start {
        for (offset, motif_mask) in motif_masks.iter().enumerate() {
            if motif_mask & read_masks[start + offset] == 0 {
                continue 'outer;
            }
        }
        hits.push(start);
    }
    hits
}

#[inline]
fn exact_match_window(window: &[u8], pattern: &[u8]) -> bool {
    if window.len() != pattern.len() {
        return false;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if window.len() >= 32 && std::is_x86_feature_detected!("avx2") {
            unsafe {
                return avx2_equal(window, pattern);
            }
        }
        if window.len() >= 16 && std::is_x86_feature_detected!("sse2") {
            unsafe {
                return sse2_equal(window, pattern);
            }
        }
    }

    window == pattern
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn avx2_equal(window: &[u8], pattern: &[u8]) -> bool {
    let mut offset = 0;
    while offset + 32 <= window.len() {
        let lhs = _mm256_loadu_si256(window.as_ptr().add(offset) as *const __m256i);
        let rhs = _mm256_loadu_si256(pattern.as_ptr().add(offset) as *const __m256i);
        let cmp = _mm256_cmpeq_epi8(lhs, rhs);
        if _mm256_movemask_epi8(cmp) != -1 {
            return false;
        }
        offset += 32;
    }
    window[offset..] == pattern[offset..]
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn sse2_equal(window: &[u8], pattern: &[u8]) -> bool {
    let mut offset = 0;
    while offset + 16 <= window.len() {
        let lhs = _mm_loadu_si128(window.as_ptr().add(offset) as *const __m128i);
        let rhs = _mm_loadu_si128(pattern.as_ptr().add(offset) as *const __m128i);
        let cmp = _mm_cmpeq_epi8(lhs, rhs);
        if _mm_movemask_epi8(cmp) != 0xFFFF {
            return false;
        }
        offset += 16;
    }
    window[offset..] == pattern[offset..]
}

#[cfg(test)]
mod tests {
    use crate::io::{Record, SourceFormat};
    use crate::motif::{compile_motifs, RawMotif};

    use super::{exact_match_window, exact_positions, iupac_positions, scan_record, ScanOptions};

    fn demo_record(id: &str, seq: &str, qual: Option<Vec<u8>>) -> Record {
        Record {
            id: id.to_string(),
            seq: seq.as_bytes().to_vec(),
            qual,
            source_format: SourceFormat::Fastq,
        }
    }

    #[test]
    fn exact_matching_finds_overlapping_hits() {
        let hits = exact_positions(b"AAAAA", b"AAA");
        assert_eq!(hits, vec![0, 1, 2]);
    }

    #[test]
    fn iupac_matching_accepts_compatible_bases() {
        let read_masks = b"ATGAA"
            .iter()
            .map(|base| crate::motif::iupac_mask(*base))
            .collect::<Vec<_>>();
        let motif_masks = b"ATGRN"
            .iter()
            .map(|base| crate::motif::iupac_mask(*base))
            .collect::<Vec<_>>();
        let hits = iupac_positions(&read_masks, &motif_masks);
        assert_eq!(hits, vec![0]);
    }

    #[test]
    fn reads_with_hit_and_total_hits_are_distinct() {
        let motifs = compile_motifs(
            &[RawMotif {
                name: "m1".to_string(),
                sequence: "AAA".to_string(),
            }],
            false,
            false,
        )
        .unwrap();
        let result = scan_record(
            &demo_record("r1", "AAAAA", Some(vec![40; 5])),
            &motifs,
            ScanOptions {
                use_iupac: false,
            },
            false,
        );
        assert_eq!(result.motif_hits[0].total_hits, 3);
        assert!(result.motif_hits[0].read_has_hit);
    }

    #[test]
    fn reverse_complement_hits_are_detected() {
        let motifs = compile_motifs(
            &[RawMotif {
                name: "m1".to_string(),
                sequence: "ATTATGAGAATAGTGTG".to_string(),
            }],
            true,
            false,
        )
        .unwrap();
        let reverse = "CACACTATTCTCATAAT";
        let result = scan_record(
            &demo_record("r1", reverse, Some(vec![40; reverse.len()])),
            &motifs,
            ScanOptions {
                use_iupac: false,
            },
            true,
        );
        assert_eq!(result.motif_hits[0].revcomp_hits, 1);
    }

    #[test]
    fn iupac_mode_does_not_treat_ambiguous_read_bases_as_matches() {
        let motifs = compile_motifs(
            &[RawMotif {
                name: "iupac".to_string(),
                sequence: "ATGRN".to_string(),
            }],
            false,
            true,
        )
        .unwrap();
        let result = scan_record(
            &demo_record("r1", "ATGNN", Some(vec![40; 5])),
            &motifs,
            ScanOptions {
                use_iupac: true,
            },
            false,
        );
        assert_eq!(result.motif_hits[0].total_hits, 0);
    }

    #[test]
    fn simd_window_match_falls_back_correctly() {
        let pattern = b"ATTATGAGAATAGTGTGATTATGAGAATAGTGTG";
        assert!(exact_match_window(pattern, pattern));
        assert!(!exact_match_window(
            pattern,
            b"ATTATGAGAATAGTGTGATTATGAGAATAGTGTA"
        ));
    }
}
