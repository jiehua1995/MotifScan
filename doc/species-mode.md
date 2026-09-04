# Species-aware fuzzy mode

`motifscan species` is designed for long-read DNA/cDNA data where a long diagnostic reference window is useful for locus specificity, but exact matching is too brittle because reads contain sequencing errors.

The mode separates two questions:

1. **Does this read come from this homologous locus?** A set of exact shared k-mer anchors finds candidates quickly, then a semi-global edit alignment of the full long reference window tolerates substitutions, insertions, and deletions. Diagnostic SNP positions are excluded from the shared-site identity calculation.
2. **Does the read support the mel or sim allele?** Diagnostic substitution SNPs are extracted automatically by globally aligning the paired references. Only covered, sufficiently high-quality SNP bases vote for species identity. Pair-specific indels are reported for QC but do not vote by default.

This is intended for data such as Oxford Nanopore full-length cDNA where 100-150 bp diagnostic windows contain multiple species-informative SNPs.

## Input files

### Motif/reference CSV

Same two-column format as exact mode:

```csv
name,sequence
18S_dmel,AATTCCGATAACGAAC...
18S_dsim,AATTCCGATAACGAAC...
```

### Pair CSV

Three columns identify the homologous locus and which reference is treated as mel/sim:

```csv
locus,mel,sim
18S,18S_dmel,18S_dsim
28S,28S_dmel,28S_dsim
FBgn0031044__LOC27207129,FBgn0031044,LOC27207129
```

The two references may have different lengths. MotifScan globally aligns them before extracting diagnostic substitution SNPs, so coordinates do not drift after an indel.

### Multiple FASTQ files

`--input-list` accepts either one path per line:

```text
/path/Ras3.fastq.gz
/path/ML82.fastq.gz
```

or an explicit sample name followed by a TAB and the path:

```text
Ras3	/path/Ras3.fastq.gz
ML82	/path/ML82.fastq.gz
HYC1_rep1	/path/HYC1_rep1.fastq.gz
```

## Example

```bash
motifscan species \
  --input-list rRNA_fastq_files.txt \
  --motifs motifs.csv \
  --pairs pairs.csv \
  -o rRNA_species_scan.csv \
  --pair-qc-output rRNA_species_scan.pairs.csv \
  -t 48 \
  --progress
```

Recommended initial defaults for Nanopore data are already the CLI defaults:

```text
--anchor-k 11
--anchors-per-locus 8
--alignment-slack 20
--min-shared-identity 0.85
--min-aligned-bases 80
--min-snp-baseq 15
--min-informative-snps 2
--species-fraction 0.75
--locus-mode best
```

For calibration, run pure parental samples first. A useful sensitivity analysis is to repeat with `--min-shared-identity 0.90` and `0.95` and verify that the mel/sim fraction is stable while read yield decreases as expected.

## Algorithm

For each locus pair, MotifScan:

1. globally aligns the mel and sim reference windows;
2. records substitution SNPs and pairwise indel columns;
3. chooses shared k-mers that do not overlap diagnostic SNP positions;
4. builds one Aho-Corasick index containing forward and reverse-complement anchors for all loci;
5. scans each read once for those anchors;
6. uses anchor coordinates to create a small candidate window;
7. runs a semi-global edit alignment of the full mel reference against that window;
8. calculates identity only over non-diagnostic positions, allowing ordinary sequencing errors;
9. maps diagnostic SNP coordinates through the alignment;
10. filters SNP observations by FASTQ Phred quality;
11. assigns `mel`, `sim`, `conflict`, `ambiguous`, or `no_snp`;
12. by default keeps only the best locus assignment for each read.

The exact anchors are **candidate retrieval only**. A read is not called a locus/species hit merely because an anchor is present.

## Output

The main CSV contains one row per `sample x locus`, including:

```text
sample,input_file,locus,mel_motif,sim_motif,mel_sequence,sim_sequence,
mel_length,sim_length,diagnostic_snp_count,pairwise_indel_columns,diagnostic_snps,
reads_with_hit,total_hits,forward_hits,revcomp_hits,
mel_reads,sim_reads,ambiguous_reads,conflict_reads,no_snp_reads,
species_assigned_reads,mel_fraction_among_assigned,sim_fraction_among_assigned,
mean_shared_identity,mean_informative_snps_per_hit
```

A second `*.snps.csv` reports each diagnostic SNP separately, including mel/sim/other counts, low-quality observations, depth, allele fractions, and mean high-quality base quality. This is important QC: independent SNPs within one reference window should normally support similar allele fractions.

`--pair-qc-output` writes the automatically extracted diagnostic SNP positions. Inspect this file before interpreting biological results, especially for reference pairs with substantial length differences.

## Interpretation and QC

- `reads_with_hit` means the long-window approximate alignment passed the locus filters; it is not an exact 150-bp hit count.
- `no_snp_reads` are valid locus hits that did not cover enough high-quality diagnostic SNPs to assign species.
- `conflict_reads` contain support for both alleles without either reaching the required species fraction.
- `ambiguous_reads` have enough informative observations but are dominated by errors/other bases rather than a clean mel/sim vote.
- Species-specific indels are intentionally not used as votes in the initial implementation because Nanopore indel errors are harder to model robustly than substitutions.
- For rDNA/rRNA, interpret within-locus mel/sim allele fractions more strongly than absolute read counts between different rRNA regions, because library preparation and reverse-transcription efficiency can differ by region.

## Validation strategy

Unit tests include synthetic mel/sim references and reads covering:

- exact mel and exact sim reads;
- shared-position substitutions;
- shared-position indels;
- reverse-complement reads;
- automatic diagnostic SNP extraction.

The CI workflow runs formatting, clippy, unit tests, and release builds. For a new experimental dataset, pure parental samples remain the most important empirical calibration: mel parents should overwhelmingly call mel and sim parents should overwhelmingly call sim.
