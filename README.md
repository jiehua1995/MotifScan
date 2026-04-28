# MotifScan

MotifScan is a streaming, low-memory, multi-threaded Rust CLI for motif scanning in FASTA and FASTQ reads.

- Exact matching only
- Optional reverse-complement scanning
- CSV motif input and CSV output
- FASTA, FASTQ, FASTA.GZ, and FASTQ.GZ support

## Installation

```bash
cargo build --release
```

Binary path:

```bash
./target/release/motifscan
```

Version:

```bash
motifscan -v
motifscan --version
```

## Citation

```bibtex
@software{motifscan,
  author = {jiehua1995},
  title = {MotifScan},
  url = {https://github.com/jiehua1995/MotifScan},
  version = {0.1.0}
}
```

## Quick Start

Multiple motifs:

```bash
motifscan count \
  -i reads.fastq \
  --motifs motifs.csv \
  --revcomp \
  -o count.csv
```

Single motif:

```bash
motifscan count \
  -i reads.fa \
  --motif ATTATGAGAATAGTGTG \
  --motif-name motif1 \
  -o count.csv
```

Read-level hits:

```bash
motifscan count \
  -i reads.fastq \
  --motifs motifs.csv \
  --report-read-hits read_hits.csv \
  -o count.csv
```

## Main Options

- `-i`, `--input <FILE>`: input reads file
- `--motif <SEQUENCE>`: one motif from the command line
- `--motif-name <NAME>`: name for `--motif`, default is `motif`
- `--motifs <FILE>`: two-column CSV motif table
- `--revcomp`: also scan reverse complement
- `-t`, `--threads <INT>`: worker threads
- `--progress`: show progress on stderr
- `-o`, `--output <FILE>`: summary CSV output
- `--report-read-hits <FILE>`: optional hit-level CSV output

## Motif CSV Format

```text
name,sequence
motif1,ATTATGAGAATAGTGTG
motif2,TTCATTCATGGTGGCAGTAAAATGTTTATTGTG
motif3,ATGAA
```

Rules:

- Comma-separated only
- Optional header row
- Exactly two columns: `name,sequence`
- Motifs must use exact bases only: `A`, `C`, `G`, `T`, `U`

## Output CSV Columns

Summary:

```text
motif,sequence,length,reads_with_hit,total_hits,forward_hits,revcomp_hits
```

Read hits:

```text
read_id,motif,strand,position,matched_sequence
```

## Notes

- Input is normalized to uppercase before matching.
- Overlapping hits are counted.
- Palindromic motifs are not double-counted in reverse-complement mode.
- If a motif is longer than a read, it is skipped for that read.
- FASTQ currently expects the standard 4-line format.