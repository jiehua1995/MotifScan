# MotifScan

MotifScan is a streaming, low-memory, multi-threaded Rust CLI for motif scanning in FASTA and FASTQ sequencing reads. It supports exact matching with a memchr fast path, IUPAC-compatible matching, reverse-complement search, gzip-compressed input, and read-level hit reports.

The first implementation targets exact motif search across Sanger, Illumina, and Nanopore/PacBio style reads while keeping I/O streaming by default.

## Installation

Build the release binary from the repository root:

```bash
cargo build --release
```

Print version and citation text:

```bash
motifscan -v
motifscan --version
```

The executable will be available at:

```bash
./target/release/motifscan
```

On Windows shells, use `motifscan.exe` when needed.

## Supported Input

- FASTA
- FASTQ
- FASTA.GZ
- FASTQ.GZ

Notes:

- Parsing is streaming by default and does not read the whole input file into memory.
- Multiline FASTA is supported.
- FASTQ currently targets the standard 4-line format.
- Gzip input is decoded with a buffered streaming `flate2` pipeline; gzip files are not memory-mapped.
- Exact matching uses a memchr candidate scan plus SIMD-accelerated byte comparison when the CPU supports SSE2 or AVX2.

## Count Mode

Count motif hits across reads:

```bash
motifscan count \
	-i reads.fastq \
	--motifs motifs.tsv \
	--revcomp \
	-t 4 \
	-o count.csv
```

Scan a single motif:

```bash
motifscan count \
	-i reads.fa \
	--motif ATTATGAGAATAGTGTG \
	--motif-name Dmel_28 \
	--revcomp \
	-t 2 \
	-o count.csv
```

Write read-level hits:

```bash
motifscan count \
	-i reads.fastq \
	--motifs motifs.tsv \
	--revcomp \
	--report-read-hits read_hits.csv \
	-o count.csv
```

## CLI Summary

```text
motifscan count \
	-i <input.fa/fq/fa.gz/fq.gz> \
	--motif <SEQUENCE> \
	--motif-name <NAME> \
	--motifs <motifs.tsv/csv/txt> \
	--revcomp \
	--iupac \
	--threads <INT> \
	--output <FILE> \
	--report-read-hits <FILE>
```

## Command And Argument Tables

### Subcommands

| Command | Usage | Purpose | Notes |
| --- | --- | --- | --- |
| `motifscan count` | Scan one motif or a motif list against reads and report hit counts. | Use when you want per-motif totals, per-read hit counts, or filtered count summaries. | Supports `--motif` or `--motifs`. |

### Global Options

| Option | Usage | Required | Explanation |
| --- | --- | --- | --- |
| `-v`, `--version` | `motifscan -v` | No | Print the current MotifScan version together with citation guidance. |

### `count` Arguments

| Argument | Usage | Required | Explanation |
| --- | --- | --- | --- |
| `-i`, `--input <FILE>` | `-i reads.fastq` | Yes | Input read file. Supports FASTA, FASTQ, FASTA.GZ, and FASTQ.GZ. |
| `--motif <SEQUENCE>` | `--motif ATTATGAGAATAGTGTG` | Conditionally required | Scan one motif sequence directly from the command line. Mutually exclusive with `--motifs`. |
| `--motif-name <NAME>` | `--motif-name Dmel_28` | No | Name used in the output when `--motif` is supplied. Defaults to `motif`. |
| `--motifs <FILE>` | `--motifs motifs.tsv` | Conditionally required | Path to a two-column motif table. Use this instead of `--motif` when scanning multiple motifs. |
| `--revcomp` | `--revcomp` | No | Also scan the reverse complement of each motif. Palindromic motifs are not double-counted. |
| `--iupac` | `--iupac` | No | Enable IUPAC matching for ambiguous bases in motifs only. If a motif contains IUPAC bases and this flag is not set, MotifScan returns an error. |
| `-t`, `--threads <INT>` | `-t 4` | No | Number of worker threads. Defaults to the detected CPU count. |
| `--progress` | `--progress` | No | Enable a live stderr progress bar with input file name, motif count, processed reads, average read length, reads per second, byte progress, and ETA. Disabled by default. |
| `-o`, `--output <FILE>` | `-o count.csv` | Yes | Output path for the count summary table. Output is always comma-separated text regardless of file extension. |
| `--report-read-hits <FILE>` | `--report-read-hits read_hits.csv` | No | Optional read-level hit report containing one row per hit. |

### Common Usage Patterns

| Goal | Command Pattern | Explanation |
| --- | --- | --- |
| Scan one known motif | `motifscan count -i <input> --motif <seq> --motif-name <name> -o <out>` | Best for validating one motif quickly without preparing a motif table. |
| Scan many motifs together | `motifscan count -i <input> --motifs <table> -o <out>` | Reads motif definitions from a file and writes one summary row per motif. |
| Search both strands | Add `--revcomp` | Required when the motif can appear on the reverse-complement strand. |
| Use ambiguous motif matching | Add `--iupac` | Required when motif tables contain symbols such as `R`, `Y`, or `N`. IUPAC applies to motifs only, not to ambiguous read bases. |
| Produce read-level hit locations | Add `--report-read-hits <file>` to `count` | Writes one line per hit with strand, position, and matched sequence. |

## Motif File Formats

Two-column motif table:

```text
name	sequence
motif1	ATTATGAGAATAGTGTG
motif2	TTCATTCATGGTGGCAGTAAAATGTTTATTGTG
iupac_test	ATGRN
```

TSV, CSV, and whitespace-delimited text are accepted.

## Matching Behavior

- Input sequences are normalized to uppercase.
- Exact mode uses a memchr-based first-base candidate scan followed by byte comparison.
- IUPAC mode converts motif bases to bitmasks and compares them against canonical read bases only.
- If a motif contains IUPAC bases and `--iupac` is not enabled, MotifScan fails fast with an error.
- Reverse-complement search can be enabled with `--revcomp`.
- Palindromic motifs are not double-counted on the reverse-complement path.
- Overlapping hits are counted.
- If a motif is longer than a read, that read is skipped for that motif.

## Reverse Complement

`--revcomp` scans both the forward motif and its reverse complement.

Supported reverse-complement mapping includes standard and IUPAC bases:

- `A <-> T`
- `C <-> G`
- `R <-> Y`
- `S <-> S`
- `W <-> W`
- `K <-> M`
- `B <-> V`
- `D <-> H`
- `N <-> N`

## IUPAC Matching

Enable IUPAC-compatible motif matching with `--iupac`.

Examples:

- `R` matches `A` or `G`
- `Y` matches `C` or `T`
- `N` matches any canonical base

In exact mode, characters are matched literally. In IUPAC mode, ambiguous symbols are allowed in motifs only. Read bases are treated as canonical `A/C/G/T` only, and ambiguous read bases do not satisfy IUPAC motif positions.

## Output Columns

Count summary:

```text
motif	sequence	length	reads_with_hit	total_hits	forward_hits	revcomp_hits
```

Read-level hit report:

```text
read_id	motif	strand	position	matched_sequence
```

Definitions:

- `reads_with_hit`: number of reads with at least one hit for that motif.
- `total_hits`: total number of matching windows, including multiple hits in one read.

## Multi-threading

Use `--threads` or `-t` to control worker count. The default is `num_cpus`.

The implementation uses chunked streaming reads and parallel chunk processing with Rayon. Aggregation happens after each chunk to avoid frequent global locking during per-read scanning.

## Limitations

- Gzip input is streaming-decoded and cannot use a true mmap fast path.
- Default operation does not search across different reads because reads are treated as separate molecules.
- FASTQ support currently targets the standard 4-line format.
- Approximate matching and mismatches are not implemented in this version.
- The current exact-matching fast path is per-motif; there is no multi-pattern Aho-Corasick path yet.

## Future Optimization Ideas

- SIMD exact-matching fast path
- Aho-Corasick fast path for large exact motif sets
- mmap parser for uncompressed FASTA/FASTQ
- Approximate matching and mismatch support