# MotifScan Design

## Goals

MotifScan is designed as a streaming, low-memory, multi-threaded Rust CLI for motif retrieval in sequencing reads. The implementation priorities are:

1. Correct exact and reverse-complement motif detection.
2. Support for FASTA, FASTQ, and gzip-compressed input without whole-file loading.
3. Scalable processing for short and long reads.
4. Stable comma-separated tabular outputs for downstream analysis.
5. Clear code structure so developers can extend matching, I/O, and reporting paths safely.

## High-Level Architecture

The project is split into five main modules.

### `src/cli.rs`

Defines the command-line surface using Clap.

Responsibilities:

- Parse subcommands and arguments.
- Validate mutual exclusions and required argument combinations.
- Provide the explicit `--progress` flag and custom `-v/--version` behavior.

### `src/io.rs`

Implements streaming input readers.

Responsibilities:

- Open uncompressed and gzip-compressed files.
- Detect FASTA vs FASTQ.
- Parse multiline FASTA.
- Parse standard 4-line FASTQ.
- Normalize sequence bytes to uppercase.
- Decode FASTQ qualities as Phred+33.
- Expose chunked `next_chunk()` iteration for parallel scanning.
- Track byte-level input progress for optional progress reporting.

### `src/motif.rs`

Owns motif parsing and compilation.

Responsibilities:

- Parse motif tables from two-column CSV input.
- Validate single-motif and motif-table inputs.
- Normalize motif sequence case.
- Precompute forward motif representation.
- Precompute reverse-complement motif representation when requested.
- Detect palindromic motifs to avoid double counting.

### `src/scanner.rs`

Contains the scan engine and result aggregation logic.

Responsibilities:

- Coordinate count flows.
- Read input in chunks and process records in parallel using Rayon.
- Run the exact-matching scan path.
- Aggregate motif hit counts.
- Emit optional progress information.

### `src/output.rs`

Owns stable table writing.

Responsibilities:

- Create writers for summary and report files.
- Emit fixed-column CSV tables.
- Serialize read-level hit reports.

## Execution Flow

### Count Mode

1. Parse CLI arguments.
2. Load motif definitions from `--motif` or `--motifs`.
3. Compile motifs and optional reverse complements.
4. Open a streaming reader for the input file.
5. Repeatedly read chunks of records.
6. Process each chunk in parallel.
7. Merge per-record motif counts into global rows.
8. Write final CSV summary and optional read-level hit output.

## Matching Design

### Exact Matching Path

The exact path is optimized for the common case of mismatch-free motif retrieval.

Algorithm:

1. Use `memchr` to find candidate positions matching the first motif byte.
2. Reject candidates early if the second byte or last byte does not match.
3. Compare the remaining window using SIMD-assisted equality when the CPU supports SSE2 or AVX2.
4. Fall back to scalar slice comparison when SIMD is unavailable or the window is short.

This path keeps branching low and avoids full-window comparison for most non-matching candidates.

### Reverse Complement Strategy

Reverse-complement motifs are precomputed during motif compilation rather than at scan time.

Benefits:

- No repeated reverse-complement construction during scanning.
- Symmetric handling between forward and reverse paths.
- Palindromic motifs can be recognized once and skipped on the reverse path.

## Streaming and Memory Model

The implementation is intentionally streaming-first.

Key decisions:

- Input is consumed chunk by chunk instead of loading all records into memory.
- Each `Record` owns its sequence and optional quality bytes.
- Chunk size is bounded so memory use scales with chunk size rather than file size.
- Gzip input is processed via a buffered decoder instead of full decompression to disk.

Progress reporting is derived from streamed input bytes rather than from a separate counting pass, which keeps startup latency low.

## Parallelism Model

MotifScan uses chunk-level parallelism.

Flow:

1. The reader loads a chunk of records.
2. Rayon processes records in that chunk with `par_iter` or `into_par_iter`.
3. Each record produces local hit summaries.
4. The chunk is reduced into global output rows.

Why this model:

- Avoids a global lock for each record.
- Keeps synchronization overhead moderate.
- Works for both short-read and long-read inputs.
- Preserves streaming behavior without a complex producer-consumer queue.

## Output Design

All outputs are written as comma-separated tables.

Rationale:

- One stable output format is easier for downstream pipelines.
- Users can choose their own file extension.
- Writers and documentation stay simpler than supporting multiple delimiters.

Current outputs:

- Count summary
- Read-level hit report

## Progress Reporting

Progress is opt-in with `--progress`.

The progress bar reports:

- Current mode and input file name
- Motif count
- Processed reads
- Average read length
- Reads per second
- Byte progress and ETA based on the streamed input source

For gzip input, ETA is based on compressed byte progress and should be treated as an approximation.

## Testing Strategy

The project currently uses three layers of testing.

### Unit Tests

Located inside the source modules.

Coverage includes:

- Reverse complement correctness
- IUPAC masks
- Exact match overlap behavior
- Reverse-complement hit detection
- FASTA parsing
- FASTQ parsing
- Gzip smoke parsing

### Large Integration Tests

Located in `tests/cli_large_end_to_end.rs`.

Coverage includes:

- End-to-end CLI count execution on a larger synthetic dataset
- Reverse-complement handling
- Gzip input behavior
- Progress-enabled execution path
- Exact output row values for key summary fields

### Manual Benchmark Validation

The `benchmark/` directory stores demo command outputs and timing records. These are useful for regression checks and user-facing examples.

## Important Data Structures

### `Record`

Represents one FASTA or FASTQ entry.

Fields:

- `id`
- `seq`
- `qual`
- `source_format`

### `CompiledMotif`

Represents one precompiled motif.

Fields include:

- motif name
- uppercase forward sequence
- optional reverse-complement pattern
- palindrome flag

## Known Constraints

- FASTQ support is currently limited to standard 4-line layout.
- Only exact motif matching is supported.
- There is no Aho-Corasick multi-motif exact fast path yet.
- Approximate matching and mismatch-tolerant matching are not implemented.
- ETA for gzip input is approximate because compressed bytes do not scale linearly with logical records.

## Extension Points

Good next development targets:

1. Add `--no-progress` or richer progress sinks for batch environments.
2. Add Aho-Corasick for many exact motifs.
3. Add more formal benchmark harnesses.
4. Add approximate matching with bounded mismatches.
5. Add richer machine-readable run metadata.

## Developer Notes

When extending the codebase:

- Keep the stream-first model unless there is a measured reason not to.
- Prefer precomputation in `motif.rs` over repeated work in the hot scan loop.
- Keep output schemas stable.
- Validate changes with both unit tests and large CLI integration tests.
- Benchmark exact-path changes separately so regressions are attributable.
