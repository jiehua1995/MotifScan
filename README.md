# MotifScan

<div style="display:flex;align-items:center;gap:18px;">
  <div style="flex:0 0 auto;">
    <img src="doc/logo.png" alt="MotifScan logo" width="120" style="border-radius:8px; box-shadow:0 6px 18px rgba(0,0,0,0.12);">
  </div>
  <div style="flex:1 1 auto;">
    <h2 style="margin:0 0 6px 0; font-weight:700;">Streaming, low-memory motif scanning CLI in Rust</h2>
    <p style="margin:0;color:#444;">Supports FASTA/FASTQ (gzipped), CSV motif tables, optional reverse-complement scanning, and Aho–Corasick acceleration for many motifs.</p>
  </div>
</div>

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

## 下载发布版（预编译二进制）与校验

如果你不想自己从源码打包编译，可以直接从 GitHub Releases 下载预编译的发布包。发布包一般以语义版本号与平台命名，例如：

- `motifscan-0.1.2-linux.tar.gz`
- `motifscan-0.1.2-macos.tar.gz`
- `motifscan-0.1.2-windows.zip`

每个发布包通常会附带一个 SHA256 校验文件（扩展名为 `.sha256`）。下载后请务必校验文件完整性：

Linux / macOS:

```bash
sha256sum -c motifscan-0.1.2-linux.sha256
# 或者
shasum -a 256 -c motifscan-0.1.2-linux.sha256
```

Windows (PowerShell)：

```powershell
Get-FileHash -Algorithm SHA256 .\motifscan-0.1.2-windows.zip
# 将输出的哈希值与 motifscan-0.1.2-windows.sha256 中记录的值进行比对
```

如果你需要从源码自己打包（例如为了自定义编译选项或特定平台优化），请参见仓库下的 `doc/release_cn.md` 获取中文打包与发布步骤说明。