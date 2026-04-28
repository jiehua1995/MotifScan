# MotifScan 中文说明

MotifScan 是一个用 Rust 编写的流式、低内存、多线程 motif 扫描命令行工具，适用于 FASTA、FASTQ、FASTA.GZ 和 FASTQ.GZ 输入。它支持精确匹配、IUPAC 简并碱基匹配、反向互补搜索，以及 read-level hit 报告。

## 安装

在仓库根目录执行：

```bash
cargo build --release
```

查看版本和引用信息：

```bash
motifscan -v
motifscan --version
```

生成的可执行文件位于：

```bash
./target/release/motifscan
```

Windows 下可执行文件通常为 `motifscan.exe`。

## 支持的输入格式

- FASTA
- FASTQ
- FASTA.GZ
- FASTQ.GZ

说明：

- 默认使用流式解析，不会一次性把整个输入文件读入内存。
- 支持 multiline FASTA。
- FASTQ 当前优先支持标准 4 行格式。
- gzip 输入使用带缓冲的流式解压。
- exact matching 路径会使用 memchr 候选扫描，并在 CPU 支持时启用 SSE2 或 AVX2 SIMD 比较。
- 运行时会在 stderr 显示进度条，实时展示已处理 reads、处理速率和 ETA。

## 子命令总览

| 子命令 | 用法 | 作用 | 说明 |
| --- | --- | --- | --- |
| `motifscan count` | 统计一个或多个 motif 在 reads 中的命中情况 | 输出每个 motif 的 reads_with_hit、total_hits、forward_hits、revcomp_hits 等指标 | 可选输出 read-level hit 明细 |

## 全局参数

| 参数 | 用法 | 是否必需 | 解释 |
| --- | --- | --- | --- |
| `-v`, `--version` | `motifscan -v` | 否 | 打印版本号和引用说明 |

## `count` 子命令参数表

| 参数 | 用法示例 | 是否必需 | 解释 |
| --- | --- | --- | --- |
| `-i`, `--input <FILE>` | `-i reads.fastq` | 是 | 输入文件路径，支持 FASTA、FASTQ、FASTA.GZ、FASTQ.GZ |
| `--motif <SEQUENCE>` | `--motif ATTATGAGAATAGTGTG` | 条件必需 | 直接在命令行指定一个 motif 序列，和 `--motifs` 二选一 |
| `--motif-name <NAME>` | `--motif-name Dmel_28` | 否 | 单 motif 模式下输出中使用的 motif 名称，默认值为 `motif` |
| `--motifs <FILE>` | `--motifs motifs.tsv` | 条件必需 | 指定两列 motif 表格文件，和 `--motif` 二选一 |
| `--revcomp` | `--revcomp` | 否 | 同时搜索 motif 的反向互补序列 |
| `--iupac` | `--iupac` | 否 | 只对 motif 启用 IUPAC 简并碱基匹配。如果 motif 中含有 IUPAC 字符但未开启该参数，程序会直接报错 |
| `-t`, `--threads <INT>` | `-t 4` | 否 | 线程数，默认使用 CPU 核心数 |
| `--progress` | `--progress` | 否 | 打开 stderr 实时进度条，显示输入文件名、motif 数量、已处理 reads、平均 read 长度、reads/s、字节进度和 ETA。默认关闭 |
| `-o`, `--output <FILE>` | `-o count.csv` | 是 | count 汇总输出文件。输出内容始终为逗号分隔文本，文件后缀由用户自行决定 |
| `--report-read-hits <FILE>` | `--report-read-hits read_hits.csv` | 否 | 额外输出 read-level hit 明细，每个 hit 一行 |

## 常见使用方式

| 目标 | 命令模板 | 说明 |
| --- | --- | --- |
| 扫描单个 motif | `motifscan count -i <input> --motif <seq> --motif-name <name> -o <out>` | 适合快速验证一个已知 motif |
| 扫描多个 motif | `motifscan count -i <input> --motifs <table> -o <out>` | 从文件读取多个 motif，一次输出所有 motif 的统计结果 |
| 搜索双链 | 在命令中增加 `--revcomp` | 当 motif 可能出现在反向互补链时应开启 |
| 启用 motif 简并碱基匹配 | 在命令中增加 `--iupac` | 当 motif 中存在 `R`、`Y`、`N` 等字符时应开启。IUPAC 只作用于 motif，不作用于 read |
| 输出每个 hit 的位置 | 在 `count` 中增加 `--report-read-hits <file>` | 输出 strand、position、matched_sequence 等明细 |

## 命令示例

### 1. 统计多个 motif

```bash
motifscan count \
  -i reads.fastq \
  --motifs motifs.tsv \
  --revcomp \
  -t 4 \
  -o count.csv
```

### 2. 统计单个 motif

```bash
motifscan count \
  -i reads.fa \
  --motif ATTATGAGAATAGTGTG \
  --motif-name Dmel_28 \
  --revcomp \
  -t 2 \
  -o count.csv
```

### 3. 输出 read-level hits

```bash
motifscan count \
  -i reads.fastq \
  --motifs motifs.tsv \
  --revcomp \
  --report-read-hits read_hits.csv \
  -o count.csv
```

## motif 表格格式

### motif 两列表格

```text
name	sequence
motif1	ATTATGAGAATAGTGTG
motif2	TTCATTCATGGTGGCAGTAAAATGTTTATTGTG
iupac_test	ATGRN
```

支持 TSV、CSV 和空白字符分隔文本。

## 输出列说明

### count 汇总输出

```text
motif	sequence	length	reads_with_hit	total_hits	forward_hits	revcomp_hits
```

### read-level hit 输出

```text
read_id	motif	strand	position	matched_sequence
```

字段解释：

- `reads_with_hit`：至少命中该 motif 一次的 read 数量
- `total_hits`：所有命中窗口的总数，同一个 read 内多次命中会累计
- `forward_hits`：正向 motif 命中次数
- `revcomp_hits`：反向互补 motif 命中次数

## 匹配和分类逻辑

- 输入序列统一转成大写后再匹配。
- exact 模式按字节精确匹配。
- `--iupac` 模式下，只允许 motif 使用简并碱基，read 仍按 canonical `A/C/G/T` 处理。
- 如果 motif 中含有 IUPAC 字符但没有开启 `--iupac`，程序会直接报错。
- `--revcomp` 开启后会同时扫描反向互补序列。
- 回文 motif 不会在 reverse complement 路径上重复计数。
- 支持 overlapping hits。
- motif 比 read 更长时会直接跳过。
## 多线程与流式处理

- 线程数通过 `--threads` 控制。
- 读取端按 chunk 流式读取 records。
- 每个 chunk 内部用 Rayon 并行处理 reads。
- 汇总在 chunk 完成后统一 reduce，避免每条 read 都锁全局结构。
- gzip 输入也走流式解压，不需要先完整解压到磁盘。

## 已知限制

- gzip 输入不能使用真正的 mmap fast path。
- 当前 FASTQ 只保证标准 4 行格式。
- SIMD 目前只加速 exact matching 路径。
- 还没有 Aho-Corasick 多 motif fast path。
- 还没有 mismatch 或 approximate matching 支持。

