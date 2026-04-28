# MotifScan 中文说明

MotifScan 是一个用 Rust 编写的流式、低内存、多线程 motif 扫描命令行工具，适用于 FASTA 和 FASTQ reads。

- 只支持 exact matching
- 可选反向互补搜索
- motif 输入和输出都使用 CSV
- 支持 FASTA、FASTQ、FASTA.GZ 和 FASTQ.GZ

## 安装

```bash
cargo build --release
```

可执行文件路径：

```bash
./target/release/motifscan
```

查看版本：

```bash
motifscan -v
motifscan --version
```

## 引用

```bibtex
@software{motifscan,
  author = {jiehua1995},
  title = {MotifScan},
  url = {https://github.com/jiehua1995/MotifScan},
  version = {0.1.0}
}
```

## 快速使用

扫描多个 motif：

```bash
motifscan count \
  -i reads.fastq \
  --motifs motifs.csv \
  --revcomp \
  -o count.csv
```

扫描单个 motif：

```bash
motifscan count \
  -i reads.fa \
  --motif ATTATGAGAATAGTGTG \
  --motif-name motif1 \
  -o count.csv
```

输出 read-level hits：

```bash
motifscan count \
  -i reads.fastq \
  --motifs motifs.csv \
  --report-read-hits read_hits.csv \
  -o count.csv
```

## 主要参数

- `-i`, `--input <FILE>`：输入 reads 文件
- `--motif <SEQUENCE>`：直接在命令行指定一个 motif
- `--motif-name <NAME>`：`--motif` 对应的名称，默认是 `motif`
- `--motifs <FILE>`：两列 CSV motif 表格
- `--revcomp`：同时搜索反向互补序列
- `-t`, `--threads <INT>`：线程数
- `--progress`：在 stderr 显示进度
- `-o`, `--output <FILE>`：汇总 CSV 输出
- `--report-read-hits <FILE>`：可选的 hit 明细 CSV 输出

## motif CSV 格式

```text
name,sequence
motif1,ATTATGAGAATAGTGTG
motif2,TTCATTCATGGTGGCAGTAAAATGTTTATTGTG
motif3,ATGAA
```

规则：

- 只支持逗号分隔
- 表头可选
- 必须严格是两列：`name,sequence`
- motif 只允许使用精确碱基：`A`、`C`、`G`、`T`、`U`

## 输出 CSV 列

汇总输出：

```text
motif,sequence,length,reads_with_hit,total_hits,forward_hits,revcomp_hits
```

read-level hit 输出：

```text
read_id,motif,strand,position,matched_sequence
```

## 说明

- 输入序列会先统一转成大写再匹配。
- 支持 overlapping hits。
- 回文 motif 在反向互补模式下不会重复计数。
- 如果 motif 比 read 更长，会跳过该 read。
- FASTQ 当前按标准 4 行格式处理。

