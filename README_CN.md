# MotifScan 中文说明

<div align="center">
  <a href="https://github.com/jiehua1995/MotifScan">
    <img src="doc/logo.png" alt="MotifScan 标志" width="420" style="border-radius:10px; box-shadow:0 8px 24px rgba(0,0,0,0.14);">
  </a>
  <p style="margin:8px 0 0 0; font-size:16px; color:#444;">流式、低内存的 motif 扫描 CLI（Rust）— 支持 FASTA/FASTQ、CSV motif、Aho–Corasick 加速</p>
</div>

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
motifscan -v
motifscan --version

```bibtex
## 下载发布版（预编译二进制）与校验

如果你不想自行从源码打包编译，可以直接从 GitHub Releases 下载预编译的发布包。发布包一般以版本号与平台命名，例如：

- `motifscan-0.1.2-linux.tar.gz`
- `motifscan-0.1.2-macos.tar.gz`
- `motifscan-0.1.2-windows.zip`

每个发布包通常会附带一个 SHA256 校验文件（扩展名为 `.sha256`）。下载后请务必校验文件完整性：

Linux / macOS：

```bash
sha256sum -c motifscan-0.1.2-linux.sha256
# 或
shasum -a 256 -c motifscan-0.1.2-linux.sha256
```

Windows (PowerShell)：

```powershell
Get-FileHash -Algorithm SHA256 .\motifscan-0.1.2-windows.zip
# 比对 motifscan-0.1.2-windows.sha256 中的值
```

如果你需要从源码自行打包发布，请参考仓库 `doc/release_cn.md` 中的中文步骤说明。
 
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

## 发布产物与校验

CI 生成的 release 包会使用 `Cargo.toml` 中的语义版本号并包含平台标识，例如：

- `motifscan-0.1.2-linux.tar.gz`
- `motifscan-0.1.2-macos.tar.gz`
- `motifscan-0.1.2-windows.zip`

每个包会同时生成一个 SHA256 校验文件（`.sha256`），下载后建议校验：

Linux / macOS：

```bash
sha256sum -c motifscan-0.1.2-linux.sha256
# 或
shasum -a 256 -c motifscan-0.1.2-linux.sha256
```

Windows (PowerShell)：

```powershell
Get-FileHash -Algorithm SHA256 .\motifscan-0.1.2-windows.zip
# 比对 motifscan-0.1.2-windows.sha256 中的值
```

