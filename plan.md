Status: approved for implementation on 2026-04-28.


## Plan: Build MotifScan Rust CLI

从零搭建一个流式、低内存、多线程的 Rust 命令行工具 motifscan，先保证功能闭环和可验证性，再在 exact matching 路径上落实 memchr fast path。实现顺序按 I/O -> motif 编译与匹配 -> 并行扫描与聚合 -> 输出与 demo/benchmark -> 测试与 README 展开，确保每一步都可单独验证。

**Steps**
1. Phase 1: 项目骨架与依赖
   1. 创建 Cargo.toml，配置 package 名 MotifScan、binary 名 motifscan，以及 anyhow、clap(derive)、flate2、memchr、rayon、csv、num_cpus 等依赖。
   2. 创建 src/main.rs、src/cli.rs、src/io.rs、src/motif.rs、src/scanner.rs、src/output.rs，先定义最小可编译接口与核心数据结构，避免后续跨模块返工。
   3. 在 cli.rs 中用 clap derive 定义 count/classify 子命令、参数互斥关系、默认线程数、输出格式枚举和参数校验入口。
2. Phase 2: 流式输入解析
   1. 在 io.rs 实现 Record 结构，至少包含 id、seq、qual、source_format 等字段；序列统一大写，quality 保留原始 phred+33 数值或字节。
   2. 实现按文件扩展名或内容分支的输入打开逻辑：未压缩文件走 BufReader，.gz 走 flate2 流式 decoder。
   3. 实现 multiline FASTA parser 和标准 4-line FASTQ parser；对空 read、格式错误给出 anyhow 错误；对 FASTA + 质量过滤参数记录 warning 而非失败。
   4. 提供按 chunk 读取 Vec<Record> 的 streaming API，供 scanner 以块级并行消费，避免整文件入内存。
3. Phase 3: motif 编译与匹配核心
   1. 在 motif.rs 定义 RawMotif、CompiledMotif、Strand 枚举，以及从 CLI 输入或 motifs.tsv/csv/txt 加载 motif 的解析逻辑。
   2. 实现 sequence uppercase、非法空 motif 拒绝、IUPAC 字符到 bitmask 的映射，以及 reads 中非法碱基在 IUPAC 模式下的一致处理策略。
   3. 实现 reverse_complement，覆盖标准碱基和 IUPAC 映射；为 CompiledMotif 预编译 forward/revcomp 版本并标记是否为 palindrome，避免双计数。
   4. 实现两条匹配路径：
      - exact path：非 IUPAC 时用 memchr 锁定首碱基候选，再逐字节比较，支持 overlapping hits。
      - iupac path：motif/read 逐位 bitmask 兼容匹配，条件为 motif_mask[i] & base_mask[j+i] != 0。
   5. 在质量过滤逻辑中实现 motif window 的 min-baseq 与 min-mean-baseq 检查，仅对 FASTQ 生效。
4. Phase 4: 扫描、聚合与分类
   1. 在 scanner.rs 设计 chunk 级处理：reader 按块读取记录，rayon 对块内 records 并行扫描，每个任务维护 local counts/local classification/local read-hit 缓冲，再 reduce。
   2. 实现 count 模式聚合，输出每个 motif 的 reads_with_hit、total_hits、forward_hits、revcomp_hits，并在需要时输出全部 read-level hits。
   3. 实现 classify 模式聚合：同一 read 在每个 group 内独立判断唯一命中、ambiguous、unclassified；不同 group 彼此独立。
   4. 明确 total_hits 与 reads_with_hit 的差异：同一 read 多次命中应增加 total_hits，但只增加一次 reads_with_hit。
   5. motif 长于 read 时直接跳过，不报错；空序列 read 不崩溃。
5. Phase 5: 输出层与稳定格式
   1. 在 output.rs 统一 TSV/CSV/TXT writer，固定列顺序，避免不同子命令输出漂移。
   2. count summary 列：motif、sequence、length、reads_with_hit、total_hits、forward_hits、revcomp_hits。
   3. read-hit 报告列：read_id、motif、strand、position、matched_sequence。
   4. classify 主输出列：group、motif、sequence、length、reads_with_hit、total_hits、forward_hits、revcomp_hits。
   5. classify summary 列：group、top_call、top_reads、second_call、second_reads、total_classified、top_ratio、ambiguous_reads、unclassified_reads。
6. Phase 6: demo 数据、单元测试与集成验证
   1. 创建 test/motifs.tsv、test/rdna_motifs.tsv、test/sanger.fa、test/illumina.fastq、test/nanopore.fastq、test/mixed.fa；数据要覆盖正向、反向互补、IUPAC、低质量过滤、长 reads、多 hits、ambiguous/unclassified。
   2. 生成 test/illumina.fastq.gz 作为 gzip smoke 输入。
   3. 编写单元/集成测试，至少覆盖：reverse complement、IUPAC bitmask、exact match、revcomp match、quality filtering、reads_with_hit vs total_hits、classify ambiguous、FASTA parser、FASTQ parser、gz smoke test。
   4. 在测试里优先复用 scanner/motif/io 的公开接口，而不是走 shell 命令，确保失败可定位。
7. Phase 7: README、benchmark 与最终验证
   1. 重写 README.md，覆盖用途、安装、count/classify 示例、motif 文件格式、FASTA/FASTQ/GZ 支持、revcomp、IUPAC、质量过滤、多线程、输出字段、D. melanogaster / D. simulans rDNA 区分示例、限制与后续优化。
   2. 创建 benchmark/README.md 说明 benchmark 采集方式；创建 benchmark/benchmark_results.tsv，记录 command、input、threads、elapsed_seconds、output。
   3. 依次运行 cargo fmt、cargo test、cargo build --release、demo count/classify 命令，并把产物写入 benchmark/。
   4. 若环境没有系统 gzip 命令，则改用 PowerShell 或 Rust/flate2 方式生成 illumina.fastq.gz，并在 benchmark/README.md 记录实际命令。

**Relevant files**
- c:\GitHub\MotifScan\Cargo.toml — 项目元数据、依赖、release profile。
- c:\GitHub\MotifScan\src\main.rs — 子命令分发、线程池配置、错误出口。
- c:\GitHub\MotifScan\src\cli.rs — clap 参数定义、互斥与默认值。
- c:\GitHub\MotifScan\src\io.rs — streaming FASTA/FASTQ/GZ 打开与 chunk reader。
- c:\GitHub\MotifScan\src\motif.rs — CompiledMotif、IUPAC bitmask、reverse complement、motif 文件加载。
- c:\GitHub\MotifScan\src\scanner.rs — exact/iupac 匹配、quality filter、并行扫描、count/classify reduce。
- c:\GitHub\MotifScan\src\output.rs — TSV/CSV/TXT writer 与固定列顺序输出。
- c:\GitHub\MotifScan\test\motifs.tsv — 两列 motif demo 与 IUPAC 样例。
- c:\GitHub\MotifScan\test\rdna_motifs.tsv — 三列 group 分类 demo。
- c:\GitHub\MotifScan\test\illumina.fastq — 短 reads + 质量过滤 + revcomp 覆盖。
- c:\GitHub\MotifScan\test\nanopore.fastq — 长 reads + multi-hit + ambiguous 覆盖。
- c:\GitHub\MotifScan\test\sanger.fa — FASTA 正向/反向互补覆盖。
- c:\GitHub\MotifScan\test\mixed.fa — multiline FASTA 覆盖。
- c:\GitHub\MotifScan\benchmark\README.md — benchmark 执行说明。
- c:\GitHub\MotifScan\benchmark\benchmark_results.tsv — demo 命令耗时记录。
- c:\GitHub\MotifScan\README.md — 用户文档与方法说明。

**Verification**
1. cargo fmt
2. cargo test
3. cargo build --release
4. ./target/release/motifscan count -i test/illumina.fastq --motifs test/motifs.tsv --revcomp -t 4 -o benchmark/illumina_count.tsv
5. ./target/release/motifscan count -i test/sanger.fa --motif ATTATGAGAATAGTGTG --motif-name Dmel_28 --revcomp -t 2 -o benchmark/sanger_count.tsv
6. ./target/release/motifscan classify -i test/nanopore.fastq --motifs test/rdna_motifs.tsv --revcomp -t 4 -o benchmark/nanopore_classify.tsv --summary benchmark/nanopore_summary.tsv
7. ./target/release/motifscan count -i test/illumina.fastq --motifs test/motifs.tsv --revcomp --min-baseq 20 -t 4 -o benchmark/illumina_q20_count.tsv
8. gzip -c test/illumina.fastq > test/illumina.fastq.gz（若不可用则使用 PowerShell 或 Rust 替代）
9. ./target/release/motifscan count -i test/illumina.fastq.gz --motifs test/motifs.tsv --revcomp -t 4 -o benchmark/illumina_gz_count.tsv
10. 检查 benchmark/benchmark_results.tsv 与各输出文件列顺序、数值逻辑、ambiguous/unclassified 统计是否符合 demo 预期。

**Decisions**
- 第一版默认采用 streaming parser；不在本轮实现 mmap fast path，README 中注明为后续优化。
- 第一版优先支持标准 4 行 FASTQ；multiline FASTQ 不作为当前实现目标，但需在 README 明确限制。
- 多线程采用 chunk + rayon reduce，而不是每条 read 发任务或使用全局锁结构。
- 输出中的 sequence 统一写大写；motif 名称保留用户输入。
- TXT 输出可视为制表分隔的纯文本稳定表格，避免引入额外格式分支复杂度。

**Further Considerations**
1. demo benchmark 记录建议由外层脚本或 PowerShell Measure-Command 采集 wall-clock 秒数，再回填到 benchmark_results.tsv；若实现期需要完全自动化，可追加一个简单脚本，但不应耦合进主程序。
2. 对 classify summary 的 top_ratio，建议定义为 top_reads / total_classified；当 total_classified 为 0 时输出 0 或 NA，需要在实现前固定为一种稳定行为并在 README 说明。