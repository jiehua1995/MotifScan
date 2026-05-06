# To-Do List

这是为 MotifScan 制定的可执行改进计划（按优先级列出）。每项将逐一完成并在实现后打勾。

## 高优先级

- [x] 流式写入 read-hit（写线程 + channel）：为 `--report-read-hits` 创建单独写入线程，接收并写出匹配批次，避免在内存中积累大量 `ReadHitRow`。
- [x] 引入 Aho–Corasick 多模式路径：当 motif 数量较多时显著提升速度。

## 中优先级

 - [x] 在 CI 中添加 `cargo fmt -- --check`、`cargo clippy`、`cargo test`，并建立 release workflow。
 - [x] 在 CI 中添加 `cargo fmt -- --check`、`cargo clippy`、`cargo test`，并建立 release workflow（产物使用 `Cargo.toml` 中的语义版本号并包含平台信息）。
 - [x] 为每个平台生成 SHA256 校验文件并随 release 上传
## 低优先级 / 可选

- [ ] 支持输出压缩（例如输出路径以 `.gz` 自动压缩）。
- [ ] 支持 `-i -` 从 stdin 读取和将输出写入 stdout（便于 pipeline）。
- [ ] 暴露 `--chunk-size` / `--batch-size` 配置并尝试自动调优默认值。
- [ ] 增加性能基准（`criterion`）并记录基线。

## 测试数据与示例

- [x] 在仓库根目录下添加 `test/` 目录，包含用于回归和集成测试的模拟 `fastq.gz` 文件与 `motifs.csv`。
- [x] 添加集成测试 `tests/integration.rs`，验证 `run_count` 在示例数据上能生成输出。

---

下一步：我会创建 `test/generate_fastq.py`（生成压缩的示例 FASTQ）和 `test/motifs.csv`，然后运行脚本生成 `test/*.fastq.gz`。随后我会开始实现第一项“流式写入 read-hit”。