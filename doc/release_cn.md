# MotifScan 打包与发布（中文）

本文件说明如何从源码在本地打包预编译的发布包，并生成 SHA256 校验文件，最后通过 GitHub Releases 发布。

注意：仓库已配置 CI，会在满足触发条件时自动打包并生成校验文件；如果你只是想下载二进制，请直接使用 Releases 页面。

## 前置条件

- 已安装 Rust（stable toolchain），建议通过 `rustup` 管理。
- 已安装 `zip`（Windows 可用 PowerShell 自带压缩或使用 `Compress-Archive`）。
- 推荐安装 `gh`（GitHub CLI）以便在命令行创建 release。也可直接通过网页创建 Release 并上传文件。

## 1. 更新版本号

在 `Cargo.toml` 中更新 `version = "x.y.z"`。遵循语义版本规则。

例如：

```toml
[package]
name = "motifscan"
version = "0.1.2"

```

提交改动：

```bash
git add Cargo.toml
git commit -m "release: bump version to 0.1.2"
```

> 注意：仓库 CI 配置中会检测提交信息是否包含 `release`（不区分大小写），或是否是 tag 推送，用于决定是否在 CI 上执行打包与发布。

## 2. 本地构建 Release 二进制

```bash
cargo build --release
```

构建完成后，二进制位于：`target/release/motifscan`（Windows 为 `target/release/motifscan.exe`）。

## 3. 为不同平台制作包

在本地你通常只能为当前运行的平台构建；交叉编译更复杂，建议在 CI 中为目标平台构建（或在对应平台机器上打包）。下面示例显示常见打包方式：

Linux / macOS：

```bash
mkdir -p release
cp target/release/motifscan release/
# 可选：包含 LICENSE 与 README
cp README.md release/
cp LICENSE release/ || true
# 创建 tar.gz
cd release
tar czf ../motifscan-0.1.2-linux.tar.gz motifscan README.md LICENSE
cd ..
```

Windows（在 Windows 主机上）：

```powershell
New-Item -ItemType Directory -Force -Path release
Copy-Item -Path target\release\motifscan.exe -Destination release\
Copy-Item README.md -Destination release\
# 使用 Compress-Archive 打包
Compress-Archive -Path release\* -DestinationPath motifscan-0.1.2-windows.zip
```

## 4. 生成 SHA256 校验文件

Linux / macOS：

```bash
sha256sum motifscan-0.1.2-linux.tar.gz > motifscan-0.1.2-linux.sha256
# 验证
sha256sum -c motifscan-0.1.2-linux.sha256
```

Windows（PowerShell）：

```powershell
Get-FileHash -Algorithm SHA256 .\motifscan-0.1.2-windows.zip | Format-List
# 手动比对输出的哈希值与 .sha256 内容
```

## 5. 在 GitHub 上创建 Release 并上传文件

推荐使用 `gh`（GitHub CLI）：

```bash
# 使用 tag 创建 release 并上传文件
git tag v0.1.2
git push origin v0.1.2
# 使用 gh 创建 release（如果尚未登录，先 gh auth login）
gh release create v0.1.2 motifscan-0.1.2-linux.tar.gz motifscan-0.1.2-windows.zip motifscan-0.1.2-linux.sha256 --title "v0.1.2" --notes "Release v0.1.2"
```

或者你可以在 GitHub 网页上创建 Release，然后手动上传打包文件和 `.sha256` 校验文件。

## 6. 验证已上传的文件

下载后用户可以使用前面提到的 `sha256sum -c` 或 `Get-FileHash` 来验证文件完整性。

## 常见问题

- 如何为 macOS/x86_64 或 macOS/arm64 打包？
  - 建议在对应硬件或 CI runner（macOS-latest）上构建并打包。交叉编译 macOS 目标通常需要额外工具链和交叉编译设置。

- 我如何在 CI 中自动化上述步骤？
  - 本仓库已包含 GitHub Actions workflow（.github/workflows/ci.yml），会在满足触发条件时构建并生成平台包、计算 SHA256，并作为 Release 附件上传。

---

若需要，我可以把上面的打包脚本加入仓库（例如 `scripts/package_linux.sh`、`scripts/package_windows.ps1`），并帮助把 `gh` 命令封装为可重用脚本。