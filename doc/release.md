# MotifScan Packaging and Release (English)

This document explains how to build release-ready binaries from source, create platform archives, generate SHA256 checksums, and publish assets to GitHub Releases.

Note: The repository already includes a CI workflow that builds and publishes release artifacts when the trigger conditions are met. If you just want binaries, download the prebuilt releases from the Releases page.

Prerequisites

- Rust stable toolchain installed (use `rustup`).
- `zip`/`tar` available for packaging (macOS/Windows have built-in tooling).
- Optional but recommended: GitHub CLI `gh` for creating/releases from the command line.

1. Update version

Edit `Cargo.toml` and update the `version` field to the desired semantic version (for example `0.1.2`).

Example:

```toml
[package]
name = "motifscan"
version = "0.1.2"
```

Commit the change:

```bash
git add Cargo.toml
git commit -m "release: bump version to 0.1.2"
```

Note: The CI checks commit messages for the word `release` (case-insensitive) and also treats tag pushes as release triggers. Those conditions decide whether CI will build and publish artifacts.

2. Build release binary locally

```bash
cargo build --release
```

The release binary is located at `target/release/motifscan` (Windows: `target/release/motifscan.exe`).

3. Create platform archives

You can only build reliably for the platform you are on without cross-compilation setup. For multi-platform builds, use the CI runners or build on the target platform.

Linux / macOS:

```bash
mkdir -p release
cp target/release/motifscan release/
# optionally include README and LICENSE
cp README.md release/ || true
cp LICENSE release/ || true
cd release
tar czf ../motifscan-0.1.2-linux.tar.gz motifscan README.md LICENSE
cd ..
```

Windows (on Windows host):

```powershell
New-Item -ItemType Directory -Force -Path release
Copy-Item -Path target\release\motifscan.exe -Destination release\
Copy-Item README.md -Destination release\
Compress-Archive -Path release\* -DestinationPath motifscan-0.1.2-windows.zip
```

4. Generate SHA256 checksum files

Linux / macOS:

```bash
sha256sum motifscan-0.1.2-linux.tar.gz > motifscan-0.1.2-linux.sha256
# verify
sha256sum -c motifscan-0.1.2-linux.sha256
```

Windows (PowerShell):

```powershell
Get-FileHash -Algorithm SHA256 .\motifscan-0.1.2-windows.zip | Format-List
# compare the printed hash with the value saved in the .sha256 file
```

5. Publish to GitHub Releases

Using `gh` (recommended):

```bash
# tag the commit and push tag
git tag v0.1.2
git push origin v0.1.2

# create release and upload binaries and checksums
gh release create v0.1.2 motifscan-0.1.2-linux.tar.gz motifscan-0.1.2-windows.zip motifscan-0.1.2-linux.sha256 --title "v0.1.2" --notes "Release v0.1.2"
```

Or create a Release via GitHub web UI and upload the archives and `.sha256` files manually.

6. Verification

Users who download artifacts can verify SHA256 using `sha256sum` (Linux/macOS) or `Get-FileHash` (PowerShell). The repository README points to these verification commands.

Notes and tips

- For macOS cross-architecture builds (x86_64 vs arm64) build on the matching architecture or use CI macOS runners.
- For automated multi-platform packaging use the repository's GitHub Actions workflow (`.github/workflows/ci.yml`).
- If you want local packaging scripts, consider adding `scripts/package_linux.sh` and `scripts/package_windows.ps1` to standardize packaging steps.

If you want, I can add packaging helper scripts to `scripts/` and a small `Makefile` to simplify local packaging and release tasks.