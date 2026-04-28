//! 中文：命令行接口定义模块，负责把用户输入解析成结构化参数。
//! English: Command-line interface module that parses user input into typed configuration.

use clap::{ArgAction, Args, Parser, Subcommand};

/// 中文：程序顶层 CLI 入口，包含全局参数和子命令。
/// English: Top-level CLI object containing global flags and subcommands.
#[derive(Debug, Clone, Parser)]
#[command(
    name = "motifscan",
    version,
    disable_version_flag = true,
    about = "Streaming motif scanner for FASTA/FASTQ reads"
)]
pub struct Cli {
    #[arg(
        short = 'v',
        long = "version",
        action = ArgAction::SetTrue,
        global = true,
        help = "Print version and citation information"
    )]
    pub version_info: bool,

    #[command(subcommand)]
    pub command: Option<Command>,
}

impl Cli {
    /// 中文：返回当前运行应使用的线程数；如果用户还没进入子命令，就回退到 CPU 核心数。
    /// English: Returns the worker-thread count for the current invocation, falling back to CPU count when no subcommand is selected.
    pub fn threads(&self) -> usize {
        match &self.command {
            Some(Command::Count(args)) => args.threads,
            None => num_cpus::get(),
        }
    }
}

/// 中文：生成版本和引用信息文本，供 `-v/--version` 直接打印。
/// English: Builds the version and citation banner printed by `-v/--version`.
pub fn version_banner() -> String {
    format!(
        "motifscan {}\n\nCitation (BibTeX):\n@software{{motifscan,\n  author = {{jiehua1995}},\n  title = {{MotifScan}},\n  url = {{https://github.com/jiehua1995/MotifScan}},\n  version = {{{}}}\n}}",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_VERSION"),
    )
}

/// 中文：当前支持的子命令集合；现在只有 `count`。
/// English: Set of supported subcommands; currently only `count` is implemented.
#[derive(Debug, Clone, Subcommand)]
pub enum Command {
    Count(CountArgs),
}

/// 中文：`count` 子命令的全部参数，描述输入、motif、线程与输出目标。
/// English: Full argument set for the `count` subcommand, including input, motif, threading, and output targets.
#[derive(Debug, Clone, Args)]
pub struct CountArgs {
    #[arg(short = 'i', long)]
    pub input: std::path::PathBuf,
    #[arg(long, conflicts_with = "motifs")]
    pub motif: Option<String>,
    #[arg(long, default_value = "motif")]
    pub motif_name: String,
    #[arg(long, conflicts_with = "motif")]
    pub motifs: Option<std::path::PathBuf>,
    #[arg(long)]
    pub revcomp: bool,
    #[arg(short = 't', long, default_value_t = num_cpus::get())]
    pub threads: usize,
    #[arg(long)]
    pub progress: bool,
    #[arg(short = 'o', long)]
    pub output: std::path::PathBuf,
    #[arg(long)]
    pub report_read_hits: Option<std::path::PathBuf>,
}

impl CountArgs {
    /// 中文：检查参数组合是否合法，例如必须提供 motif，线程数也不能为 0。
    /// English: Validates argument combinations, ensuring a motif source exists and thread count is non-zero.
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.motif.is_none() && self.motifs.is_none() {
            anyhow::bail!("one of --motif or --motifs is required")
        }
        if self.threads == 0 {
            anyhow::bail!("--threads must be greater than 0")
        }
        Ok(())
    }
}
