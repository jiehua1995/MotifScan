use clap::{ArgAction, Args, Parser, Subcommand};

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
    pub fn threads(&self) -> usize {
        match &self.command {
            Some(Command::Count(args)) => args.threads,
            None => num_cpus::get(),
        }
    }
}

pub fn version_banner() -> String {
    format!(
        "motifscan {}\nCitation: MotifScan {}, streaming Rust CLI for motif scanning. When used in analysis or manuscripts, report the software name, version, and the exact source revision or archive you ran.",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_VERSION")
    )
}

#[derive(Debug, Clone, Subcommand)]
pub enum Command {
    Count(CountArgs),
}

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
    #[arg(long)]
    pub iupac: bool,
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
