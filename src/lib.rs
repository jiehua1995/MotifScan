//! Library entry module that wires together CLI parsing, thread-pool setup, and scan execution.

pub mod cli;
pub mod io;
pub mod motif;
pub mod output;
pub mod scanner;
#[allow(dead_code, clippy::needless_range_loop)]
pub mod species;

use anyhow::Result;
use clap::CommandFactory;
use clap::Parser;
use tracing_subscriber::EnvFilter;

fn print_top_help() {
    println!("motifscan — Streaming motif scanner for FASTA/FASTQ reads\n");
    println!("Usage: motifscan <COMMAND> [OPTIONS]\n");
    println!("Commands:");
    println!("  count      Count exact motif hits in reads");
    println!("  species    Fuzzy long-read locus matching with automatic SNP voting\n");
    println!("Info:");
    println!("  -v, --version    Print version and citation information");
    println!("  -h, --help       Print help\n");
}

fn print_count_help() {
    println!("Count exact motif hits in reads\n");
    println!("Usage: count [OPTIONS]\n");
    println!("Required:");
    println!("  -i, --input <INPUT>    Input read file in FASTA, FASTQ, FASTA.GZ, or FASTQ.GZ format  (required)");
    println!("  -o, --output <OUTPUT>  Output CSV file for motif summary counts  (required)\n");
    println!("Motif:");
    println!("      --motif <MOTIF>            Single motif sequence provided on the command line");
    println!("      --motif-name <MOTIF_NAME>  Output name used with --motif (requires --motif)");
    println!("      --motifs <MOTIFS>          Two-column CSV file containing motif name and sequence");
    println!("      --revcomp                  Also scan the reverse complement of each motif\n");
    println!("Performance:");
    println!("  -t, --threads <THREADS>  Number of worker threads to use [default: auto]\n");
    println!("Behavior:");
    println!("      --progress                             Show a live progress display on stderr");
    println!("      --report-read-hits <REPORT_READ_HITS>  Optional CSV file for read-level hit details\n");
    println!("Info:");
    println!("  -h, --help  Print help");
}

fn print_species_help() {
    println!("Fuzzy long-read locus matching with automatic species-diagnostic SNP voting\n");
    println!("Usage: motifscan species [OPTIONS]\n");
    println!("Core inputs:");
    println!("  -i, --input <INPUT>              Single FASTA/FASTQ(.gz) file");
    println!("      --input-list <INPUT_LIST>    TXT/TSV: path, or sample<TAB>path");
    println!("      --motifs <MOTIFS>            Two-column CSV: name,sequence");
    println!("      --pairs <PAIRS>              Three-column CSV: locus,mel,sim");
    println!("  -o, --output <OUTPUT>            Sample-by-locus summary CSV\n");
    println!("Main thresholds:");
    println!("      --min-shared-identity <F>     Non-diagnostic-site identity [default: 0.85]");
    println!("      --min-aligned-bases <N>       Minimum aligned reference bases [default: 80]");
    println!("      --min-snp-baseq <Q>           Minimum FASTQ SNP Phred [default: 15]");
    println!("      --min-informative-snps <N>    SNPs required for mel/sim call [default: 2]");
    println!("      --species-fraction <F>        Required supporting SNP fraction [default: 0.75]\n");
    println!("Performance / behavior:");
    println!("  -t, --threads <THREADS>           Rayon worker threads [default: auto]");
    println!("      --progress                    Show per-file progress");
    println!("      --anchor-k <K>                Shared anchor length [default: 11]");
    println!("      --anchors-per-locus <N>       Shared anchors/locus [default: 8]");
    println!("      --locus-mode <best|all>       One best locus/read or all passing loci\n");
    println!("Examples:");
    println!("  motifscan species --input-list samples.txt --motifs motifs.csv --pairs pairs.csv -o result.csv -t 48 --progress");
}

pub fn run() -> Result<()> {
    let cli = cli::Cli::parse();
    if cli.version_info {
        println!("{}", cli::version_banner());
        return Ok(());
    }
    if cli.help {
        print_top_help();
        return Ok(());
    }

    init_logging(&cli);
    let thread_count = cli.threads();

    let Some(command) = cli.command else {
        let mut help = cli::Cli::command();
        help.print_help()?;
        println!();
        return Ok(());
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .ok();

    match command {
        cli::Command::Count(args) => {
            if args.help {
                print_count_help();
                return Ok(());
            }
            scanner::run_count(&args)
        }
        cli::Command::Species(args) => {
            if args.help {
                print_species_help();
                return Ok(());
            }
            species::run_species(&args)
        }
    }
}

fn init_logging(_cli: &cli::Cli) {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("warn"));
    let _ = tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .try_init();
}
