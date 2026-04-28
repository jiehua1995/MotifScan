pub mod cli;
pub mod io;
pub mod motif;
pub mod output;
pub mod scanner;

use anyhow::Result;
use clap::CommandFactory;
use clap::Parser;

pub fn run() -> Result<()> {
    let cli = cli::Cli::parse();
    if cli.version_info {
        println!("{}", cli::version_banner());
        return Ok(());
    }

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
        cli::Command::Count(args) => scanner::run_count(&args),
    }
}
