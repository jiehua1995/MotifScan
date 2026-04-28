//! 中文：库入口模块，负责串联 CLI、线程池初始化和扫描执行流程。
//! English: Library entry module that wires together CLI parsing, thread-pool setup, and scan execution.

pub mod cli;
pub mod io;
pub mod motif;
pub mod output;
pub mod scanner;

use anyhow::Result;
use clap::CommandFactory;
use clap::Parser;

/// 中文：运行程序主流程。
/// English: Runs the main application flow.
///
/// 中文：这个函数先解析命令行参数，再处理版本输出、帮助文本、线程池初始化，最后把控制权交给具体子命令。
/// English: This function parses CLI arguments, handles version/help output, initializes the Rayon thread pool, and finally dispatches to the selected subcommand.
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
