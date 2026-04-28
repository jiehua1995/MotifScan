/// 中文：二进制程序入口，只负责调用库层主流程并统一处理错误退出码。
/// English: Binary entrypoint that simply calls the library runner and normalizes fatal error handling.
fn main() {
    if let Err(error) = motifscan::run() {
        eprintln!("error: {error:#}");
        std::process::exit(1);
    }
}
