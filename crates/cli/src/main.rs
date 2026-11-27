#[cfg(target_env = "musl")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> std::process::ExitCode {
    bammetrics_cli::remove_files_on_signal();
    bammetrics_cli::main()
}
