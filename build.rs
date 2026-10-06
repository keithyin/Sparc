use std::num::NonZero;
use std::path::Path;
use std::process::Command;
use std::{env, fs};

/// 递归复制目录（替代 `sh -c cp -r`，不依赖 shell，保证可移植）。
fn copy_dir_recursive(src: &Path, dst: &Path) -> std::io::Result<()> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let dst_path = dst.join(entry.file_name());
        if entry.file_type()?.is_dir() {
            copy_dir_recursive(&entry.path(), &dst_path)?;
        } else {
            fs::copy(entry.path(), &dst_path)?;
        }
    }
    Ok(())
}

/// 运行外部命令，非零退出码立刻失败并带出 stdout/stderr，
/// 避免 C++ 编译错误被延迟成莫名的链接错误。
fn run_or_panic(cmd: &mut Command, what: &str) {
    let output = cmd
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn {what}: {e}"));
    if !output.status.success() {
        panic!(
            "{what} failed with {output_status}\nstdout:\n{stdout}\nstderr:\n{stderr}",
            output_status = output.status,
            stdout = String::from_utf8_lossy(&output.stdout),
            stderr = String::from_utf8_lossy(&output.stderr),
        );
    }
}

fn main() {
    let num_cpus: usize = std::thread::available_parallelism()
        .unwrap_or(NonZero::new(10).unwrap())
        .into();
    let source_code_dir = "sparc-source-code";
    let out_path_str = env::var("OUT_DIR").unwrap();
    let build_out_path = Path::new(&out_path_str);

    let build_src = build_out_path.join("sparc-source-code");

    println!("cargo:rerun-if-changed={}", source_code_dir);

    if build_src.exists() {
        fs::remove_dir_all(&build_src).unwrap();
    }

    copy_dir_recursive(Path::new(source_code_dir), &build_src)
        .unwrap_or_else(|e| panic!("copy {source_code_dir} failed: {e}"));

    let mut make = Command::new("make");
    make.current_dir(&build_src).arg(format!("-j{}", num_cpus));
    run_or_panic(&mut make, "make");

    /* ---------- 告诉 Rust 去哪里 link ---------- */
    println!("cargo:rustc-link-search=native={}", build_src.display());
    println!("cargo:rustc-link-lib=static=sparc");

    /* ---------- C++ 标准库（非常重要） ---------- */
    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=c++");
    } else {
        println!("cargo:rustc-link-lib=stdc++");
    }
}
