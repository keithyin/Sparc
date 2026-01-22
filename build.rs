use std::num::NonZero;
use std::path::Path;
use std::process::Command;
use std::{env, fs};

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

    Command::new("sh")
        .arg("-c")
        .arg(&format!(
            "cp -r {} {}",
            source_code_dir,
            build_out_path.to_str().unwrap()
        ))
        .output()
        .expect(&format!("cp {} error", source_code_dir));

    Command::new("make")
        .current_dir(&build_src)
        .arg(format!("-j{}", num_cpus))
        .output()
        .expect("make error");

    /* ---------- 6. 告诉 Rust 去哪里 link ---------- */
    println!(
        "cargo:rustc-link-search=native={}",
        build_out_path.display()
    );
    println!("cargo:rustc-link-lib=static=sparc");

    /* ---------- 7. C++ 标准库（非常重要） ---------- */
    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=c++");
    } else {
        println!("cargo:rustc-link-lib=stdc++");
    }
}
