use std::path::Path;
use std::process::Command;
use std::{env, fs};
fn main() {
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

    let mut objects = Vec::new();

    for entry in fs::read_dir(&build_src).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().and_then(|s| s.to_str()) == Some("cpp") {
            if path.to_str().unwrap().ends_with("main.cpp") {
                continue;
            }

            let obj = path.with_extension("o");
            objects.push(obj.clone());

            // println!("cargo:rerun-if-changed={}", path.display());

            /* ---------- 4. g++ 编译 ---------- */
            let status = Command::new("g++")
                .current_dir(&build_src)
                .args([
                    "-O3",
                    "-fPIC",
                    "-std=c++14",
                    "-c",
                    path.file_name().unwrap().to_str().unwrap(),
                    "-o",
                    obj.file_name().unwrap().to_str().unwrap(),
                ])
                .status()
                .expect("failed to run g++");

            if !status.success() {
                panic!("g++ failed on {:?}", path);
            }
        }
    }

    let lib_path = build_out_path.join("libsparc.a");

    let mut ar = Command::new("ar");
    ar.arg("rcs").arg(&lib_path);

    for obj in &objects {
        ar.arg(obj);
    }

    let status = ar.status().expect("failed to run ar");
    if !status.success() {
        panic!("ar failed");
    }

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
