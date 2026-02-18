fn main() {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let workspace_root = std::path::PathBuf::from(&manifest_dir)
        .parent()  // crates/
        .unwrap()
        .parent()  // workspace root
        .unwrap()
        .to_path_buf();
    let third_party = workspace_root.join("third_party");
    let cuda_enabled = std::env::var("CARGO_FEATURE_CUDA").is_ok();

    // C/C++ PGO support:
    //   QZ_PGO_GENERATE=/path/to/dir  → instrument for profile collection
    //   QZ_PGO_USE=/path/to/dir       → optimize using collected profiles
    let pgo_generate = std::env::var("QZ_PGO_GENERATE").ok();
    let pgo_use = std::env::var("QZ_PGO_USE").ok();

    // Compile libbsc C++ library together with libsais in one build
    let libbsc = third_party.join("libbsc/libbsc");
    let mut bsc_build = cc::Build::new();
    bsc_build
        .cpp(true)
        .files(&[
            // libsais (C source - dependency of libbsc)
            libbsc.join("bwt/libsais/libsais.c"),
            // libbsc C++ sources
            libbsc.join("libbsc/libbsc.cpp"),
            libbsc.join("adler32/adler32.cpp"),
            libbsc.join("bwt/bwt.cpp"),
            libbsc.join("coder/coder.cpp"),
            libbsc.join("coder/qlfc/qlfc.cpp"),
            libbsc.join("coder/qlfc/qlfc_model.cpp"),
            libbsc.join("filters/detectors.cpp"),
            libbsc.join("filters/preprocessing.cpp"),
            libbsc.join("lzp/lzp.cpp"),
            libbsc.join("platform/platform.cpp"),
            libbsc.join("st/st.cpp"),
        ])
        .flag_if_supported("-O3")
        .flag_if_supported("-march=native")
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-fopenmp")
        .define("LIBBSC_OPENMP_SUPPORT", None)
        .define("LIBSAIS_OPENMP", None)
        .include(&libbsc)
        .include(libbsc.join("bwt/libsais"))
        .warnings(false);

    if cuda_enabled {
        bsc_build.define("LIBBSC_CUDA_SUPPORT", None);
    }

    if let Some(ref dir) = pgo_generate {
        bsc_build.flag(&format!("-fprofile-generate={}", dir));
    }
    if let Some(ref dir) = pgo_use {
        bsc_build.flag(&format!("-fprofile-use={}", dir));
        bsc_build.flag("-fprofile-correction");
    }
    bsc_build.compile("libbsc");

    // Compile CUDA sources (libcubwt + st) with nvcc when cuda feature is enabled
    if cuda_enabled {
        compile_cuda_sources(&libbsc);
    }

    println!("cargo:rerun-if-changed={}", third_party.join("libbsc").display());

    // Compile htscodecs (fqzcomp_qual + utils) as static C lib
    let htscodecs = third_party.join("htscodecs/htscodecs");
    let mut hts_build = cc::Build::new();
    hts_build
        .files(&[
            htscodecs.join("fqzcomp_qual.c"),
            htscodecs.join("utils.c"),
        ])
        .flag_if_supported("-O3")
        .include(&htscodecs)
        .warnings(false);

    if let Some(ref dir) = pgo_generate {
        hts_build.flag(&format!("-fprofile-generate={}", dir));
    }
    if let Some(ref dir) = pgo_use {
        hts_build.flag(&format!("-fprofile-use={}", dir));
        hts_build.flag("-fprofile-correction");
    }
    hts_build.compile("htscodecs");

    println!("cargo:rerun-if-changed={}", htscodecs.display());

    // Link OpenMP for multithreading
    // Use rustc-link-search + rustc-link-lib so these propagate to dependent crates
    for entry in std::fs::read_dir("/usr/lib/gcc/x86_64-linux-gnu").into_iter().flatten() {
        if let Ok(e) = entry {
            if e.path().join("libgomp.so").exists() {
                println!("cargo:rustc-link-search=native={}", e.path().display());
            }
        }
    }
    println!("cargo:rustc-link-lib=gomp");
    println!("cargo:rustc-link-lib=stdc++");
    // Link math lib for htscodecs (log, etc.)
    println!("cargo:rustc-link-lib=m");

    // Link CUDA runtime when cuda feature is enabled
    if cuda_enabled {
        // Find CUDA toolkit
        let cuda_path = std::env::var("CUDA_PATH")
            .unwrap_or_else(|_| "/usr/local/cuda".to_string());
        println!("cargo:rustc-link-search=native={}/lib64", cuda_path);
        println!("cargo:rustc-link-lib=cudart");
    }

    // Link gcov runtime for C/C++ PGO instrumentation
    if pgo_generate.is_some() {
        println!("cargo:rustc-link-arg=-Wl,--whole-archive");
        println!("cargo:rustc-link-arg=/usr/lib/gcc/x86_64-linux-gnu/13/libgcov.a");
        println!("cargo:rustc-link-arg=-Wl,--no-whole-archive");
    }
}

/// Compile CUDA .cu sources with nvcc and produce a static library.
fn compile_cuda_sources(libbsc: &std::path::Path) {
    use std::process::Command;

    let out_dir = std::env::var("OUT_DIR").unwrap();
    let out_path = std::path::PathBuf::from(&out_dir);

    let cuda_path = std::env::var("CUDA_PATH")
        .unwrap_or_else(|_| "/usr/local/cuda".to_string());
    let nvcc = std::path::PathBuf::from(&cuda_path).join("bin/nvcc");

    if !nvcc.exists() {
        panic!(
            "CUDA feature enabled but nvcc not found at {}. \
             Set CUDA_PATH or install CUDA toolkit.",
            nvcc.display()
        );
    }

    let cu_files = [
        libbsc.join("bwt/libcubwt/libcubwt.cu"),
        libbsc.join("st/st.cu"),
    ];

    let mut obj_files = Vec::new();

    for cu_file in &cu_files {
        let stem = cu_file.file_stem().unwrap().to_str().unwrap();
        let obj_file = out_path.join(format!("{}.o", stem));

        let status = Command::new(&nvcc)
            .args([
                "-O3",
                "-c",
                "-Xcompiler", "-fPIC",
                "-DLIBBSC_CUDA_SUPPORT",
                "-DLIBBSC_OPENMP_SUPPORT",
                "-I", libbsc.to_str().unwrap(),
                "-o", obj_file.to_str().unwrap(),
                cu_file.to_str().unwrap(),
            ])
            .status()
            .unwrap_or_else(|e| panic!("Failed to run nvcc: {}", e));

        if !status.success() {
            panic!("nvcc compilation failed for {}", cu_file.display());
        }

        obj_files.push(obj_file);
    }

    // Create static library from CUDA object files
    let cuda_lib = out_path.join("liblibbsc_cuda.a");
    let status = Command::new("ar")
        .args(["rcs", cuda_lib.to_str().unwrap()])
        .args(obj_files.iter().map(|p| p.to_str().unwrap()))
        .status()
        .unwrap_or_else(|e| panic!("Failed to run ar: {}", e));

    if !status.success() {
        panic!("ar failed to create CUDA static library");
    }

    println!("cargo:rustc-link-search=native={}", out_dir);
    println!("cargo:rustc-link-lib=static=libbsc_cuda");
}
