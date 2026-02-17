fn main() {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let workspace_root = std::path::PathBuf::from(&manifest_dir)
        .parent()  // crates/
        .unwrap()
        .parent()  // workspace root
        .unwrap()
        .to_path_buf();
    let third_party = workspace_root.join("third_party");

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

    if let Some(ref dir) = pgo_generate {
        bsc_build.flag(&format!("-fprofile-generate={}", dir));
    }
    if let Some(ref dir) = pgo_use {
        bsc_build.flag(&format!("-fprofile-use={}", dir));
        bsc_build.flag("-fprofile-correction");
    }
    bsc_build.compile("libbsc");

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

    // Link gcov runtime for C/C++ PGO instrumentation
    if pgo_generate.is_some() {
        println!("cargo:rustc-link-arg=-Wl,--whole-archive");
        println!("cargo:rustc-link-arg=/usr/lib/gcc/x86_64-linux-gnu/13/libgcov.a");
        println!("cargo:rustc-link-arg=-Wl,--no-whole-archive");
    }
}
