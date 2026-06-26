fn main() {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let workspace_root = std::path::PathBuf::from(&manifest_dir)
        .parent() // crates/
        .unwrap()
        .parent() // workspace root
        .unwrap()
        .to_path_buf();
    let third_party = workspace_root.join("third_party");
    let cuda_enabled = std::env::var("CARGO_FEATURE_CUDA").is_ok();

    // Friendly error if the vendored libbsc source is missing. Without this,
    // cc errors out mid-compile with an opaque "file not found" message.
    // (libbsc is now vendored in-tree under third_party/libbsc — see
    // third_party/libbsc/VENDORED.md — so a missing entrypoint means a corrupt
    // or partial checkout, not an un-init'd submodule.)
    let libbsc_entrypoint = third_party.join("libbsc/libbsc/libbsc/libbsc.cpp");
    if !libbsc_entrypoint.exists() {
        panic!(
            "third_party/libbsc is missing or empty — the vendored libbsc source \
             appears to be an incomplete checkout. Restore it from qz's history.\n\
             (looked for {})",
            libbsc_entrypoint.display()
        );
    }
    // Security patch guard: the QZ-SECURITY-PATCH that bounds the BSC decode
    // write against the output-buffer capacity is now tracked in-tree (see
    // third_party/libbsc/VENDORED.md + QZ-SECURITY-PATCH.diff), but this guard
    // stays as defense-in-depth: a future re-vendor/upgrade of libbsc could
    // silently drop the patch. An unpatched libbsc can heap-overflow when
    // decompressing a *crafted* archive (the embedded run/length fields are
    // independent of the header dataSize that sizes the buffer). Fail the build
    // loudly rather than silently shipping the vulnerable codec.
    //
    // The patch spans BOTH decode stages, so guard BOTH:
    //   - coder.cpp/qlfc.cpp — the QLFC entropy stage: re-thread `outputSize` into
    //     bsc_coder_decompress / *_decode_block / bsc_qlfc_*_decode and reject
    //     `n > outputSize` after each `DecodeWord()`.
    //   - lzp/lzp.cpp        — the LZP stage: thread `outputCapacity` into
    //     bsc_lzp_decompress / bsc_lzp_decode_block and bound every match copy +
    //     speculative fast-path store to the output buffer. (LZP blocks are live in
    //     production: columnar headers + fqzcomp sort-keys go through the adaptive
    //     LZP path, and decode is content-agnostic, so a crafted block can select
    //     LZP mode and forge a large match length.)
    for (rel, what) in [
        (
            "libbsc/libbsc/coder/coder.cpp",
            "the QLFC entropy decode stage (coder.cpp/qlfc.cpp)",
        ),
        (
            "libbsc/libbsc/lzp/lzp.cpp",
            "the LZP decode stage (lzp.cpp)",
        ),
    ] {
        let path = third_party.join(rel);
        let patched = std::fs::read_to_string(&path)
            .map(|s| s.contains("QZ-SECURITY-PATCH (bsc-decode-output-bound)"))
            .unwrap_or(false);
        if !patched {
            panic!(
                "third_party/libbsc is missing the QZ-SECURITY-PATCH (bsc-decode-output-bound) \
                 in {what}: {} carries no sentinel. An unpatched libbsc is vulnerable to a heap \
                 buffer overflow on crafted/untrusted archives — re-apply the decode-output-bound \
                 patch before building.",
                path.display()
            );
        }
    }

    let htscodecs_entrypoint = third_party.join("htscodecs/htscodecs/fqzcomp_qual.c");
    if !htscodecs_entrypoint.exists() {
        panic!(
            "third_party/htscodecs (vendored minimal subset) is missing {}. \
             It is tracked in-tree, so a clean checkout should have it — see \
             third_party/htscodecs/VENDORED.md.",
            htscodecs_entrypoint.display()
        );
    }

    // NO_THREADS tls-alloc fallback guard (parallels the libbsc QZ-SECURITY-PATCH
    // guard above). htscodecs is a vendored minimal subset pinned to 69185ce
    // (third_party/htscodecs/VENDORED.md); a future re-vendor could silently change
    // it. We compile htscodecs with `-DNO_THREADS` (see hts_build below) so
    // that `htscodecs_tls_alloc`/`htscodecs_tls_free` in utils.c become the plain
    // `calloc(1, size)` / `free(ptr)` fallback INSTEAD of the per-thread TLS arena
    // that never returns the ~25 MB fqzcomp qual model to the OS (the bounded-RAM
    // fqz decode leak). That fix is ONLY correct as long as the NO_THREADS branch
    // really is calloc/free. If a future re-vendor changes that branch's semantics
    // (e.g. a different pooling scheme), `-DNO_THREADS` would still compile but the
    // per-thread leak could return — peak decode RSS would silently regress from
    // flat to O(workers-touched). Require BOTH the `NO_THREADS` token AND the exact
    // `calloc(1, size)` alloc body + `free(ptr)` free body so a benign reformat is
    // tolerated but a semantic change to the fallback trips the build.
    let utils_c = third_party.join("htscodecs/htscodecs/utils.c");
    let no_threads_calloc_fallback = std::fs::read_to_string(&utils_c)
        .map(|s| {
            s.contains("NO_THREADS")
                && s.contains("return calloc(1, size);")
                && s.contains("free(ptr);")
        })
        .unwrap_or(false);
    if !no_threads_calloc_fallback {
        panic!(
            "third_party/htscodecs/htscodecs/utils.c no longer has the NO_THREADS \
             calloc/free tls-alloc fallback that the bounded-RAM fqz decode fix \
             depends on (htscodecs pinned to 69185ce). A re-vendor may have changed \
             htscodecs_tls_alloc semantics and could reintroduce the ~25MB-per-thread \
             fqz model leak — re-verify reference/single-end decode peak RSS is \
             constant in read count before updating this guard. (looked in {})",
            utils_c.display()
        );
    }

    // CPU target for the bundled C/C++ code. Defaults to `native` (best
    // performance on the build host, matching prior behavior). For portable
    // *distribution* builds — where the binary may run on an older CPU than it
    // was built on and `-march=native` would SIGILL — set QZ_TARGET_CPU to a
    // baseline such as `x86-64-v2` or `x86-64-v3`.
    println!("cargo:rerun-if-env-changed=QZ_TARGET_CPU");
    let target_cpu = std::env::var("QZ_TARGET_CPU").unwrap_or_else(|_| "native".to_string());
    let march_flag = format!("-march={target_cpu}");

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
        .flag_if_supported(&march_flag)
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

    bsc_build.compile("libbsc");

    // Compile CUDA sources (libcubwt + st) with nvcc when cuda feature is enabled
    if cuda_enabled {
        compile_cuda_sources(&libbsc);
    }

    println!(
        "cargo:rerun-if-changed={}",
        third_party.join("libbsc").display()
    );

    // Compile htscodecs (fqzcomp_qual + utils) as static C lib
    let htscodecs = third_party.join("htscodecs/htscodecs");
    let mut hts_build = cc::Build::new();
    hts_build
        .files(&[htscodecs.join("fqzcomp_qual.c"), htscodecs.join("utils.c")])
        .flag_if_supported("-O3")
        // NO_THREADS: disable htscodecs' per-thread TLS arena (utils.c) for the
        // fqzcomp model buffers, falling back to plain calloc/free.
        //
        // WHY: `fqz_create_models` allocates a ~25 MB qual model
        // (sizeof(SIMPLE_MODEL) * CTX_SIZE, CTX_SIZE = 1<<16) via
        // `htscodecs_tls_alloc`. With the pthread TLS arena enabled, that buffer is
        // returned to a PER-THREAD pool on `fqz_destroy_models` and is NEVER
        // reclaimed to the OS — it lives until the owning thread exits. QZ decodes
        // each quality sub-block inside a rayon `par_iter`, so over many chunks the
        // model gets allocated on O(rayon-workers-ever-touched) distinct threads and
        // retained for the whole process. Reference decode (many small per-chunk
        // quality blobs) made decode peak RSS grow ~linearly with chunk/read count
        // (~25 MB × workers-touched), breaking the bounded-RAM property; single-end
        // and paired fqz decode share the same latent leak. The arena's
        // reuse-within-a-thread saving does not survive QZ's per-call rayon thread
        // hops anyway, so disabling it costs ~nothing while making peak fqz RAM
        // bounded by the number of CONCURRENTLY-decoding sub-blocks (≤ pool size),
        // constant in read/chunk count. Measured: reference decode peak RSS goes from
        // O(chunks) (64 chunks ≈ 4 GB) to flat (~80 MB at 4/16/64 chunks).
        .define("NO_THREADS", None)
        .include(&htscodecs)
        .warnings(false);

    hts_build.compile("htscodecs");

    println!("cargo:rerun-if-changed={}", htscodecs.display());

    // Link OpenMP for multithreading.
    // Use rustc-link-search + rustc-link-lib so these propagate to dependent crates.
    // Prefer asking the C compiler where libgomp lives (portable across distros
    // and GCC versions); fall back to the Debian/Ubuntu GCC layout.
    if let Some(dir) = cc_file_dir("libgomp.so") {
        println!("cargo:rustc-link-search=native={}", dir.display());
    } else {
        for entry in std::fs::read_dir("/usr/lib/gcc/x86_64-linux-gnu")
            .into_iter()
            .flatten()
            .flatten()
        {
            if entry.path().join("libgomp.so").exists() {
                println!("cargo:rustc-link-search=native={}", entry.path().display());
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
        let cuda_path =
            std::env::var("CUDA_PATH").unwrap_or_else(|_| "/usr/local/cuda".to_string());
        println!("cargo:rustc-link-search=native={}/lib64", cuda_path);
        println!("cargo:rustc-link-lib=cudart");
    }
}

/// Ask the active C compiler for the absolute path of a runtime file (e.g.
/// `libgomp.so`, `libgcov.a`) via `-print-file-name`. Returns `None` if the
/// compiler can't be run or reports only the bare name (file not found).
fn cc_file_path(file: &str) -> Option<std::path::PathBuf> {
    let compiler = cc::Build::new().get_compiler();
    let output = std::process::Command::new(compiler.path())
        .arg(format!("-print-file-name={file}"))
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let path = std::path::PathBuf::from(String::from_utf8_lossy(&output.stdout).trim());
    // GCC echoes the bare filename when it can't locate the file.
    if path.is_absolute() && path.exists() {
        Some(path)
    } else {
        None
    }
}

/// Directory containing a compiler runtime file, via [`cc_file_path`].
fn cc_file_dir(file: &str) -> Option<std::path::PathBuf> {
    cc_file_path(file).and_then(|p| p.parent().map(|d| d.to_path_buf()))
}

/// Compile CUDA .cu sources with nvcc and produce a static library.
fn compile_cuda_sources(libbsc: &std::path::Path) {
    use std::process::Command;

    let out_dir = std::env::var("OUT_DIR").unwrap();
    let out_path = std::path::PathBuf::from(&out_dir);

    let cuda_path = std::env::var("CUDA_PATH").unwrap_or_else(|_| "/usr/local/cuda".to_string());
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
                "-Xcompiler",
                "-fPIC",
                "-DLIBBSC_CUDA_SUPPORT",
                "-DLIBBSC_OPENMP_SUPPORT",
                "-I",
                libbsc.to_str().unwrap(),
                "-o",
                obj_file.to_str().unwrap(),
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
