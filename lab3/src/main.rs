use faer::Mat;
use hmatrix_svd::{
    calculate_psnr, calculate_rgb_metrics, compute_thin_svd, create_tree, decompress_hmatrix,
    draw_hmatrix_structure, frobenius_error, hmatrix_matvec, save_cases_csv, save_matrix_as_image,
    save_results_json, save_rgb_image, save_scaling_csv, save_small_tests_json, CaseMetrics,
    FullExperimentResults, ScalingPoint, SingularValueData, SmallTestMatrixData, SweepPoint,
};
use image::GenericImageView;
use std::fs;
use std::path::Path;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Start Hierarchical Matrix Compression (H-Matrix / rSVD) ===");

    fs::create_dir_all("results")?;
    fs::create_dir_all("typst/imgs")?;

    let img_path = "assets/input_512.png";
    if !Path::new(img_path).exists() {
        println!("Input image {} not found!", img_path);
    }

    let dyn_img = image::open(img_path)?;
    let (width, height) = dyn_img.dimensions();
    println!("Loaded image: {}x{} pixels", width, height);
    assert_eq!(width, 512);
    assert_eq!(height, 512);

    let n = 512;
    let mut r_mat = Mat::<f64>::zeros(n, n);
    let mut g_mat = Mat::<f64>::zeros(n, n);
    let mut b_mat = Mat::<f64>::zeros(n, n);

    for (x, y, pixel) in dyn_img.pixels() {
        let (i, j) = (y as usize, x as usize);
        r_mat.write(i, j, pixel[0] as f64);
        g_mat.write(i, j, pixel[1] as f64);
        b_mat.write(i, j, pixel[2] as f64);
    }

    // Save original RGB image and individual channels
    println!("Saving uncompressed images...");
    dyn_img.save("typst/imgs/orig_rgb.png")?;
    save_matrix_as_image(&r_mat, "typst/imgs/orig_r.png")?;
    save_matrix_as_image(&g_mat, "typst/imgs/orig_g.png")?;
    save_matrix_as_image(&b_mat, "typst/imgs/orig_b.png")?;

    // Theoretical Peak GFLOP/s for Intel Core i5-1235U (1 core AVX2 FMA 4.40 GHz: 35.2 GFLOP/s)
    let peak_gflops = 35.2;
    let cpu_model = "12th Gen Intel Core i5-1235U (4.40 GHz, AVX2/FMA)".to_string();

    // 1. Full SVD of individual channels
    println!("Computing full SVD for R, G, B channels...");
    let svd_r = compute_thin_svd(&r_mat);
    let svd_g = compute_thin_svd(&g_mat);
    let svd_b = compute_thin_svd(&b_mat);

    let extract_sv = |svd: &faer::linalg::solvers::ThinSvd<f64>, name: &str| -> SingularValueData {
        let s = svd.s_diagonal();
        let mut vals = Vec::with_capacity(n);
        for i in 0..n {
            vals.push(s.read(i));
        }
        SingularValueData {
            channel: name.to_string(),
            values: vals.clone(),
            sigma_1: vals[0],
            sigma_256: vals[255],
            sigma_512: vals[511],
        }
    };

    let sv_data_r = extract_sv(&svd_r, "Red");
    let sv_data_g = extract_sv(&svd_g, "Green");
    let sv_data_b = extract_sv(&svd_b, "Blue");

    println!(
        "Red channel:   sigma_1 = {:.3}, sigma_256 = {:.3}, sigma_512 = {:.3}",
        sv_data_r.sigma_1, sv_data_r.sigma_256, sv_data_r.sigma_512
    );
    println!(
        "Green channel: sigma_1 = {:.3}, sigma_256 = {:.3}, sigma_512 = {:.3}",
        sv_data_g.sigma_1, sv_data_g.sigma_256, sv_data_g.sigma_512
    );
    println!(
        "Blue channel:  sigma_1 = {:.3}, sigma_256 = {:.3}, sigma_512 = {:.3}",
        sv_data_b.sigma_1, sv_data_b.sigma_256, sv_data_b.sigma_512
    );

    // 2. Eight Cases (6 mandatory from instruction + 2 optimal)
    let case_configs = vec![
        ("case1", 1, sv_data_r.sigma_1, sv_data_g.sigma_1, sv_data_b.sigma_1),
        ("case2", 1, sv_data_r.sigma_512, sv_data_g.sigma_512, sv_data_b.sigma_512),
        ("case3", 1, sv_data_r.sigma_256, sv_data_g.sigma_256, sv_data_b.sigma_256),
        ("case4", 4, sv_data_r.sigma_1, sv_data_g.sigma_1, sv_data_b.sigma_1),
        ("case5", 4, sv_data_r.sigma_512, sv_data_g.sigma_512, sv_data_b.sigma_512),
        ("case6", 4, sv_data_r.sigma_256, sv_data_g.sigma_256, sv_data_b.sigma_256),
        ("opt1", 2, 3.0 * sv_data_r.sigma_256, 3.0 * sv_data_g.sigma_256, 3.0 * sv_data_b.sigma_256),
        ("opt2", 4, 4.0 * sv_data_r.sigma_256, 4.0 * sv_data_g.sigma_256, 4.0 * sv_data_b.sigma_256),
    ];

    let mut case_results = Vec::new();

    for (cname, r_rank, delta_r, delta_g, delta_b) in case_configs {
        println!(
            "--> Running {}: r = {}, delta = [R: {:.2}, G: {:.2}, B: {:.2}]",
            cname, r_rank, delta_r, delta_g, delta_b
        );

        let mut comp_flops = 0u64;
        let t_comp_start = Instant::now();

        let tree_r = create_tree(&r_mat, 0, n, 0, n, r_rank, delta_r, &mut comp_flops);
        let tree_g = create_tree(&g_mat, 0, n, 0, n, r_rank, delta_g, &mut comp_flops);
        let tree_b = create_tree(&b_mat, 0, n, 0, n, r_rank, delta_b, &mut comp_flops);

        let t_comp_elapsed = t_comp_start.elapsed().as_secs_f64();
        let comp_time_ms = t_comp_elapsed * 1000.0;
        let comp_gflops = (comp_flops as f64 / 1e9) / t_comp_elapsed.max(1e-9);
        let comp_peak_percent = (comp_gflops / peak_gflops) * 100.0;

        // Tree metrics
        let (z_r, l_r, _) = tree_r.node_counts();
        let (z_g, l_g, _) = tree_g.node_counts();
        let (z_b, l_b, _) = tree_b.node_counts();

        let d_r = tree_r.depth();
        let d_g = tree_g.depth();
        let d_b = tree_b.depth();

        let stored_floats_r = tree_r.stored_floats();
        let stored_floats_g = tree_g.stored_floats();
        let stored_floats_b = tree_b.stored_floats();
        let stored_floats_total = stored_floats_r + stored_floats_g + stored_floats_b;

        let dense_floats_total = 3 * n * n;
        let compression_ratio = dense_floats_total as f64 / stored_floats_total.max(1) as f64;
        let memory_savings =
            (1.0 - (stored_floats_total as f64 / dense_floats_total as f64)) * 100.0;

        let tree_mem_bytes = tree_r.memory_bytes() + tree_g.memory_bytes() + tree_b.memory_bytes();
        let tree_mem_kb = tree_mem_bytes as f64 / 1024.0;
        let dense_mem_kb = (dense_floats_total * 8) as f64 / 1024.0;

        // Decompression
        let mut decomp_flops = 0u64;
        let mut rec_r = Mat::<f64>::zeros(n, n);
        let mut rec_g = Mat::<f64>::zeros(n, n);
        let mut rec_b = Mat::<f64>::zeros(n, n);

        let t_decomp_start = Instant::now();
        decompress_hmatrix(&tree_r, &mut rec_r, &mut decomp_flops);
        decompress_hmatrix(&tree_g, &mut rec_g, &mut decomp_flops);
        decompress_hmatrix(&tree_b, &mut rec_b, &mut decomp_flops);
        let decomp_time_ms = t_decomp_start.elapsed().as_secs_f64() * 1000.0;

        // Reconstruction metrics
        let (_, _, rel_frob_r) = frobenius_error(&r_mat, &rec_r);
        let (_, _, rel_frob_g) = frobenius_error(&g_mat, &rec_g);
        let (_, _, rel_frob_b) = frobenius_error(&b_mat, &rec_b);

        let psnr_r = calculate_psnr(&r_mat, &rec_r);
        let psnr_g = calculate_psnr(&g_mat, &rec_g);
        let psnr_b = calculate_psnr(&b_mat, &rec_b);

        let (rel_frob_rgb, psnr_rgb) =
            calculate_rgb_metrics(&r_mat, &g_mat, &b_mat, &rec_r, &rec_g, &rec_b);

        // Save H-matrix structure images (typst/imgs/{case}_tree_{r,g,b}.png)
        draw_hmatrix_structure(&tree_r, n, n, 512, &format!("typst/imgs/{}_tree_r.png", cname))?;
        draw_hmatrix_structure(&tree_g, n, n, 512, &format!("typst/imgs/{}_tree_g.png", cname))?;
        draw_hmatrix_structure(&tree_b, n, n, 512, &format!("typst/imgs/{}_tree_b.png", cname))?;

        // Save reconstructed channel bitmaps
        save_matrix_as_image(&rec_r, &format!("typst/imgs/{}_rec_r.png", cname))?;
        save_matrix_as_image(&rec_g, &format!("typst/imgs/{}_rec_g.png", cname))?;
        save_matrix_as_image(&rec_b, &format!("typst/imgs/{}_rec_b.png", cname))?;

        // Save composite reconstructed RGB image
        save_rgb_image(&rec_r, &rec_g, &rec_b, &format!("typst/imgs/{}_rec_rgb.png", cname))?;

        println!(
            "   CR = {:.2}x ({:.1}% memory saved), PSNR RGB = {:.2} dB, Time = {:.1} ms",
            compression_ratio, memory_savings, psnr_rgb, comp_time_ms
        );

        case_results.push(CaseMetrics {
            name: cname.to_string(),
            r: r_rank,
            delta_r,
            delta_g,
            delta_b,
            zero_leaves_r: z_r,
            rank_leaves_r: l_r,
            depth_r: d_r,
            zero_leaves_g: z_g,
            rank_leaves_g: l_g,
            depth_g: d_g,
            zero_leaves_b: z_b,
            rank_leaves_b: l_b,
            depth_b: d_b,
            stored_floats_total,
            dense_floats_total,
            compression_ratio,
            memory_savings_percent: memory_savings,
            tree_memory_kb: tree_mem_kb,
            dense_memory_kb: dense_mem_kb,
            rel_frob_r,
            rel_frob_g,
            rel_frob_b,
            rel_frob_rgb,
            psnr_r,
            psnr_g,
            psnr_b,
            psnr_rgb,
            compress_time_ms: comp_time_ms,
            decompress_time_ms: decomp_time_ms,
            compress_flops: comp_flops,
            decompress_flops: decomp_flops,
            compress_gflops: comp_gflops,
            compress_peak_percent: comp_peak_percent,
        });
    }

    // 3. Rate-Distortion Parameter Sweep (Memory & Quality vs Delta & Rank)
    println!("Performing Parameter Sweep for Rate-Distortion curves...");
    let ranks = vec![1, 2, 4, 8];
    let multipliers = vec![0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0, 40.0];
    let mut sweep_results = Vec::new();

    for &r_rank in &ranks {
        for &mult in &multipliers {
            let delta = mult * sv_data_r.sigma_256;
            let mut flops = 0u64;
            let tree = create_tree(&r_mat, 0, n, 0, n, r_rank, delta, &mut flops);
            let mut rec = Mat::<f64>::zeros(n, n);
            let mut d_flops = 0u64;
            decompress_hmatrix(&tree, &mut rec, &mut d_flops);

            let stored = tree.stored_floats();
            let cr = (n * n) as f64 / stored.max(1) as f64;
            let mem_kb = tree.memory_bytes() as f64 / 1024.0;
            let (_, _, rel_frob) = frobenius_error(&r_mat, &rec);
            let psnr = calculate_psnr(&r_mat, &rec);

            sweep_results.push(SweepPoint {
                r: r_rank,
                delta_mult: mult,
                delta,
                stored_floats: stored,
                compression_ratio: cr,
                memory_kb: mem_kb,
                rel_frob,
                psnr,
            });
        }
    }

    // 4. Scaling with Matrix Dimension N: time, FLOPs, GFLOP/s, % Peak CPU
    println!("Benchmarking scaling with dimension N in [32, 64, 128, 256, 512, 1024]...");
    let test_ns = vec![32, 64, 128, 256, 512, 1024];
    let mut scaling_results = Vec::new();

    for &sz in &test_ns {
        let sub_mat = if sz <= 512 {
            Mat::from_fn(sz, sz, |i, j| r_mat.read(i, j))
        } else {
            Mat::from_fn(sz, sz, |i, j| {
                let base = r_mat.read(i % 512, j % 512);
                let factor = 1.0 + 0.05 * ((i / 512 + j / 512) as f64);
                base * factor
            })
        };

        // Full SVD
        let svd_flops = (4 * sz * sz * sz + 8 * sz * sz * sz + 9 * sz * sz * sz) as u64;
        let t0 = Instant::now();
        let svd_full = compute_thin_svd(&sub_mat);
        let t_svd = t0.elapsed().as_secs_f64();
        let gflops_svd = (svd_flops as f64 / 1e9) / t_svd.max(1e-9);
        let peak_svd = (gflops_svd / peak_gflops) * 100.0;

        let mid_sigma = svd_full.s_diagonal().read(sz / 2);

        // H-Matrix r = 1
        let mut flops_r1 = 0u64;
        let t1 = Instant::now();
        let tree_r1 = create_tree(&sub_mat, 0, sz, 0, sz, 1, mid_sigma, &mut flops_r1);
        let t_r1 = t1.elapsed().as_secs_f64();
        let gflops_r1 = (flops_r1 as f64 / 1e9) / t_r1.max(1e-9);
        let peak_r1 = (gflops_r1 / peak_gflops) * 100.0;

        // H-Matrix r = 4
        let mut flops_r4 = 0u64;
        let t2 = Instant::now();
        let tree_r4 = create_tree(&sub_mat, 0, sz, 0, sz, 4, 3.0 * mid_sigma, &mut flops_r4);
        let t_r4 = t2.elapsed().as_secs_f64();
        let gflops_r4 = (flops_r4 as f64 / 1e9) / t_r4.max(1e-9);
        let peak_r4 = (gflops_r4 / peak_gflops) * 100.0;

        // Matvec benchmark
        let x_vec = vec![1.0; sz];
        let mut y_dense = vec![0.0; sz];
        let mut y_hmat = vec![0.0; sz];

        let iters = (10_000 / sz).max(10);
        let t_dm0 = Instant::now();
        for _ in 0..iters {
            for i in 0..sz {
                let mut sum = 0.0;
                for j in 0..sz {
                    sum += sub_mat.read(i, j) * x_vec[j];
                }
                y_dense[i] = sum;
            }
        }
        let dense_matvec_time_us = (t_dm0.elapsed().as_secs_f64() / iters as f64) * 1e6;
        let dense_matvec_flops = (2 * sz * sz) as u64;

        let mut hmat_flops = 0u64;
        let t_hm0 = Instant::now();
        for _ in 0..iters {
            hmat_flops = 0;
            y_hmat.fill(0.0);
            hmatrix_matvec(&tree_r4, &x_vec, &mut y_hmat, &mut hmat_flops);
        }
        let hmat_matvec_time_us = (t_hm0.elapsed().as_secs_f64() / iters as f64) * 1e6;
        let speedup = dense_matvec_time_us / hmat_matvec_time_us.max(1e-6);

        println!(
            "   N = {:4}: Full SVD = {:6.2} GFLOP/s ({:4.1}% peak) | H-Mat r=4 = {:6.2} GFLOP/s ({:4.1}% peak) | Matvec speedup: {:5.2}x",
            sz, gflops_svd, peak_svd, gflops_r4, peak_r4, speedup
        );

        scaling_results.push(ScalingPoint {
            n: sz,
            full_svd_time_s: t_svd,
            full_svd_flops: svd_flops,
            full_svd_gflops: gflops_svd,
            full_svd_peak_percent: peak_svd,
            hmat_r1_time_s: t_r1,
            hmat_r1_flops: flops_r1,
            hmat_r1_gflops: gflops_r1,
            hmat_r1_peak_percent: peak_r1,
            hmat_r1_stored_floats: tree_r1.stored_floats(),
            hmat_r4_time_s: t_r4,
            hmat_r4_flops: flops_r4,
            hmat_r4_gflops: gflops_r4,
            hmat_r4_peak_percent: peak_r4,
            hmat_r4_stored_floats: tree_r4.stored_floats(),
            dense_matvec_time_us,
            dense_matvec_flops,
            hmat_matvec_time_us,
            hmat_matvec_flops: hmat_flops,
            matvec_speedup: speedup,
        });
    }

    // 5. Small 4x4 and 16x16 Test Verification
    println!("Generating small 4x4 and 16x16 verification datasets for NumPy / GNU Octave...");
    let mut small_tests = Vec::new();

    // 4x4 matrix from hand example:
    let m4_raw = vec![
        vec![4.0, 4.0, 1.0, 1.0],
        vec![4.0, 4.0, 1.0, 1.0],
        vec![2.0, 2.0, 8.0, 8.0],
        vec![2.0, 2.0, 8.0, 8.0],
    ];
    let m4 = Mat::from_fn(4, 4, |i, j| m4_raw[i][j]);
    let svd_4 = compute_thin_svd(&m4);
    let s4_diag = svd_4.s_diagonal();
    let mut s4_vals = Vec::new();
    for i in 0..4 {
        s4_vals.push(s4_diag.read(i));
    }

    let mut flops4 = 0u64;
    let tree4 = create_tree(&m4, 0, 4, 0, 4, 1, 0.5, &mut flops4);
    let mut rec4 = Mat::<f64>::zeros(4, 4);
    let mut dflops4 = 0u64;
    decompress_hmatrix(&tree4, &mut rec4, &mut dflops4);

    let (frob4, _, rel4) = frobenius_error(&m4, &rec4);
    let (z4, l4, _) = tree4.node_counts();
    let mut rec4_raw = vec![vec![0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            rec4_raw[i][j] = rec4.read(i, j);
        }
    }

    small_tests.push(SmallTestMatrixData {
        name: "hand_example_4x4".to_string(),
        rows: 4,
        cols: 4,
        matrix: m4_raw,
        singular_values: s4_vals,
        r: 1,
        delta: 0.5,
        reconstructed: rec4_raw,
        frob_diff: frob4,
        rel_frob_diff: rel4,
        leaves_rank: l4,
        leaves_zero: z4,
    });

    // 16x16 matrix slice from Red channel
    let mut m16_raw = vec![vec![0.0; 16]; 16];
    let m16 = Mat::from_fn(16, 16, |i, j| {
        let v = r_mat.read(i, j);
        m16_raw[i][j] = v;
        v
    });
    let svd_16 = compute_thin_svd(&m16);
    let s16_diag = svd_16.s_diagonal();
    let mut s16_vals = Vec::new();
    for i in 0..16 {
        s16_vals.push(s16_diag.read(i));
    }
    let mut flops16 = 0u64;
    let tree16 = create_tree(&m16, 0, 16, 0, 16, 2, s16_vals[8], &mut flops16);
    let mut rec16 = Mat::<f64>::zeros(16, 16);
    let mut dflops16 = 0u64;
    decompress_hmatrix(&tree16, &mut rec16, &mut dflops16);
    let (frob16, _, rel16) = frobenius_error(&m16, &rec16);
    let (z16, l16, _) = tree16.node_counts();
    let mut rec16_raw = vec![vec![0.0; 16]; 16];
    for i in 0..16 {
        for j in 0..16 {
            rec16_raw[i][j] = rec16.read(i, j);
        }
    }

    small_tests.push(SmallTestMatrixData {
        name: "image_slice_16x16".to_string(),
        rows: 16,
        cols: 16,
        matrix: m16_raw,
        singular_values: s16_vals,
        r: 2,
        delta: s16_diag.read(8),
        reconstructed: rec16_raw,
        frob_diff: frob16,
        rel_frob_diff: rel16,
        leaves_rank: l16,
        leaves_zero: z16,
    });

    // 6. Save JSON and CSV outputs via io_stats
    println!("Writing results to results/results.json and results/small_test.json...");
    let full_results = FullExperimentResults {
        image_dim: n,
        singular_values_r: sv_data_r,
        singular_values_g: sv_data_g,
        singular_values_b: sv_data_b,
        cases: case_results.clone(),
        memory_sweep: sweep_results,
        scaling: scaling_results.clone(),
        cpu_model,
        peak_gflops,
    };

    save_results_json("results/results.json", &full_results)?;
    save_small_tests_json("results/small_test.json", &small_tests)?;

    // Also export CSV summaries
    save_cases_csv("results/cases.csv", &case_results)?;
    save_scaling_csv("results/scaling.csv", &scaling_results)?;

    println!("=== All Rust Computations and Visualizations Completed Successfully! ===");
    Ok(())
}
