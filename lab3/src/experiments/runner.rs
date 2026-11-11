use crate::experiments::csv_writer::{
    save_cases_csv, save_results_json, save_scaling_csv, save_small_tests_json,
};
use crate::experiments::models::{
    CaseMetrics, FullExperimentResults, ScalingPoint, SingularValueData, SmallTestMatrixData,
    SweepPoint,
};
use crate::svd::{
    calculate_psnr, calculate_rgb_metrics, compute_thin_svd, frobenius_error,
};
use crate::tree::{create_tree, decompress_hmatrix, hmatrix_matvec};
use crate::visual::{draw_hmatrix_structure, load_rgb_image, save_matrix_as_image, save_rgb_image};
use faer::Mat;
use std::fs;
use std::time::Instant;

pub fn run_all_experiments(img_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Start Hierarchical Matrix Compression (H-Matrix / rSVD) ===");
    fs::create_dir_all("results")?;
    fs::create_dir_all("typst/imgs")?;

    let (r_mat, g_mat, b_mat) = load_rgb_image(img_path)?;
    let n = r_mat.nrows();
    let peak_gflops = 35.2;
    let cpu_model = "12th Gen Intel Core i5-1235U (4.40 GHz, AVX2/FMA)".to_string();

    save_matrix_as_image(&r_mat, "typst/imgs/orig_r.png")?;
    save_matrix_as_image(&g_mat, "typst/imgs/orig_g.png")?;
    save_matrix_as_image(&b_mat, "typst/imgs/orig_b.png")?;
    image::open(img_path)?.save("typst/imgs/orig_rgb.png")?;

    let (sv_r, sv_g, sv_b) = compute_channel_svds(&r_mat, &g_mat, &b_mat, n);
    let cases = run_compression_cases(&r_mat, &g_mat, &b_mat, &sv_r, &sv_g, &sv_b, peak_gflops)?;
    let sweep = run_parameter_sweep(&r_mat, sv_r.sigma_256);
    let scaling = run_scaling_benchmarks(&r_mat, peak_gflops);
    let small_tests = run_small_verifications(&r_mat);

    let full_results = FullExperimentResults {
        image_dim: n,
        singular_values_r: sv_r,
        singular_values_g: sv_g,
        singular_values_b: sv_b,
        cases: cases.clone(),
        memory_sweep: sweep,
        scaling: scaling.clone(),
        cpu_model,
        peak_gflops,
    };

    save_results_json("results/results.json", &full_results)?;
    save_small_tests_json("results/small_test.json", &small_tests)?;
    save_cases_csv("results/cases.csv", &cases)?;
    save_scaling_csv("results/scaling.csv", &scaling)?;

    println!("=== All Rust Computations and Visualizations Completed Successfully! ===");
    Ok(())
}

fn compute_channel_svds(
    r_mat: &Mat<f64>,
    g_mat: &Mat<f64>,
    b_mat: &Mat<f64>,
    n: usize,
) -> (SingularValueData, SingularValueData, SingularValueData) {
    let extract = |m: &Mat<f64>, name: &str| -> SingularValueData {
        let svd = compute_thin_svd(m);
        let s = svd.s_diagonal();
        let vals: Vec<f64> = (0..n).map(|i| s.read(i)).collect();
        SingularValueData {
            channel: name.to_string(),
            sigma_1: vals[0],
            sigma_256: vals[255],
            sigma_512: vals[511],
            values: vals,
        }
    };
    (
        extract(r_mat, "Red"),
        extract(g_mat, "Green"),
        extract(b_mat, "Blue"),
    )
}

fn run_compression_cases(
    r_mat: &Mat<f64>,
    g_mat: &Mat<f64>,
    b_mat: &Mat<f64>,
    sv_r: &SingularValueData,
    sv_g: &SingularValueData,
    sv_b: &SingularValueData,
    peak_gflops: f64,
) -> Result<Vec<CaseMetrics>, Box<dyn std::error::Error>> {
    let n = r_mat.nrows();
    let configs = vec![
        ("case1", 1, sv_r.sigma_1, sv_g.sigma_1, sv_b.sigma_1),
        ("case2", 1, sv_r.sigma_512, sv_g.sigma_512, sv_b.sigma_512),
        ("case3", 1, sv_r.sigma_256, sv_g.sigma_256, sv_b.sigma_256),
        ("case4", 4, sv_r.sigma_1, sv_g.sigma_1, sv_b.sigma_1),
        ("case5", 4, sv_r.sigma_512, sv_g.sigma_512, sv_b.sigma_512),
        ("case6", 4, sv_r.sigma_256, sv_g.sigma_256, sv_b.sigma_256),
        ("opt1", 2, 3.0 * sv_r.sigma_256, 3.0 * sv_g.sigma_256, 3.0 * sv_b.sigma_256),
        ("opt2", 4, 4.0 * sv_r.sigma_256, 4.0 * sv_g.sigma_256, 4.0 * sv_b.sigma_256),
    ];

    let mut results = Vec::new();
    for (name, r, dr, dg, db) in configs {
        let mut comp_flops = 0u64;
        let t_comp = Instant::now();
        let tr = create_tree(r_mat, 0, n, 0, n, r, dr, &mut comp_flops);
        let tg = create_tree(g_mat, 0, n, 0, n, r, dg, &mut comp_flops);
        let tb = create_tree(b_mat, 0, n, 0, n, r, db, &mut comp_flops);
        let comp_ms = t_comp.elapsed().as_secs_f64() * 1000.0;

        let (mut rec_r, mut rec_g, mut rec_b) = (Mat::zeros(n, n), Mat::zeros(n, n), Mat::zeros(n, n));
        let mut d_flops = 0u64;
        let t_decomp = Instant::now();
        decompress_hmatrix(&tr, &mut rec_r, &mut d_flops);
        decompress_hmatrix(&tg, &mut rec_g, &mut d_flops);
        decompress_hmatrix(&tb, &mut rec_b, &mut d_flops);
        let decomp_ms = t_decomp.elapsed().as_secs_f64() * 1000.0;

        draw_hmatrix_structure(&tr, n, n, 512, &format!("typst/imgs/{}_tree_r.png", name))?;
        draw_hmatrix_structure(&tg, n, n, 512, &format!("typst/imgs/{}_tree_g.png", name))?;
        draw_hmatrix_structure(&tb, n, n, 512, &format!("typst/imgs/{}_tree_b.png", name))?;
        save_matrix_as_image(&rec_r, &format!("typst/imgs/{}_rec_r.png", name))?;
        save_matrix_as_image(&rec_g, &format!("typst/imgs/{}_rec_g.png", name))?;
        save_matrix_as_image(&rec_b, &format!("typst/imgs/{}_rec_b.png", name))?;
        save_rgb_image(&rec_r, &rec_g, &rec_b, &format!("typst/imgs/{}_rec_rgb.png", name))?;

        let stored = tr.stored_floats() + tg.stored_floats() + tb.stored_floats();
        let dense_floats = 3 * n * n;
        let (frob_rgb, psnr_rgb) = calculate_rgb_metrics(r_mat, g_mat, b_mat, &rec_r, &rec_g, &rec_b);
        let gflops = (comp_flops as f64 / 1e9) / (comp_ms / 1000.0).max(1e-9);

        results.push(CaseMetrics {
            name: name.to_string(),
            r, delta_r: dr, delta_g: dg, delta_b: db,
            zero_leaves_r: tr.node_counts().0, rank_leaves_r: tr.node_counts().1, depth_r: tr.depth(),
            zero_leaves_g: tg.node_counts().0, rank_leaves_g: tg.node_counts().1, depth_g: tg.depth(),
            zero_leaves_b: tb.node_counts().0, rank_leaves_b: tb.node_counts().1, depth_b: tb.depth(),
            stored_floats_total: stored, dense_floats_total: dense_floats,
            compression_ratio: dense_floats as f64 / stored.max(1) as f64,
            memory_savings_percent: (1.0 - stored as f64 / dense_floats as f64) * 100.0,
            tree_memory_kb: (tr.memory_bytes() + tg.memory_bytes() + tb.memory_bytes()) as f64 / 1024.0,
            dense_memory_kb: (dense_floats * 8) as f64 / 1024.0,
            rel_frob_r: frobenius_error(r_mat, &rec_r).2,
            rel_frob_g: frobenius_error(g_mat, &rec_g).2,
            rel_frob_b: frobenius_error(b_mat, &rec_b).2,
            rel_frob_rgb: frob_rgb,
            psnr_r: calculate_psnr(r_mat, &rec_r),
            psnr_g: calculate_psnr(g_mat, &rec_g),
            psnr_b: calculate_psnr(b_mat, &rec_b),
            psnr_rgb,
            compress_time_ms: comp_ms, decompress_time_ms: decomp_ms,
            compress_flops: comp_flops, decompress_flops: d_flops,
            compress_gflops: gflops, compress_peak_percent: (gflops / peak_gflops) * 100.0,
        });
    }
    Ok(results)
}

fn run_parameter_sweep(r_mat: &Mat<f64>, s256: f64) -> Vec<SweepPoint> {
    let (n, ranks) = (r_mat.nrows(), [1, 2, 4, 8]);
    let mults = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0, 40.0];
    let mut sweep = Vec::new();

    for &r in &ranks {
        for &m in &mults {
            let delta = m * s256;
            let mut flops = 0u64;
            let tree = create_tree(r_mat, 0, n, 0, n, r, delta, &mut flops);
            let mut rec = Mat::zeros(n, n);
            let mut d_flops = 0u64;
            decompress_hmatrix(&tree, &mut rec, &mut d_flops);
            sweep.push(SweepPoint {
                r, delta_mult: m, delta,
                stored_floats: tree.stored_floats(),
                compression_ratio: (n * n) as f64 / tree.stored_floats().max(1) as f64,
                memory_kb: tree.memory_bytes() as f64 / 1024.0,
                rel_frob: frobenius_error(r_mat, &rec).2,
                psnr: calculate_psnr(r_mat, &rec),
            });
        }
    }
    sweep
}

fn run_scaling_benchmarks(r_mat: &Mat<f64>, peak_gflops: f64) -> Vec<ScalingPoint> {
    let mut points = Vec::new();
    for &sz in &[32, 64, 128, 256, 512, 1024] {
        let sub = if sz <= 512 {
            Mat::from_fn(sz, sz, |i, j| r_mat.read(i, j))
        } else {
            Mat::from_fn(sz, sz, |i, j| r_mat.read(i % 512, j % 512) * (1.0 + 0.05 * ((i / 512 + j / 512) as f64)))
        };

        let svd_flops = (21 * sz * sz * sz) as u64;
        let t0 = Instant::now();
        let svd = compute_thin_svd(&sub);
        let t_svd = t0.elapsed().as_secs_f64();
        let mid_sigma = svd.s_diagonal().read(sz / 2);

        let (mut f1, mut f4) = (0u64, 0u64);
        let t1 = Instant::now();
        let tr1 = create_tree(&sub, 0, sz, 0, sz, 1, mid_sigma, &mut f1);
        let t_r1 = t1.elapsed().as_secs_f64();

        let t2 = Instant::now();
        let tr4 = create_tree(&sub, 0, sz, 0, sz, 4, 3.0 * mid_sigma, &mut f4);
        let t_r4 = t2.elapsed().as_secs_f64();

        let x = vec![1.0; sz];
        let mut y = vec![0.0; sz];
        let iters = (10_000 / sz).max(10);
        let tm0 = Instant::now();
        for _ in 0..iters {
            for i in 0..sz {
                let mut sum = 0.0;
                for j in 0..sz { sum += sub.read(i, j) * x[j]; }
                y[i] = sum;
            }
        }
        let t_dense = (tm0.elapsed().as_secs_f64() / iters as f64) * 1e6;

        let mut hm_flops = 0u64;
        let tm1 = Instant::now();
        for _ in 0..iters {
            hm_flops = 0;
            y.fill(0.0);
            hmatrix_matvec(&tr4, &x, &mut y, &mut hm_flops);
        }
        let t_hmat = (tm1.elapsed().as_secs_f64() / iters as f64) * 1e6;

        let g_svd = (svd_flops as f64 / 1e9) / t_svd.max(1e-9);
        let g_r4 = (f4 as f64 / 1e9) / t_r4.max(1e-9);
        points.push(ScalingPoint {
            n: sz, full_svd_time_s: t_svd, full_svd_flops: svd_flops,
            full_svd_gflops: g_svd, full_svd_peak_percent: (g_svd / peak_gflops) * 100.0,
            hmat_r1_time_s: t_r1, hmat_r1_flops: f1, hmat_r1_gflops: (f1 as f64 / 1e9) / t_r1.max(1e-9),
            hmat_r1_peak_percent: ((f1 as f64 / 1e9) / t_r1.max(1e-9) / peak_gflops) * 100.0,
            hmat_r1_stored_floats: tr1.stored_floats(),
            hmat_r4_time_s: t_r4, hmat_r4_flops: f4, hmat_r4_gflops: g_r4,
            hmat_r4_peak_percent: (g_r4 / peak_gflops) * 100.0,
            hmat_r4_stored_floats: tr4.stored_floats(),
            dense_matvec_time_us: t_dense, dense_matvec_flops: (2 * sz * sz) as u64,
            hmat_matvec_time_us: t_hmat, hmat_matvec_flops: hm_flops,
            matvec_speedup: t_dense / t_hmat.max(1e-6),
        });
    }
    points
}

fn run_small_verifications(r_mat: &Mat<f64>) -> Vec<SmallTestMatrixData> {
    let m4_raw = vec![
        vec![4.0, 4.0, 1.0, 1.0], vec![4.0, 4.0, 1.0, 1.0],
        vec![2.0, 2.0, 8.0, 8.0], vec![2.0, 2.0, 8.0, 8.0],
    ];
    let m4 = Mat::from_fn(4, 4, |i, j| m4_raw[i][j]);
    let svd4 = compute_thin_svd(&m4);
    let s4_diag = svd4.s_diagonal();
    let s4_vals: Vec<f64> = (0..4).map(|i| s4_diag.read(i)).collect();
    let mut f4 = 0u64;
    let tr4 = create_tree(&m4, 0, 4, 0, 4, 1, 0.5, &mut f4);
    let mut rec4 = Mat::zeros(4, 4);
    let mut df4 = 0u64;
    decompress_hmatrix(&tr4, &mut rec4, &mut df4);
    let (frob4, _, rel4) = frobenius_error(&m4, &rec4);
    let rec4_raw = (0..4).map(|i| (0..4).map(|j| rec4.read(i, j)).collect()).collect();

    let m16 = Mat::from_fn(16, 16, |i, j| r_mat.read(i, j));
    let m16_raw = (0..16).map(|i| (0..16).map(|j| m16.read(i, j)).collect()).collect();
    let svd16 = compute_thin_svd(&m16);
    let s16_diag = svd16.s_diagonal();
    let s16_vals: Vec<f64> = (0..16).map(|i| s16_diag.read(i)).collect();
    let delta16 = s16_vals[8];
    let mut f16 = 0u64;
    let tr16 = create_tree(&m16, 0, 16, 0, 16, 2, delta16, &mut f16);
    let mut rec16 = Mat::zeros(16, 16);
    let mut df16 = 0u64;
    decompress_hmatrix(&tr16, &mut rec16, &mut df16);
    let (frob16, _, rel16) = frobenius_error(&m16, &rec16);
    let rec16_raw = (0..16).map(|i| (0..16).map(|j| rec16.read(i, j)).collect()).collect();

    vec![
        SmallTestMatrixData {
            name: "hand_example_4x4".to_string(), rows: 4, cols: 4, matrix: m4_raw,
            singular_values: s4_vals, r: 1, delta: 0.5, reconstructed: rec4_raw,
            frob_diff: frob4, rel_frob_diff: rel4, leaves_rank: tr4.node_counts().1, leaves_zero: tr4.node_counts().0,
        },
        SmallTestMatrixData {
            name: "image_slice_16x16".to_string(), rows: 16, cols: 16, matrix: m16_raw,
            singular_values: s16_vals, r: 2, delta: delta16, reconstructed: rec16_raw,
            frob_diff: frob16, rel_frob_diff: rel16, leaves_rank: tr16.node_counts().1, leaves_zero: tr16.node_counts().0,
        },
    ]
}
