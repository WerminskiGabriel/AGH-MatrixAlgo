use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CaseMetrics {
    pub name: String,
    pub r: usize,
    pub delta_r: f64,
    pub delta_g: f64,
    pub delta_b: f64,
    // Leaves and depth
    pub zero_leaves_r: usize,
    pub rank_leaves_r: usize,
    pub depth_r: usize,
    pub zero_leaves_g: usize,
    pub rank_leaves_g: usize,
    pub depth_g: usize,
    pub zero_leaves_b: usize,
    pub rank_leaves_b: usize,
    pub depth_b: usize,
    // Stored floats and memory
    pub stored_floats_total: usize,
    pub dense_floats_total: usize,
    pub compression_ratio: f64,
    pub memory_savings_percent: f64,
    pub tree_memory_kb: f64,
    pub dense_memory_kb: f64,
    // Reconstruction quality
    pub rel_frob_r: f64,
    pub rel_frob_g: f64,
    pub rel_frob_b: f64,
    pub rel_frob_rgb: f64,
    pub psnr_r: f64,
    pub psnr_g: f64,
    pub psnr_b: f64,
    pub psnr_rgb: f64,
    // Performance
    pub compress_time_ms: f64,
    pub decompress_time_ms: f64,
    pub compress_flops: u64,
    pub decompress_flops: u64,
    pub compress_gflops: f64,
    pub compress_peak_percent: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SingularValueData {
    pub channel: String,
    pub values: Vec<f64>,
    pub sigma_1: f64,
    pub sigma_256: f64,
    pub sigma_512: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SweepPoint {
    pub r: usize,
    pub delta_mult: f64,
    pub delta: f64,
    pub stored_floats: usize,
    pub compression_ratio: f64,
    pub memory_kb: f64,
    pub rel_frob: f64,
    pub psnr: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ScalingPoint {
    pub n: usize,
    // Full SVD
    pub full_svd_time_s: f64,
    pub full_svd_flops: u64,
    pub full_svd_gflops: f64,
    pub full_svd_peak_percent: f64,
    // H-matrix compress r=1
    pub hmat_r1_time_s: f64,
    pub hmat_r1_flops: u64,
    pub hmat_r1_gflops: f64,
    pub hmat_r1_peak_percent: f64,
    pub hmat_r1_stored_floats: usize,
    // H-matrix compress r=4
    pub hmat_r4_time_s: f64,
    pub hmat_r4_flops: u64,
    pub hmat_r4_gflops: f64,
    pub hmat_r4_peak_percent: f64,
    pub hmat_r4_stored_floats: usize,
    // Matvec comparison
    pub dense_matvec_time_us: f64,
    pub dense_matvec_flops: u64,
    pub hmat_matvec_time_us: f64,
    pub hmat_matvec_flops: u64,
    pub matvec_speedup: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct FullExperimentResults {
    pub image_dim: usize,
    pub singular_values_r: SingularValueData,
    pub singular_values_g: SingularValueData,
    pub singular_values_b: SingularValueData,
    pub cases: Vec<CaseMetrics>,
    pub memory_sweep: Vec<SweepPoint>,
    pub scaling: Vec<ScalingPoint>,
    pub cpu_model: String,
    pub peak_gflops: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SmallTestMatrixData {
    pub name: String,
    pub rows: usize,
    pub cols: usize,
    pub matrix: Vec<Vec<f64>>,
    pub singular_values: Vec<f64>,
    pub r: usize,
    pub delta: f64,
    pub reconstructed: Vec<Vec<f64>>,
    pub frob_diff: f64,
    pub rel_frob_diff: f64,
    pub leaves_rank: usize,
    pub leaves_zero: usize,
}

/// Save complete experiment results to formatted JSON file
pub fn save_results_json(
    path: &str,
    results: &FullExperimentResults,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    serde_json::to_writer_pretty(&mut file, results)?;
    Ok(())
}

/// Save small test verification matrices to JSON file
pub fn save_small_tests_json(
    path: &str,
    tests: &[SmallTestMatrixData],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    serde_json::to_writer_pretty(&mut file, tests)?;
    Ok(())
}

/// Save summary of compression cases to CSV format
pub fn save_cases_csv(path: &str, cases: &[CaseMetrics]) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "name,r,delta_r,delta_g,delta_b,compression_ratio,memory_savings_percent,tree_memory_kb,dense_memory_kb,rel_frob_rgb,psnr_rgb,compress_time_ms,decompress_time_ms,compress_flops,compress_gflops,compress_peak_percent"
    )?;
    for c in cases {
        writeln!(
            file,
            "{},{},{:.4},{:.4},{:.4},{:.4},{:.2},{:.2},{:.2},{:.6},{:.2},{:.2},{:.2},{},{:.2},{:.1}",
            c.name,
            c.r,
            c.delta_r,
            c.delta_g,
            c.delta_b,
            c.compression_ratio,
            c.memory_savings_percent,
            c.tree_memory_kb,
            c.dense_memory_kb,
            c.rel_frob_rgb,
            c.psnr_rgb,
            c.compress_time_ms,
            c.decompress_time_ms,
            c.compress_flops,
            c.compress_gflops,
            c.compress_peak_percent
        )?;
    }
    Ok(())
}

/// Save scaling benchmark results to CSV format
pub fn save_scaling_csv(
    path: &str,
    scaling: &[ScalingPoint],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "n,full_svd_time_s,full_svd_flops,full_svd_gflops,full_svd_peak_percent,hmat_r4_time_s,hmat_r4_flops,hmat_r4_gflops,hmat_r4_peak_percent,matvec_speedup"
    )?;
    for s in scaling {
        writeln!(
            file,
            "{},{:.6},{},{:.2},{:.1},{:.6},{},{:.2},{:.1},{:.2}",
            s.n,
            s.full_svd_time_s,
            s.full_svd_flops,
            s.full_svd_gflops,
            s.full_svd_peak_percent,
            s.hmat_r4_time_s,
            s.hmat_r4_flops,
            s.hmat_r4_gflops,
            s.hmat_r4_peak_percent,
            s.matvec_speedup
        )?;
    }
    Ok(())
}
