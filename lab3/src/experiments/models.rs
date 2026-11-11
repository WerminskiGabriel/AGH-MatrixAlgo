use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CaseMetrics {
    pub name: String,
    pub r: usize,
    pub delta_r: f64,
    pub delta_g: f64,
    pub delta_b: f64,
    pub zero_leaves_r: usize,
    pub rank_leaves_r: usize,
    pub depth_r: usize,
    pub zero_leaves_g: usize,
    pub rank_leaves_g: usize,
    pub depth_g: usize,
    pub zero_leaves_b: usize,
    pub rank_leaves_b: usize,
    pub depth_b: usize,
    pub stored_floats_total: usize,
    pub dense_floats_total: usize,
    pub compression_ratio: f64,
    pub memory_savings_percent: f64,
    pub tree_memory_kb: f64,
    pub dense_memory_kb: f64,
    pub rel_frob_r: f64,
    pub rel_frob_g: f64,
    pub rel_frob_b: f64,
    pub rel_frob_rgb: f64,
    pub psnr_r: f64,
    pub psnr_g: f64,
    pub psnr_b: f64,
    pub psnr_rgb: f64,
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
    pub full_svd_time_s: f64,
    pub full_svd_flops: u64,
    pub full_svd_gflops: f64,
    pub full_svd_peak_percent: f64,
    pub hmat_r1_time_s: f64,
    pub hmat_r1_flops: u64,
    pub hmat_r1_gflops: f64,
    pub hmat_r1_peak_percent: f64,
    pub hmat_r1_stored_floats: usize,
    pub hmat_r4_time_s: f64,
    pub hmat_r4_flops: u64,
    pub hmat_r4_gflops: f64,
    pub hmat_r4_peak_percent: f64,
    pub hmat_r4_stored_floats: usize,
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
