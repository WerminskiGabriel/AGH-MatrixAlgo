use crate::experiments::models::{CaseMetrics, FullExperimentResults, ScalingPoint, SmallTestMatrixData};
use std::fs::File;
use std::io::Write;

pub fn save_results_json(
    path: &str,
    results: &FullExperimentResults,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    serde_json::to_writer_pretty(&mut file, results)?;
    Ok(())
}

pub fn save_small_tests_json(
    path: &str,
    tests: &[SmallTestMatrixData],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    serde_json::to_writer_pretty(&mut file, tests)?;
    Ok(())
}

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
            c.name, c.r, c.delta_r, c.delta_g, c.delta_b, c.compression_ratio,
            c.memory_savings_percent, c.tree_memory_kb, c.dense_memory_kb,
            c.rel_frob_rgb, c.psnr_rgb, c.compress_time_ms, c.decompress_time_ms,
            c.compress_flops, c.compress_gflops, c.compress_peak_percent
        )?;
    }
    Ok(())
}

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
            s.n, s.full_svd_time_s, s.full_svd_flops, s.full_svd_gflops,
            s.full_svd_peak_percent, s.hmat_r4_time_s, s.hmat_r4_flops,
            s.hmat_r4_gflops, s.hmat_r4_peak_percent, s.matvec_speedup
        )?;
    }
    Ok(())
}
