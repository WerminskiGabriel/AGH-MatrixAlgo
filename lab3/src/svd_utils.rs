use faer::linalg::solvers::ThinSvd;
use faer::Mat;

/// Estimates FLOPs count for thin SVD of m x n matrix (Golub-Van Loan: 4*m^2*n + 8*m*n^2 + 9*n^3)
pub fn estimate_svd_flops(m: usize, n: usize) -> u64 {
    let (dim_max, dim_min) = if m >= n {
        (m as u64, n as u64)
    } else {
        (n as u64, m as u64)
    };
    4 * dim_max * dim_max * dim_min
        + 8 * dim_max * dim_min * dim_min
        + 9 * dim_min * dim_min * dim_min
}

/// Compute thin SVD wrapper using faer
pub fn compute_thin_svd(mat: &Mat<f64>) -> ThinSvd<f64> {
    mat.thin_svd()
}

/// Compute Frobenius norm of difference ||A - B||_F, ||A||_F and relative error
pub fn frobenius_error(orig: &Mat<f64>, rec: &Mat<f64>) -> (f64, f64, f64) {
    let rows = orig.nrows();
    let cols = orig.ncols();
    let mut diff_sum_sq = 0.0;
    let mut orig_sum_sq = 0.0;

    for i in 0..rows {
        for j in 0..cols {
            let o = orig.read(i, j);
            let r = rec.read(i, j);
            let d = o - r;
            diff_sum_sq += d * d;
            orig_sum_sq += o * o;
        }
    }

    let diff_frob = diff_sum_sq.sqrt();
    let orig_frob = orig_sum_sq.sqrt();
    let rel_frob = if orig_frob > 0.0 {
        diff_frob / orig_frob
    } else {
        0.0
    };
    (diff_frob, orig_frob, rel_frob)
}

/// Compute PSNR [dB] between original and reconstructed matrix
pub fn calculate_psnr(orig: &Mat<f64>, rec: &Mat<f64>) -> f64 {
    let rows = orig.nrows();
    let cols = orig.ncols();
    let mut mse = 0.0;
    for i in 0..rows {
        for j in 0..cols {
            let d = orig.read(i, j) - rec.read(i, j);
            mse += d * d;
        }
    }
    mse /= (rows * cols) as f64;
    if mse < 1e-15 {
        return 99.0;
    }
    10.0 * (255.0 * 255.0 / mse).log10()
}

/// Compute combined RGB relative Frobenius error and PSNR
pub fn calculate_rgb_metrics(
    r_orig: &Mat<f64>,
    g_orig: &Mat<f64>,
    b_orig: &Mat<f64>,
    r_rec: &Mat<f64>,
    g_rec: &Mat<f64>,
    b_rec: &Mat<f64>,
) -> (f64, f64) {
    let n = r_orig.nrows();
    let mut total_diff_sq = 0.0;
    let mut total_orig_sq = 0.0;
    for i in 0..n {
        for j in 0..n {
            let dr = r_orig.read(i, j) - r_rec.read(i, j);
            let dg = g_orig.read(i, j) - g_rec.read(i, j);
            let db = b_orig.read(i, j) - b_rec.read(i, j);
            total_diff_sq += dr * dr + dg * dg + db * db;
            let or = r_orig.read(i, j);
            let og = g_orig.read(i, j);
            let ob = b_orig.read(i, j);
            total_orig_sq += or * or + og * og + ob * ob;
        }
    }
    let rel_frob_rgb = (total_diff_sq / total_orig_sq.max(1e-15)).sqrt();
    let mse_rgb = total_diff_sq / (3.0 * (n * n) as f64);
    let psnr_rgb = if mse_rgb < 1e-15 {
        99.0
    } else {
        10.0 * (255.0 * 255.0 / mse_rgb).log10()
    };
    (rel_frob_rgb, psnr_rgb)
}

/// Calculate compression ratio and savings percentage
pub fn calculate_compression_stats(stored_floats: usize, dense_floats: usize) -> (f64, f64) {
    let cr = dense_floats as f64 / stored_floats.max(1) as f64;
    let savings = (1.0 - (stored_floats as f64 / dense_floats as f64)) * 100.0;
    (cr, savings)
}
