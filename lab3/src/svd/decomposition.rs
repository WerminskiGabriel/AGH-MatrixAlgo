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

/// Thin SVD wrapper using faer
#[inline]
pub fn compute_thin_svd(mat: &Mat<f64>) -> ThinSvd<f64> {
    mat.thin_svd()
}
