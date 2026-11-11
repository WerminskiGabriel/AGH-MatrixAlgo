use crate::tree::node::{HNode, HNodeKind};
use faer::Mat;

/// Rapid matrix-vector multiplication y += H * x
pub fn hmatrix_matvec(node: &HNode, x: &[f64], y: &mut [f64], flops: &mut u64) {
    match &node.kind {
        HNodeKind::Zero => {}
        HNodeKind::Leaf { rank, u, v, .. } => {
            let (m, n, k) = (node.nrows(), node.ncols(), *rank);
            if k == 0 {
                return;
            }
            // Step 1: z = V * x[s_min..s_max] of length k
            let mut z = vec![0.0; k];
            for p in 0..k {
                let mut sum = 0.0;
                for j in 0..n {
                    sum += v[p * n + j] * x[node.s_min + j];
                }
                z[p] = sum;
            }
            *flops += (k * (2 * n - 1)) as u64;

            // Step 2: y[t_min..t_max] += U * z
            for i in 0..m {
                let mut sum = 0.0;
                for p in 0..k {
                    sum += u[i * k + p] * z[p];
                }
                y[node.t_min + i] += sum;
            }
            *flops += (m * (2 * k)) as u64;
        }
        HNodeKind::Internal { children } => {
            for child in children.iter() {
                hmatrix_matvec(child, x, y, flops);
            }
        }
    }
}

/// Dense matrix-vector multiplication y = A * x for benchmark comparison
pub fn dense_matvec(a: &Mat<f64>, x: &[f64], y: &mut [f64]) {
    let (m, n) = (a.nrows(), a.ncols());
    for i in 0..m {
        let mut sum = 0.0;
        for j in 0..n {
            sum += a.read(i, j) * x[j];
        }
        y[i] = sum;
    }
}
