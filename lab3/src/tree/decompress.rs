use crate::tree::node::{HNode, HNodeKind};
use faer::Mat;

pub fn decompress_hmatrix(node: &HNode, out: &mut Mat<f64>, flops: &mut u64) {
    match &node.kind {
        HNodeKind::Zero => {
            for i in node.t_min..node.t_max {
                for j in node.s_min..node.s_max {
                    out.write(i, j, 0.0);
                }
            }
        }
        HNodeKind::Leaf { rank, u, v, .. } => {
            let (m, n, k) = (node.nrows(), node.ncols(), *rank);
            for i in 0..m {
                for j in 0..n {
                    let mut sum = 0.0;
                    for p in 0..k {
                        sum += u[i * k + p] * v[p * n + j];
                    }
                    out.write(node.t_min + i, node.s_min + j, sum);
                }
            }
            if k > 0 {
                *flops += (m * n * (2 * k - 1)) as u64;
            }
        }
        HNodeKind::Internal { children } => {
            for child in children.iter() {
                decompress_hmatrix(child, out, flops);
            }
        }
    }
}
