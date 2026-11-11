use crate::svd::estimate_svd_flops;
use crate::tree::node::{HNode, HNodeKind};
use faer::linalg::solvers::ThinSvd;
use faer::Mat;

pub fn create_tree(
    a: &Mat<f64>,
    t_min: usize,
    t_max: usize,
    s_min: usize,
    s_max: usize,
    r: usize,
    delta: f64,
    flops: &mut u64,
) -> HNode {
    let m = t_max - t_min;
    let n = s_max - s_min;
    assert!(m > 0 && n > 0);

    let mut max_val: f64 = 0.0;
    for i in t_min..t_max {
        for j in s_min..s_max {
            max_val = max_val.max(a.read(i, j).abs());
        }
    }
    if max_val < 1e-12 {
        return HNode {
            t_min,
            t_max,
            s_min,
            s_max,
            kind: HNodeKind::Zero,
        };
    }

    let block = Mat::from_fn(m, n, |i, j| a.read(t_min + i, s_min + j));
    let svd = block.thin_svd();
    let s_diag = svd.s_diagonal();
    *flops += estimate_svd_flops(m, n);

    let is_admissible = if m.min(n) <= r {
        true
    } else {
        s_diag.read(r) < delta
    };

    if is_admissible || m <= 2 || n <= 2 {
        compress_matrix(t_min, t_max, s_min, s_max, &svd, r, delta, flops)
    } else {
        let (t_mid, s_mid) = (t_min + m / 2, s_min + n / 2);
        let tl = create_tree(a, t_min, t_mid, s_min, s_mid, r, delta, flops);
        let tr = create_tree(a, t_min, t_mid, s_mid, s_max, r, delta, flops);
        let bl = create_tree(a, t_mid, t_max, s_min, s_mid, r, delta, flops);
        let br = create_tree(a, t_mid, t_max, s_mid, s_max, r, delta, flops);
        HNode {
            t_min,
            t_max,
            s_min,
            s_max,
            kind: HNodeKind::Internal {
                children: Box::new([tl, tr, bl, br]),
            },
        }
    }
}

pub fn compress_matrix(
    t_min: usize,
    t_max: usize,
    s_min: usize,
    s_max: usize,
    svd: &ThinSvd<f64>,
    r: usize,
    delta: f64,
    flops: &mut u64,
) -> HNode {
    let (m, n) = (t_max - t_min, s_max - s_min);
    let s_diag = svd.s_diagonal();
    if s_diag.read(0) < delta.min(1e-12) {
        return HNode {
            t_min,
            t_max,
            s_min,
            s_max,
            kind: HNodeKind::Zero,
        };
    }

    let k = r.min(m.min(n));
    let (u_mat, v_mat) = (svd.u(), svd.v());

    let mut u_vec = Vec::with_capacity(m * k);
    for i in 0..m {
        for p in 0..k {
            u_vec.push(u_mat.read(i, p));
        }
    }

    let mut v_vec = Vec::with_capacity(k * n);
    for p in 0..k {
        let sigma_p = s_diag.read(p);
        for j in 0..n {
            v_vec.push(sigma_p * v_mat.read(j, p));
        }
    }
    *flops += (k * n) as u64;

    let singular_values = (0..k).map(|p| s_diag.read(p)).collect();
    HNode {
        t_min,
        t_max,
        s_min,
        s_max,
        kind: HNodeKind::Leaf {
            rank: k,
            singular_values,
            u: u_vec,
            v: v_vec,
        },
    }
}
