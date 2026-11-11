use crate::svd_utils::estimate_svd_flops;
use faer::Mat;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum HNodeKind {
    /// Zero block (all values below threshold or zero)
    Zero,
    /// Low-rank compressed leaf storing U (m x rank) and V = diag(S) * V^T (rank x n)
    Leaf {
        rank: usize,
        singular_values: Vec<f64>,
        /// Row-major storage of U: dimension m * rank
        u: Vec<f64>,
        /// Row-major storage of V: dimension rank * n
        v: Vec<f64>,
    },
    /// Internal quadtree node with 4 children: [TL, TR, BL, BR]
    Internal {
        children: Box<[HNode; 4]>,
    },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HNode {
    pub t_min: usize,
    pub t_max: usize,
    pub s_min: usize,
    pub s_max: usize,
    pub kind: HNodeKind,
}

impl HNode {
    pub fn nrows(&self) -> usize {
        self.t_max - self.t_min
    }

    pub fn ncols(&self) -> usize {
        self.s_max - self.s_min
    }

    /// Recursively counts (zero_leaves, rank_leaves, internal_nodes)
    pub fn node_counts(&self) -> (usize, usize, usize) {
        match &self.kind {
            HNodeKind::Zero => (1, 0, 0),
            HNodeKind::Leaf { .. } => (0, 1, 0),
            HNodeKind::Internal { children } => {
                let mut zeros = 0;
                let mut ranks = 0;
                let mut internals = 1;
                for child in children.iter() {
                    let (z, r, i) = child.node_counts();
                    zeros += z;
                    ranks += r;
                    internals += i;
                }
                (zeros, ranks, internals)
            }
        }
    }

    /// Returns max depth of the quadtree
    pub fn depth(&self) -> usize {
        match &self.kind {
            HNodeKind::Zero | HNodeKind::Leaf { .. } => 1,
            HNodeKind::Internal { children } => {
                1 + children.iter().map(|c| c.depth()).max().unwrap_or(0)
            }
        }
    }

    /// Number of f64 numbers stored in all leaf nodes
    pub fn stored_floats(&self) -> usize {
        match &self.kind {
            HNodeKind::Zero => 0,
            HNodeKind::Leaf { rank, .. } => {
                let m = self.nrows();
                let n = self.ncols();
                *rank * (m + n) + *rank // U + V + singular values
            }
            HNodeKind::Internal { children } => {
                children.iter().map(|c| c.stored_floats()).sum()
            }
        }
    }

    /// Total memory consumption in bytes, including struct overhead and heap allocations
    pub fn memory_bytes(&self) -> usize {
        let node_overhead = std::mem::size_of::<HNode>();
        match &self.kind {
            HNodeKind::Zero => node_overhead,
            HNodeKind::Leaf {
                u,
                v,
                singular_values,
                ..
            } => {
                node_overhead
                    + u.len() * std::mem::size_of::<f64>()
                    + v.len() * std::mem::size_of::<f64>()
                    + singular_values.len() * std::mem::size_of::<f64>()
            }
            HNodeKind::Internal { children } => {
                node_overhead
                    + std::mem::size_of::<[HNode; 4]>()
                    + children.iter().map(|c| c.memory_bytes()).sum::<usize>()
            }
        }
    }
}

/// Computes thin SVD of submatrix A[t_min..t_max, s_min..s_max] and builds H-matrix tree
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

    // 1. Check if block consists of zeros (or near zero)
    let mut max_val: f64 = 0.0;
    for i in t_min..t_max {
        for j in s_min..s_max {
            let v = a.read(i, j).abs();
            if v > max_val {
                max_val = v;
            }
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

    // 2. Extract submatrix block
    let block = Mat::from_fn(m, n, |i, j| a.read(t_min + i, s_min + j));

    // Compute SVD of block using faer
    let svd = block.thin_svd();
    let min_mn = m.min(n);
    let s_diag = svd.s_diagonal();

    // Standard flop count for thin SVD of m x n matrix (Golub-Van Loan)
    *flops += estimate_svd_flops(m, n);

    // 3. Admissibility condition:
    // If min(m, n) <= r, or (r < min(m, n) and sigma_{r+1} < delta), compress to leaf!
    let is_admissible = if min_mn <= r {
        true
    } else {
        let sigma_r_plus_1 = s_diag.read(r);
        sigma_r_plus_1 < delta
    };

    if is_admissible || m <= 2 || n <= 2 {
        compress_matrix(t_min, t_max, s_min, s_max, &svd, r, delta, flops)
    } else {
        // Step 5: Internal node, subdivide into 4 quadrants
        let t_mid = t_min + m / 2;
        let s_mid = s_min + n / 2;

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

/// Compress admissible block into Leaf or Zero node (Slide 34)
pub fn compress_matrix(
    t_min: usize,
    t_max: usize,
    s_min: usize,
    s_max: usize,
    svd: &faer::linalg::solvers::ThinSvd<f64>,
    r: usize,
    delta: f64,
    flops: &mut u64,
) -> HNode {
    let m = t_max - t_min;
    let n = s_max - s_min;
    let min_mn = m.min(n);
    let s_diag = svd.s_diagonal();

    // If first singular value is already smaller than delta (or near zero), entire block is zero
    if s_diag.read(0) < delta.min(1e-12) {
        return HNode {
            t_min,
            t_max,
            s_min,
            s_max,
            kind: HNodeKind::Zero,
        };
    }

    let k = r.min(min_mn);
    let u_mat = svd.u();
    let v_mat = svd.v();

    let mut u_vec = Vec::with_capacity(m * k);
    for i in 0..m {
        for p in 0..k {
            u_vec.push(u_mat.read(i, p));
        }
    }

    // V_leaf = diag(S(0..k)) * V^T(0..k, :)
    // v_mat in faer has dimension n x min_mn, so V^T(p, j) = v_mat.read(j, p)
    let mut v_vec = Vec::with_capacity(k * n);
    for p in 0..k {
        let sigma_p = s_diag.read(p);
        for j in 0..n {
            v_vec.push(sigma_p * v_mat.read(j, p));
        }
    }
    *flops += (k * n) as u64;

    let mut singular_values = Vec::with_capacity(k);
    for p in 0..k {
        singular_values.push(s_diag.read(p));
    }

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

/// Reconstruct dense matrix from H-matrix tree
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
            let m = node.nrows();
            let n = node.ncols();
            let k = *rank;
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

/// Rapid matrix-vector multiplication y += H * x
pub fn hmatrix_matvec(node: &HNode, x: &[f64], y: &mut [f64], flops: &mut u64) {
    match &node.kind {
        HNodeKind::Zero => {
            // Nothing to add
        }
        HNodeKind::Leaf { rank, u, v, .. } => {
            let m = node.nrows();
            let n = node.ncols();
            let k = *rank;
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
