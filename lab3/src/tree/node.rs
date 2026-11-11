use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum HNodeKind {
    /// Zero block (all values below threshold or zero)
    Zero,
    /// Low-rank compressed leaf storing U (m x rank) and V = diag(S) * V^T (rank x n)
    Leaf {
        rank: usize,
        singular_values: Vec<f64>,
        u: Vec<f64>, // Row-major: m * rank
        v: Vec<f64>, // Row-major: rank * n
    },
    /// Internal quadtree node with 4 children: [TL, TR, BL, BR]
    Internal { children: Box<[HNode; 4]> },
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
    #[inline]
    pub fn nrows(&self) -> usize {
        self.t_max - self.t_min
    }

    #[inline]
    pub fn ncols(&self) -> usize {
        self.s_max - self.s_min
    }

    pub fn node_counts(&self) -> (usize, usize, usize) {
        match &self.kind {
            HNodeKind::Zero => (1, 0, 0),
            HNodeKind::Leaf { .. } => (0, 1, 0),
            HNodeKind::Internal { children } => {
                let (mut z, mut r, mut i) = (0, 0, 1);
                for c in children.iter() {
                    let (cz, cr, ci) = c.node_counts();
                    z += cz;
                    r += cr;
                    i += ci;
                }
                (z, r, i)
            }
        }
    }

    pub fn depth(&self) -> usize {
        match &self.kind {
            HNodeKind::Zero | HNodeKind::Leaf { .. } => 1,
            HNodeKind::Internal { children } => {
                1 + children.iter().map(|c| c.depth()).max().unwrap_or(0)
            }
        }
    }

    pub fn stored_floats(&self) -> usize {
        match &self.kind {
            HNodeKind::Zero => 0,
            HNodeKind::Leaf { rank, .. } => *rank * (self.nrows() + self.ncols()) + *rank,
            HNodeKind::Internal { children } => children.iter().map(|c| c.stored_floats()).sum(),
        }
    }

    pub fn memory_bytes(&self) -> usize {
        let base = std::mem::size_of::<HNode>();
        match &self.kind {
            HNodeKind::Zero => base,
            HNodeKind::Leaf {
                u,
                v,
                singular_values,
                ..
            } => {
                base + (u.len() + v.len() + singular_values.len()) * std::mem::size_of::<f64>()
            }
            HNodeKind::Internal { children } => {
                base + std::mem::size_of::<[HNode; 4]>()
                    + children.iter().map(|c| c.memory_bytes()).sum::<usize>()
            }
        }
    }
}
