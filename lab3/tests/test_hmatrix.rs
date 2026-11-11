use faer::Mat;
use hmatrix_svd::svd::{calculate_psnr, frobenius_error};
use hmatrix_svd::tree::{create_tree, decompress_hmatrix, hmatrix_matvec};

#[test]
fn test_hmatrix_compress_decompress_4x4() {
    let m4_raw = vec![
        vec![4.0, 4.0, 1.0, 1.0],
        vec![4.0, 4.0, 1.0, 1.0],
        vec![2.0, 2.0, 8.0, 8.0],
        vec![2.0, 2.0, 8.0, 8.0],
    ];
    let m4 = Mat::from_fn(4, 4, |i, j| m4_raw[i][j]);

    let mut flops = 0u64;
    let tree = create_tree(&m4, 0, 4, 0, 4, 1, 0.5, &mut flops);
    assert_eq!(tree.nrows(), 4);
    assert_eq!(tree.ncols(), 4);

    let mut rec = Mat::zeros(4, 4);
    let mut d_flops = 0u64;
    decompress_hmatrix(&tree, &mut rec, &mut d_flops);

    let (diff, _, rel_diff) = frobenius_error(&m4, &rec);
    assert!(diff < 1e-10, "Frobenius difference too high: {}", diff);
    assert!(rel_diff < 1e-10, "Relative difference too high: {}", rel_diff);

    let psnr = calculate_psnr(&m4, &rec);
    assert!(psnr >= 99.0, "PSNR should be 99.0 (capped max) for exact reconstruction, got {}", psnr);
}

#[test]
fn test_hmatrix_matvec_equivalence() {
    let m = Mat::from_fn(8, 8, |i, j| (i + j + 1) as f64);
    let mut flops = 0u64;
    let tree = create_tree(&m, 0, 8, 0, 8, 2, 1e-4, &mut flops);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let mut y_dense = vec![0.0; 8];
    for i in 0..8 {
        for j in 0..8 {
            y_dense[i] += m.read(i, j) * x[j];
        }
    }

    let mut y_hmat = vec![0.0; 8];
    let mut hm_flops = 0u64;
    hmatrix_matvec(&tree, &x, &mut y_hmat, &mut hm_flops);

    for i in 0..8 {
        assert!(
            (y_dense[i] - y_hmat[i]).abs() < 1e-9,
            "Mismatch at index {}: dense={}, hmat={}",
            i,
            y_dense[i],
            y_hmat[i]
        );
    }
}
