use faer::Mat;
use image::{GenericImageView, ImageBuffer, Rgb};

/// Load RGB image from file and split into three separate f64 matrices
pub fn load_rgb_image(
    path: &str,
) -> Result<(Mat<f64>, Mat<f64>, Mat<f64>), Box<dyn std::error::Error>> {
    let img = image::open(path)?;
    let (width, height) = img.dimensions();
    let mut r_mat = Mat::<f64>::zeros(height as usize, width as usize);
    let mut g_mat = Mat::<f64>::zeros(height as usize, width as usize);
    let mut b_mat = Mat::<f64>::zeros(height as usize, width as usize);

    for (x, y, pixel) in img.pixels() {
        let (i, j) = (y as usize, x as usize);
        r_mat.write(i, j, pixel[0] as f64);
        g_mat.write(i, j, pixel[1] as f64);
        b_mat.write(i, j, pixel[2] as f64);
    }
    Ok((r_mat, g_mat, b_mat))
}

pub fn save_matrix_as_image(mat: &Mat<f64>, path: &str) -> Result<(), image::ImageError> {
    let (rows, cols) = (mat.nrows() as u32, mat.ncols() as u32);
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(cols, rows);
    for i in 0..rows {
        for j in 0..cols {
            let val = mat.read(i as usize, j as usize).clamp(0.0, 255.0).round() as u8;
            img.put_pixel(j, i, Rgb([val, val, val]));
        }
    }
    img.save(path)
}

pub fn save_rgb_image(
    r_mat: &Mat<f64>,
    g_mat: &Mat<f64>,
    b_mat: &Mat<f64>,
    path: &str,
) -> Result<(), image::ImageError> {
    let (rows, cols) = (r_mat.nrows() as u32, r_mat.ncols() as u32);
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(cols, rows);
    for i in 0..rows {
        for j in 0..cols {
            let r = r_mat.read(i as usize, j as usize).clamp(0.0, 255.0).round() as u8;
            let g = g_mat.read(i as usize, j as usize).clamp(0.0, 255.0).round() as u8;
            let b = b_mat.read(i as usize, j as usize).clamp(0.0, 255.0).round() as u8;
            img.put_pixel(j, i, Rgb([r, g, b]));
        }
    }
    img.save(path)
}
