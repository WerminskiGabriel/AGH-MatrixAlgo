use crate::tree::{HNode, HNodeKind};
use faer::Mat;
use image::{ImageBuffer, Rgb, Rgba, RgbaImage};

/// Draw the hierarchical block structure of the H-matrix into an image file
/// Leaves are outlined, and rank-k leaves have a top horizontal band and left vertical band
/// forming a "Gamma" shape representing the factors U (width k) and V (height k).
pub fn draw_hmatrix_structure(
    node: &HNode,
    total_rows: usize,
    total_cols: usize,
    img_size: u32,
    out_path: &str,
) -> Result<(), image::ImageError> {
    let mut img = RgbaImage::from_pixel(img_size, img_size, Rgba([255, 255, 255, 255]));

    fn draw_recursive(
        n: &HNode,
        total_rows: usize,
        total_cols: usize,
        img_size: u32,
        img: &mut RgbaImage,
    ) {
        let x0 = ((n.s_min as f64 / total_cols as f64) * img_size as f64).round() as u32;
        let x1 = ((n.s_max as f64 / total_cols as f64) * img_size as f64).round() as u32;
        let y0 = ((n.t_min as f64 / total_rows as f64) * img_size as f64).round() as u32;
        let y1 = ((n.t_max as f64 / total_rows as f64) * img_size as f64).round() as u32;

        let bw = x1.saturating_sub(x0).max(1);
        let bh = y1.saturating_sub(y0).max(1);

        match &n.kind {
            HNodeKind::Zero => {
                for y in y0..y1.min(img_size) {
                    for x in x0..x1.min(img_size) {
                        img.put_pixel(x, y, Rgba([245, 245, 245, 255]));
                    }
                }
                draw_rect_border(img, x0, y0, bw, bh, Rgba([180, 180, 180, 255]));
            }
            HNodeKind::Leaf { rank, .. } => {
                let k = *rank;
                for y in y0..y1.min(img_size) {
                    for x in x0..x1.min(img_size) {
                        img.put_pixel(x, y, Rgba([255, 255, 255, 255]));
                    }
                }

                let thickness = ((k as u32 * 2).max(2)).min(bw.min(bh) / 2);
                let strip_color = Rgba([30, 30, 30, 255]);

                // Top horizontal strip (V factor)
                for y in y0..(y0 + thickness).min(y1) {
                    for x in x0..x1.min(img_size) {
                        img.put_pixel(x, y, strip_color);
                    }
                }

                // Left vertical strip (U factor)
                for y in y0..y1.min(img_size) {
                    for x in x0..(x0 + thickness).min(x1) {
                        img.put_pixel(x, y, strip_color);
                    }
                }

                draw_rect_border(img, x0, y0, bw, bh, Rgba([20, 20, 20, 255]));
            }
            HNodeKind::Internal { children } => {
                for child in children.iter() {
                    draw_recursive(child, total_rows, total_cols, img_size, img);
                }
                draw_rect_border(img, x0, y0, bw, bh, Rgba([60, 60, 60, 255]));
            }
        }
    }

    fn draw_rect_border(img: &mut RgbaImage, x0: u32, y0: u32, w: u32, h: u32, color: Rgba<u8>) {
        let x1 = (x0 + w).saturating_sub(1).min(img.width() - 1);
        let y1 = (y0 + h).saturating_sub(1).min(img.height() - 1);
        for x in x0..=x1 {
            img.put_pixel(x, y0, color);
            img.put_pixel(x, y1, color);
        }
        for y in y0..=y1 {
            img.put_pixel(x0, y, color);
            img.put_pixel(x1, y, color);
        }
    }

    draw_recursive(node, total_rows, total_cols, img_size, &mut img);
    img.save(out_path)
}

/// Save matrix as grayscale image
pub fn save_matrix_as_image(mat: &Mat<f64>, path: &str) -> Result<(), image::ImageError> {
    let rows = mat.nrows() as u32;
    let cols = mat.ncols() as u32;
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(cols, rows);
    for i in 0..rows {
        for j in 0..cols {
            let val = mat.read(i as usize, j as usize).clamp(0.0, 255.0).round() as u8;
            img.put_pixel(j, i, Rgb([val, val, val]));
        }
    }
    img.save(path)
}

/// Save 3 matrices (R, G, B) as a single RGB image
pub fn save_rgb_image(
    r_mat: &Mat<f64>,
    g_mat: &Mat<f64>,
    b_mat: &Mat<f64>,
    path: &str,
) -> Result<(), image::ImageError> {
    let rows = r_mat.nrows() as u32;
    let cols = r_mat.ncols() as u32;
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
