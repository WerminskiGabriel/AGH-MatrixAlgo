use crate::tree::{HNode, HNodeKind};
use image::{Rgba, RgbaImage};

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
                let thickness = ((k as u32 * 2).max(2)).min(bw.min(bh) / 2);
                let strip_color = Rgba([30, 30, 30, 255]);

                for y in y0..(y0 + thickness).min(y1) {
                    for x in x0..x1.min(img_size) {
                        img.put_pixel(x, y, strip_color);
                    }
                }
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

    draw_recursive(node, total_rows, total_cols, img_size, &mut img);
    img.save(out_path)
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
