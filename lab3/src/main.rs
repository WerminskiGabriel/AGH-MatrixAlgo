use hmatrix_svd::experiments::runner::run_all_experiments;
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    let img_path = if args.len() > 1 {
        args[1].as_str()
    } else {
        "assets/input_512.png"
    };

    println!("Hierarchical Matrix Compression CLI");
    println!("Input image path: {}", img_path);

    run_all_experiments(img_path)?;

    println!("All experiments and visualizations finished successfully.");
    Ok(())
}
