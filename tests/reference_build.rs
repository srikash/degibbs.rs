//! Compares the current build with a reference build (e.g. 0.2.1) on a corrupted phantom volume. Run by hand:
//!   DEGIBBS_REF_BIN=/path/to/reference/degibbs cargo test --release --test reference_build -- --ignored --nocapture
mod common;

use common::{cli, corrupted_and_truth, save};
use ndarray::{ArrayD, IxDyn, ShapeBuilder};

/// F-ordered `[nx, ny, nz, nt]` volume whose planes are scaled copies of one corrupted phantom.
fn volume(nx: usize, ny: usize, nz: usize, nt: usize) -> ArrayD<f32> {
    let (plane, _) = corrupted_and_truth(nx, ny, 8, 0.75);
    let mut v = ArrayD::<f32>::zeros(IxDyn(&[nx, ny, nz, nt]).f());
    for z in 0..nz {
        for t in 0..nt {
            for y in 0..ny {
                for x in 0..nx {
                    v[[x, y, z, t]] = (plane.data[y * nx + x]
                        * (1.0 + 0.1 * z as f64 + 0.3 * t as f64)
                        + 0.01 * ((x * 7 + y * 13 + z * 3 + t) % 17) as f64)
                        as f32;
                }
            }
        }
    }
    v
}

#[test]
#[ignore = "needs DEGIBBS_REF_BIN"]
fn output_matches_the_reference_build_within_tolerance() {
    let reference = std::env::var("DEGIBBS_REF_BIN").expect("set DEGIBBS_REF_BIN");
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.nii");
    for (nx, ny) in [(64usize, 64usize), (63, 61)] {
        save(&input, volume(nx, ny, 4, 3), None);
        for pf in ["1", "7/8", "6/8", "5/8"] {
            let (old, new) = (dir.path().join("old.nii"), dir.path().join("new.nii"));
            let _ = (
                std::fs::remove_file(&old),
                std::fs::remove_file(&new),
                std::fs::remove_file(dir.path().join("new.json")),
            );
            let args = ["--pe-axis", "j", "--pf", pf];
            let st = std::process::Command::new(&reference)
                .arg(&input)
                .arg(&old)
                .args(args)
                .status()
                .unwrap();
            assert!(st.success());
            cli().arg(&input).arg(&new).args(args).assert().success();
            let (a, b) = (
                degibbs::io::read(&old).unwrap().data,
                degibbs::io::read(&new).unwrap().data,
            );
            let range = a.iter().cloned().fold(f32::MIN, f32::max)
                - a.iter().cloned().fold(f32::MAX, f32::min);
            let mut worst = 0f32;
            let mut differing = 0usize;
            for (x, y) in a.iter().zip(b.iter()) {
                let d = (x - y).abs();
                if d > 0.0 {
                    differing += 1;
                }
                worst = worst.max(d);
            }
            println!(
                "{nx}x{ny} pf {pf}: {differing} of {} voxels differ, largest difference {worst:e} ({:e} of the range {range})",
                a.len(),
                worst / range
            );
            assert!(worst <= 1e-5 * range, "pf {pf}: largest difference {worst}");
        }
    }
}
