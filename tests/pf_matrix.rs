//! Acquired-PF x correction-PF matrix on the analytic phantom: which `--pf` setting works on which data.
//! Run with `--nocapture` to print the far-field RMSE table.
mod common;

use common::*;
use degibbs::{fft::Engine, pf::Pf, plane::Plane, shift::ShiftParams};

const R: usize = 8;
const PFS: [Pf; 4] = [Pf::Full, Pf::P78, Pf::P68, Pf::P58];
const NAMES: [&str; 4] = ["1", "7/8", "6/8", "5/8"];

fn nu(pf: Pf) -> f64 {
    match pf {
        Pf::Full => 1.0,
        Pf::P78 => 0.75,
        Pf::P68 => 0.5,
        Pf::P58 => 0.25,
    }
}

fn run(cor: &Plane, pf: Pf) -> Plane {
    let mut p = cor.clone();
    Engine::new()
        .rpg_plane(&mut p, pf, &ShiftParams::default())
        .unwrap();
    p
}

/// `[acquired][0] = uncorrected`, `[acquired][1 + correction]` = far-field RMSE after correction.
fn matrix(nx: usize, ny: usize, side: i64) -> [[f64; 5]; 4] {
    let mut m = [[0.0; 5]; 4];
    for (i, &acq) in PFS.iter().enumerate() {
        let (cor, truth) = corrupted_and_truth_with(nx, ny, R, nu(acq), side, true);
        m[i][0] = rmse_far(&cor, &truth, 2);
        for (j, &c) in PFS.iter().enumerate() {
            m[i][1 + j] = rmse_far(&run(&cor, c), &truth, 2);
        }
    }
    eprintln!(
        "{nx}x{ny} side {side}: far-field RMSE (rows: acquired PF, cols: none, then correction PF)"
    );
    for i in 0..4 {
        eprintln!(
            "  {:>3} {:6.3} | {:6.3} {:6.3} {:6.3} {:6.3}",
            NAMES[i], m[i][0], m[i][1], m[i][2], m[i][3], m[i][4]
        );
    }
    m
}

#[test]
fn matching_correction_helps_and_is_best_for_7_8_and_6_8() {
    for (nx, ny) in [(96, 96), (64, 64)] {
        let m = matrix(nx, ny, -1);
        // full Fourier: clear gain
        assert!(m[0][1] < 0.8 * m[0][0], "{nx}x{ny} full: {:?}", m[0]);
        for i in [1, 2] {
            let (none, matching) = (m[i][0], m[i][1 + i]);
            assert!(matching < 0.95 * none, "{nx}x{ny} {}: {:?}", NAMES[i], m[i]);
            // the matching setting is at least as good as treating the data as full Fourier
            assert!(
                matching <= m[i][1] * 1.005,
                "{nx}x{ny} {}: {:?}",
                NAMES[i],
                m[i]
            );
        }
    }
}

#[test]
fn five_eighths_setting_does_not_help_other_data() {
    let m = matrix(96, 96, -1);
    for i in 0..3 {
        // 5/8 setting on 1, 7/8, 6/8 data is worse than doing nothing
        assert!(m[i][4] > m[i][0], "{}: {:?}", NAMES[i], m[i]);
    }
    // 5/8 data: no setting gives a large gain
    assert!(m[3][1..].iter().all(|&v| v > 0.9 * m[3][0]), "{:?}", m[3]);
}

#[test]
fn unsampled_side_does_not_change_the_picture() {
    let (neg, pos) = (matrix(96, 96, -1), matrix(96, 96, 1));
    for i in 0..4 {
        for j in 0..5 {
            let (a, b) = (neg[i][j], pos[i][j]);
            assert!((a - b).abs() < 0.15 * a.max(b), "[{i}][{j}] {a} vs {b}");
        }
    }
}
