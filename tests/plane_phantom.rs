mod common;

use approx::assert_abs_diff_eq;
use common::*;
use degibbs::{fft::Engine, pf::Pf, plane::Plane, shift::ShiftParams};

const R: usize = 8;

fn run(corrupted: &Plane, pf: Pf) -> Plane {
    let mut e = Engine::new();
    let mut p = corrupted.clone();
    e.rpg_plane(&mut p, pf, &ShiftParams::default()).unwrap();
    p
}

fn nu(pf: Pf) -> f64 {
    match pf {
        Pf::Full => 1.0,
        Pf::P78 => 0.75,
        Pf::P68 => 0.5,
        Pf::P58 => 0.25,
    }
}

#[test]
fn full_fourier_correction_reduces_error_on_all_plane_shapes() {
    // even x even (incl. the 0/0 weight corner), odd, and mixed sizes
    for (nx, ny) in [(64, 64), (96, 96), (65, 63), (50, 64), (63, 50)] {
        let (cor, truth) = corrupted_and_truth(nx, ny, R, 1.0);
        let out = run(&cor, Pf::Full);
        assert!(
            out.data.iter().all(|v| v.is_finite()),
            "non-finite output for {nx}x{ny}"
        );
        let (before, after) = (rmse_far(&cor, &truth, 2), rmse_far(&out, &truth, 2));
        assert!(after < 0.85 * before, "{nx}x{ny}: rmse {before} -> {after}");
    }
}

#[test]
fn pf_7_8_and_6_8_reduce_ringing_away_from_edges() {
    for pf in [Pf::P78, Pf::P68] {
        for (nx, ny) in [(64, 64), (96, 96), (65, 63), (50, 64)] {
            let (cor, truth) = corrupted_and_truth(nx, ny, R, nu(pf));
            let out = run(&cor, pf);
            assert!(out.data.iter().all(|v| v.is_finite()));
            let (before, after) = (rmse_far(&cor, &truth, 2), rmse_far(&out, &truth, 2));
            assert!(
                after < 0.9 * before,
                "{pf:?} {nx}x{ny}: far rmse {before} -> {after}"
            );
        }
    }
}

/// Lee 2021 reports that 5/8 keeps residual ringing and blurs; only require sane, bounded output.
#[test]
fn pf_5_8_is_finite_and_not_grossly_worse() {
    for (nx, ny) in [(64, 64), (96, 96), (65, 63)] {
        let (cor, truth) = corrupted_and_truth(nx, ny, R, nu(Pf::P58));
        let out = run(&cor, Pf::P58);
        assert!(out.data.iter().all(|v| v.is_finite()));
        assert!(rmse_far(&out, &truth, 2) < 1.4 * rmse_far(&cor, &truth, 2));
    }
}

#[test]
fn rpg_beats_full_fourier_correction_on_pf_data_for_7_8_and_6_8() {
    for pf in [Pf::P78, Pf::P68] {
        let (cor, truth) = corrupted_and_truth(96, 96, R, nu(pf));
        let rpg = rmse_far(&run(&cor, pf), &truth, 2);
        let full = rmse_far(&run(&cor, Pf::Full), &truth, 2);
        assert!(rpg < full, "{pf:?}: rpg {rpg} vs full-fourier {full}");
    }
}

#[test]
fn constant_plane_is_unchanged_for_every_mode() {
    for pf in [Pf::Full, Pf::P78, Pf::P68, Pf::P58] {
        let mut p = Plane::new(40, 36);
        p.data.fill(2.0);
        let out = run(&p, pf);
        for v in out.data {
            assert_abs_diff_eq!(v, 2.0, epsilon = 1e-9);
        }
    }
}

#[test]
fn transposing_the_plane_transposes_the_full_fourier_result() {
    let (cor, _) = corrupted_and_truth(60, 44, R, 1.0);
    let a = transpose(&run(&cor, Pf::Full));
    let b = run(&transpose(&cor), Pf::Full);
    for (x, y) in a.data.iter().zip(&b.data) {
        assert_abs_diff_eq!(x, y, epsilon = 1e-9);
    }
}

#[test]
fn too_small_plane_is_rejected() {
    let mut e = Engine::new();
    let mut p = Plane::new(40, 6);
    assert!(
        e.rpg_plane(&mut p, Pf::P78, &ShiftParams::default())
            .is_err()
    );
}

/// Guards against a pipeline that "reduces ringing" simply by blurring: edge-pixel error must stay
/// close to the uncorrected data and to full-Fourier correction, and nothing may shift along y.
#[test]
fn pf_7_8_and_6_8_do_not_degrade_edges_or_shift_the_image() {
    for pf in [Pf::P78, Pf::P68] {
        let (cor, truth) = corrupted_and_truth(96, 96, R, nu(pf));
        let out = run(&cor, pf);
        let full = run(&cor, Pf::Full);
        let (raw, rpg, fc) = (
            rmse_near(&cor, &truth, 2),
            rmse_near(&out, &truth, 2),
            rmse_near(&full, &truth, 2),
        );
        assert!(rpg < 1.3 * raw, "{pf:?}: edge rmse raw {raw} -> rpg {rpg}");
        assert!(
            rpg < 1.3 * fc,
            "{pf:?}: edge rmse full-fourier {fc} vs rpg {rpg}"
        );
        assert!((centroid_y(&out) - centroid_y(&cor)).abs() < 0.05);
    }
}

/// Real data has phase: with a smooth phase map, zeroing either side of k-space must still be improved.
#[test]
fn pf_with_smooth_phase_is_improved_for_either_unsampled_side() {
    for pf in [Pf::P78, Pf::P68] {
        for side in [-1, 1] {
            let (cor, truth) = corrupted_and_truth_with(96, 96, R, nu(pf), side, true);
            let out = run(&cor, pf);
            assert!(out.data.iter().all(|v| v.is_finite()));
            let (before, after) = (rmse_far(&cor, &truth, 2), rmse_far(&out, &truth, 2));
            assert!(
                after < 0.9 * before,
                "{pf:?} side {side}: far rmse {before} -> {after}"
            );
        }
    }
}
