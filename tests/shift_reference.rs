//! The fast kernel is checked against a deliberately naive reference:
//! O(N^2) DFT shifting, TV recomputed from scratch with modular indexing, no caching.

use approx::assert_abs_diff_eq;
use degibbs::{fft::Engine, shift::ShiftParams};
use std::f64::consts::PI;

fn shifted_naive(line: &[f64], sigma: f64, zero_nyquist: bool) -> Vec<f64> {
    let n = line.len();
    let freqs: Vec<i64> = (0..n as i64)
        .map(|j| {
            if j <= (n as i64 - 1) / 2 {
                j
            } else {
                j - n as i64
            }
        })
        .collect();
    // forward DFT
    let coef: Vec<(f64, f64)> = (0..n)
        .map(|j| {
            let (mut re, mut im) = (0.0, 0.0);
            for (x, &v) in line.iter().enumerate() {
                let ph = -2.0 * PI * ((j * x) % n) as f64 / n as f64;
                re += v * ph.cos();
                im += v * ph.sin();
            }
            (re, im)
        })
        .collect();
    (0..n)
        .map(|x| {
            let mut acc = 0.0;
            for j in 0..n {
                if zero_nyquist && n.is_multiple_of(2) && j == n / 2 {
                    continue;
                }
                let ph = 2.0
                    * PI
                    * (((j * x) % n) as f64 / n as f64 + freqs[j] as f64 * sigma / n as f64);
                acc += coef[j].0 * ph.cos() - coef[j].1 * ph.sin();
            }
            acc / n as f64
        })
        .collect()
}

fn reference(line: &[f64], p: &ShiftParams) -> Vec<f64> {
    let n = line.len() as i64;
    let m = p.nshifts as i64;
    let imgs: Vec<Vec<f64>> = (-m..=m)
        .map(|s| {
            if s == 0 {
                line.to_vec()
            } else {
                shifted_naive(line, s as f64 / (2 * m) as f64, true)
            }
        })
        .collect();
    let at = |img: &Vec<f64>, i: i64| img[i.rem_euclid(n) as usize];
    (0..n)
        .map(|x| {
            let mut best = (f64::INFINITY, 0usize);
            for (si, img) in imgs.iter().enumerate() {
                let dp: f64 = (p.min_w as i64..=p.max_w as i64)
                    .map(|k| (at(img, x + k) - at(img, x + k + 1)).abs())
                    .sum();
                let dm: f64 = (p.min_w as i64..=p.max_w as i64)
                    .map(|k| (at(img, x - k) - at(img, x - k - 1)).abs())
                    .sum();
                if dp < best.0 {
                    best = (dp, si);
                }
                if dm < best.0 {
                    best = (dm, si);
                }
            }
            let s = best.1 as i64 - m;
            let img = &imgs[best.1];
            if s == 0 {
                at(img, x)
            } else {
                let a = s.abs() as f64 / (2 * m) as f64;
                let nb = if s > 0 { x - 1 } else { x + 1 };
                (1.0 - a) * at(img, x) + a * at(img, nb)
            }
        })
        .collect()
}

fn pseudo_random_line(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    let mut edge = 0.0;
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let r = ((s >> 33) as f64) / (1u64 << 31) as f64;
            if r > 0.9 {
                edge = 1.0 - edge; // occasional steps so TV selection is non-trivial
            }
            edge + 0.2 * r
        })
        .collect()
}

#[test]
fn fast_kernel_matches_naive_reference_on_odd_and_even_lines() {
    let p = ShiftParams {
        nshifts: 6,
        min_w: 1,
        max_w: 3,
    };
    let mut e = Engine::new();
    for (n, seed) in [(16, 1), (17, 2), (24, 3), (31, 4), (40, 5)] {
        let line = pseudo_random_line(n, seed);
        let mut out = vec![0.0; n];
        e.unring_line(&line, &mut out, &p).unwrap();
        let want = reference(&line, &p);
        for (a, b) in out.iter().zip(&want) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-10);
        }
    }
}

#[test]
fn constant_line_is_unchanged() {
    let mut e = Engine::new();
    for n in [8, 9, 32] {
        let line = vec![3.5; n];
        let mut out = vec![0.0; n];
        e.unring_line(&line, &mut out, &ShiftParams::default())
            .unwrap();
        for v in out {
            assert_abs_diff_eq!(v, 3.5, epsilon = 1e-12);
        }
    }
}

#[test]
fn too_short_line_is_an_error() {
    let mut e = Engine::new();
    let line = vec![1.0; 5];
    let mut out = vec![0.0; 5];
    assert!(
        e.unring_line(&line, &mut out, &ShiftParams::default())
            .is_err()
    );
}

/// A box edge sampled from truncated k-space rings; after unringing the error must shrink.
#[test]
fn truncated_box_ringing_is_reduced() {
    use rustfft::num_complex::Complex64 as C;
    let n = 64usize;
    let (lo, hi) = (20.3, 43.6);
    // Fourier series of the indicator of [lo, hi), truncated symmetrically at |k| < n/2.
    let coef = |k: i64| -> C {
        if k == 0 {
            return C::new((hi - lo) / n as f64, 0.0);
        }
        let w = 2.0 * PI * k as f64 / n as f64;
        (C::from_polar(1.0, -w * lo) - C::from_polar(1.0, -w * hi)) / C::new(0.0, w * n as f64)
    };
    let kmax = n as i64 / 2 - 1;
    let line: Vec<f64> = (0..n)
        .map(|x| {
            (-kmax..=kmax)
                .map(|k| {
                    (coef(k) * C::from_polar(1.0, 2.0 * PI * k as f64 * x as f64 / n as f64)).re
                })
                .sum()
        })
        .collect();
    let truth: Vec<f64> = (0..n)
        .map(|x| {
            if (x as f64) >= lo && (x as f64) < hi {
                1.0
            } else {
                0.0
            }
        })
        .collect();
    // Compare away from the edges themselves (within 1.5 voxels the point-sampled truth is ambiguous).
    let far = |x: usize| (x as f64 - lo).abs() > 1.5 && (x as f64 - hi).abs() > 1.5;
    let rmse = |v: &[f64]| {
        let xs: Vec<usize> = (0..n).filter(|&x| far(x)).collect();
        (xs.iter().map(|&x| (v[x] - truth[x]).powi(2)).sum::<f64>() / xs.len() as f64).sqrt()
    };
    let mut e = Engine::new();
    let mut out = vec![0.0; n];
    e.unring_line(&line, &mut out, &ShiftParams::default())
        .unwrap();
    assert!(
        rmse(&out) < 0.5 * rmse(&line),
        "rmse before {} after {}",
        rmse(&line),
        rmse(&out)
    );
}
