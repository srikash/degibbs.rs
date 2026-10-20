//! Test phantoms: a high-resolution analytic-geometry image is truncated in k-space to make
//! exactly Gibbs-corrupted low-resolution data, and block-averaged to make the ground truth.
#![allow(dead_code)]

use degibbs::{fft::Engine, plane::Plane};
use rustfft::num_complex::Complex64 as C;

struct Ellipse {
    cx: f64,
    cy: f64,
    rx: f64,
    ry: f64,
    angle: f64,
    value: f64,
}

/// Phantom in unit coordinates (0..1 across the FOV); later shapes overwrite earlier ones.
/// `shift` moves everything by a fraction of a low-res voxel so edges sit at varied sub-voxel positions.
fn shapes() -> Vec<Ellipse> {
    vec![
        Ellipse {
            cx: 0.5,
            cy: 0.5,
            rx: 0.36,
            ry: 0.42,
            angle: 0.3,
            value: 1.0,
        },
        Ellipse {
            cx: 0.42,
            cy: 0.58,
            rx: 0.12,
            ry: 0.2,
            angle: -0.5,
            value: 0.35,
        },
        Ellipse {
            cx: 0.58,
            cy: 0.38,
            rx: 0.1,
            ry: 0.07,
            angle: 0.8,
            value: 0.7,
        },
    ]
}

pub fn phantom_hi(nx: usize, ny: usize, r: usize) -> Vec<f64> {
    let (hx, hy) = (nx * r, ny * r);
    let shapes = shapes();
    let mut img = vec![0.0; hx * hy];
    for y in 0..hy {
        for x in 0..hx {
            let (u, v) = ((x as f64 + 0.5) / hx as f64, (y as f64 + 0.5) / hy as f64);
            for s in &shapes {
                let (dx, dy) = (u - s.cx, v - s.cy);
                let (c, sn) = (s.angle.cos(), s.angle.sin());
                let (a, b) = (dx * c + dy * sn, -dx * sn + dy * c);
                if (a / s.rx).powi(2) + (b / s.ry).powi(2) <= 1.0 {
                    img[y * hx + x] = s.value;
                }
            }
        }
    }
    img
}

/// Signed frequency of FFT index `j` for length `n` (same convention as the library).
pub fn freq(j: usize, n: usize) -> i64 {
    if j <= (n - 1) / 2 {
        j as i64
    } else {
        j as i64 - n as i64
    }
}

/// Returns `(corrupted magnitude, truth)` on an `nx x ny` grid. `nu = 2*pf - 1` (1 for full Fourier):
/// ky samples with `-ky > nu*ny/2` (negative side) are zeroed.
pub fn corrupted_and_truth(nx: usize, ny: usize, r: usize, nu: f64) -> (Plane, Plane) {
    corrupted_and_truth_with(nx, ny, r, nu, -1, false)
}

/// General form: `side = -1` zeroes the negative-ky side, `+1` the positive side; `phase` multiplies the
/// object by a smooth phase map (the ground truth stays the magnitude).
pub fn corrupted_and_truth_with(
    nx: usize,
    ny: usize,
    r: usize,
    nu: f64,
    side: i64,
    phase: bool,
) -> (Plane, Plane) {
    let (hx, hy) = (nx * r, ny * r);
    let hi = phantom_hi(nx, ny, r);
    let mut e = Engine::new();
    let mut h: Vec<C> = hi
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            if !phase {
                return C::new(v, 0.0);
            }
            let (u, w) = (
                ((i % hx) as f64 + 0.5) / hx as f64,
                ((i / hx) as f64 + 0.5) / hy as f64,
            );
            let phi = 2.0
                * std::f64::consts::PI
                * (0.6 * u + 0.4 * (-((u - 0.5).powi(2) + (w - 0.5).powi(2)) / 0.1).exp());
            C::from_polar(v, phi)
        })
        .collect();
    e.fft2(&mut h, hy, hx, false);
    let mut low = vec![C::new(0.0, 0.0); nx * ny];
    for jy in 0..ny {
        let fy = freq(jy, ny);
        if fy * side > 0 && (fy.abs() as f64) > nu * ny as f64 / 2.0 {
            continue;
        }
        for jx in 0..nx {
            let fx = freq(jx, nx);
            let (iy, ix) = (
                fy.rem_euclid(hy as i64) as usize,
                fx.rem_euclid(hx as i64) as usize,
            );
            low[jy * nx + jx] = h[iy * hx + ix] / (r * r) as f64;
        }
    }
    // The hi-res forward FFT is unnormalised and `fft2(inverse)` divides by nx*ny, matching the scaling above.
    e.fft2(&mut low, ny, nx, true);
    let corrupted = Plane {
        data: low.iter().map(|c| c.norm()).collect(),
        nx,
        ny,
    };

    let mut truth = Plane::new(nx, ny);
    for y in 0..ny {
        for x in 0..nx {
            let mut acc = 0.0;
            for dy in 0..r {
                for dx in 0..r {
                    let yy = (y * r + dy + hy - r / 2) % hy;
                    let xx = (x * r + dx + hx - r / 2) % hx;
                    acc += hi[yy * hx + xx];
                }
            }
            truth.data[y * nx + x] = acc / (r * r) as f64;
        }
    }
    (corrupted, truth)
}

pub fn rmse(a: &Plane, truth: &Plane) -> f64 {
    (a.data
        .iter()
        .zip(&truth.data)
        .map(|(x, t)| (x - t).powi(2))
        .sum::<f64>()
        / a.data.len() as f64)
        .sqrt()
}

pub fn transpose(p: &Plane) -> Plane {
    let mut t = Plane::new(p.ny, p.nx);
    for y in 0..p.ny {
        for x in 0..p.nx {
            t.data[x * p.ny + y] = p.data[y * p.nx + x];
        }
    }
    t
}

/// True where the `(2w+1)^2` neighbourhood of `(y, x)` in `truth` is constant (no edge nearby).
fn is_flat(truth: &Plane, y: usize, x: usize, w: usize) -> bool {
    let (nx, ny) = (truth.nx, truth.ny);
    let t0 = truth.data[y * nx + x];
    (0..=2 * w).all(|dy| {
        (0..=2 * w).all(|dx| {
            let (yy, xx) = ((y + ny + dy - w) % ny, (x + nx + dx - w) % nx);
            (truth.data[yy * nx + xx] - t0).abs() < 1e-9
        })
    })
}

fn rmse_where(a: &Plane, truth: &Plane, w: usize, want_flat: bool) -> f64 {
    let (mut acc, mut cnt) = (0.0, 0usize);
    for y in 0..a.ny {
        for x in 0..a.nx {
            if is_flat(truth, y, x, w) == want_flat {
                acc += (a.data[y * a.nx + x] - truth.data[y * a.nx + x]).powi(2);
                cnt += 1;
            }
        }
    }
    (acc / cnt as f64).sqrt()
}

/// RMSE over pixels far from any edge: pure ringing, edge blur excluded.
pub fn rmse_far(a: &Plane, truth: &Plane, w: usize) -> f64 {
    rmse_where(a, truth, w, true)
}

/// RMSE over edge pixels (blur and bias).
pub fn rmse_near(a: &Plane, truth: &Plane, w: usize) -> f64 {
    rmse_where(a, truth, w, false)
}

/// Intensity-weighted mean y position, to detect systematic shifts along the PE axis.
pub fn centroid_y(a: &Plane) -> f64 {
    let (mut s, mut w) = (0.0, 0.0);
    for y in 0..a.ny {
        for x in 0..a.nx {
            let v = a.data[y * a.nx + x];
            s += y as f64 * v;
            w += v;
        }
    }
    s / w
}
