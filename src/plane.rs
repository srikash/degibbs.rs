//! 2D plane operations (SPEC sections 2-3): axis-wise unringing and the full-Fourier "SuShi" pipeline.
//! Plane layout is row-major with y (the phase-encoding axis) as the slow axis: `data[y * nx + x]`.

use crate::{fft::Engine, shift::ShiftParams};
use anyhow::Result;
use rustfft::num_complex::Complex64 as C;
use std::f64::consts::PI;

#[derive(Clone, Debug, PartialEq)]
pub struct Plane {
    pub data: Vec<f64>,
    pub nx: usize,
    pub ny: usize,
}

impl Plane {
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            data: vec![0.0; nx * ny],
            nx,
            ny,
        }
    }
}

impl Plane {
    /// `self = a + b` (same dimensions).
    pub fn assign_sum(&mut self, a: &Plane, b: &Plane) {
        for ((o, x), y) in self.data.iter_mut().zip(&a.data).zip(&b.data) {
            *o = x + y;
        }
    }
}

/// Complex counterpart of [`Plane`], same layout: `data[y * nx + x]`.
#[derive(Clone, Debug, PartialEq)]
pub struct CPlane {
    pub data: Vec<C>,
    pub nx: usize,
    pub ny: usize,
}

impl CPlane {
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            data: vec![C::new(0.0, 0.0); nx * ny],
            nx,
            ny,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
}

/// `1 + cos(2 pi m / n)`, written as `2 cos^2(pi m / n)` so it is exactly representable near zero.
fn one_plus_cos(m: usize, n: usize) -> f64 {
    let c = (PI * m as f64 / n as f64).cos();
    2.0 * c * c
}

impl Engine {
    /// 1D-unring every line along `axis`, restricted for `Axis::Y` to rows `offset, offset+step, ...`
    /// (treated as an independent sub-image of those rows).
    pub fn unring_axis_sub(
        &mut self,
        plane: &mut Plane,
        axis: Axis,
        offset: usize,
        step: usize,
        p: &ShiftParams,
    ) -> Result<()> {
        let (nx, ny) = (plane.nx, plane.ny);
        match axis {
            Axis::X => {
                let mut out = vec![0.0; nx];
                for row in plane.data.chunks_exact_mut(nx) {
                    self.unring_line(row, &mut out, p)?;
                    row.copy_from_slice(&out);
                }
            }
            Axis::Y => {
                let rows: Vec<usize> = (offset..ny).step_by(step).collect();
                let len = rows.len();
                let mut col = vec![0.0; len];
                let mut out = vec![0.0; len];
                for x in 0..nx {
                    for (i, &y) in rows.iter().enumerate() {
                        col[i] = plane.data[y * nx + x];
                    }
                    self.unring_line(&col, &mut out, p)?;
                    for (i, &y) in rows.iter().enumerate() {
                        plane.data[y * nx + x] = out[i];
                    }
                }
            }
        }
        Ok(())
    }

    pub fn unring_axis(&mut self, plane: &mut Plane, axis: Axis, p: &ShiftParams) -> Result<()> {
        self.unring_axis_sub(plane, axis, 0, 1, p)
    }

    /// 2D spectrum of a plane (row-major, same layout).
    fn spectrum(&mut self, plane: &Plane) -> Vec<C> {
        let mut buf: Vec<C> = plane.data.iter().map(|&v| C::new(v, 0.0)).collect();
        self.fft2(&mut buf, plane.ny, plane.nx, false);
        buf
    }

    /// `IFT[spec * G]` with `Gx` (for `Axis::X`) or `Gy` (for `Axis::Y`) of SPEC section 3,
    /// where `spec` comes from [`Engine::spectrum`] for an `ny x nx` plane. Keeps the imaginary part.
    fn weighted_complex(&mut self, spec: &[C], ny: usize, nx: usize, which: Axis) -> Vec<C> {
        let cx: Vec<f64> = (0..nx).map(|m| one_plus_cos(m, nx)).collect();
        let cy: Vec<f64> = (0..ny).map(|m| one_plus_cos(m, ny)).collect();
        let mut buf = spec.to_vec();
        for my in 0..ny {
            for mx in 0..nx {
                let g = if 2 * mx == nx && 2 * my == ny {
                    0.5 // both cosines are -1: 0/0, split evenly
                } else {
                    let denom = cx[mx] + cy[my];
                    match which {
                        Axis::X => cy[my] / denom,
                        Axis::Y => cx[mx] / denom,
                    }
                };
                buf[my * nx + mx] *= g;
            }
        }
        self.fft2(&mut buf, ny, nx, true);
        buf
    }

    /// Real part of [`Engine::weighted_complex`], for real input.
    fn weighted(&mut self, spec: &[C], ny: usize, nx: usize, which: Axis) -> Plane {
        let buf = self.weighted_complex(spec, ny, nx, which);
        Plane {
            data: buf.iter().map(|c| c.re).collect(),
            nx,
            ny,
        }
    }

    /// `(I_x, I_y)` of SPEC section 3. Since `Gx + Gy = 1`, `I_y = I - I_x`, so only one inverse FFT is needed.
    pub fn weighted_pair(&mut self, plane: &Plane) -> (Plane, Plane) {
        let spec = self.spectrum(plane);
        let ix = self.weighted(&spec, plane.ny, plane.nx, Axis::X);
        let iy = Plane {
            data: plane
                .data
                .iter()
                .zip(&ix.data)
                .map(|(i, x)| i - x)
                .collect(),
            nx: plane.nx,
            ny: plane.ny,
        };
        (ix, iy)
    }

    /// Full-Fourier 2D unringing (Kellner's weighted x/y scheme): `J = unring_x(I_x) + unring_y(I_y)`.
    pub fn sushi(&mut self, plane: &mut Plane, p: &ShiftParams) -> Result<()> {
        let (mut jx, mut jy) = self.weighted_pair(plane);
        self.unring_axis(&mut jx, Axis::X, p)?;
        self.unring_axis(&mut jy, Axis::Y, p)?;
        plane.assign_sum(&jx, &jy);
        Ok(())
    }

    /// 1D-unring every line of a complex plane along `axis`.
    pub fn unring_axis_complex(
        &mut self,
        plane: &mut CPlane,
        axis: Axis,
        p: &ShiftParams,
    ) -> Result<()> {
        let (nx, ny) = (plane.nx, plane.ny);
        let len = match axis {
            Axis::X => nx,
            Axis::Y => ny,
        };
        let (mut re, mut im) = (vec![0.0; len], vec![0.0; len]);
        let (mut out_re, mut out_im) = (vec![0.0; len], vec![0.0; len]);
        let lines = match axis {
            Axis::X => ny,
            Axis::Y => nx,
        };
        for l in 0..lines {
            let at = |i: usize| match axis {
                Axis::X => l * nx + i,
                Axis::Y => i * nx + l,
            };
            for i in 0..len {
                let v = plane.data[at(i)];
                (re[i], im[i]) = (v.re, v.im);
            }
            self.unring_line_complex(&re, &im, &mut out_re, &mut out_im, p)?;
            for i in 0..len {
                plane.data[at(i)] = C::new(out_re[i], out_im[i]);
            }
        }
        Ok(())
    }

    /// Complex counterpart of [`Engine::sushi`]: the x/y weighting is linear, so it keeps both parts.
    pub fn sushi_complex(&mut self, plane: &mut CPlane, p: &ShiftParams) -> Result<()> {
        let (nx, ny) = (plane.nx, plane.ny);
        let mut spec = plane.data.clone();
        self.fft2(&mut spec, ny, nx, false);
        let ix = self.weighted_complex(&spec, ny, nx, Axis::X);
        let mut jx = CPlane { data: ix, nx, ny };
        let mut jy = CPlane {
            data: plane
                .data
                .iter()
                .zip(&jx.data)
                .map(|(i, x)| i - x)
                .collect(),
            nx,
            ny,
        };
        self.unring_axis_complex(&mut jx, Axis::X, p)?;
        self.unring_axis_complex(&mut jy, Axis::Y, p)?;
        for ((o, x), y) in plane.data.iter_mut().zip(&jx.data).zip(&jy.data) {
            *o = x + y;
        }
        Ok(())
    }
}
