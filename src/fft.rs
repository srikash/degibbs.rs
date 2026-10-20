//! FFT engine: cached rustfft plans plus reusable scratch buffers.
//! One `Engine` per worker thread; the algorithm modules add `impl Engine` blocks.

use rustfft::{Fft, FftPlanner, num_complex::Complex64};
use std::{collections::HashMap, sync::Arc};

type C = Complex64;

pub struct Engine {
    planner: FftPlanner<f64>,
    fwd: HashMap<usize, Arc<dyn Fft<f64>>>,
    inv: HashMap<usize, Arc<dyn Fft<f64>>>,
    scratch: Vec<C>,
    column: Vec<C>,
    /// Scratch owned by the 1D kernel (see `shift.rs`).
    pub(crate) shifted: Vec<f64>,
    pub(crate) absdiff: Vec<f64>,
    pub(crate) spec: Vec<C>,
    pub(crate) work: Vec<C>,
}

impl Default for Engine {
    fn default() -> Self {
        Self::new()
    }
}

impl Engine {
    pub fn new() -> Self {
        Self {
            planner: FftPlanner::new(),
            fwd: HashMap::new(),
            inv: HashMap::new(),
            scratch: Vec::new(),
            column: Vec::new(),
            shifted: Vec::new(),
            absdiff: Vec::new(),
            spec: Vec::new(),
            work: Vec::new(),
        }
    }

    fn plan(&mut self, n: usize, inverse: bool) -> Arc<dyn Fft<f64>> {
        let (map, planner) = if inverse {
            (&mut self.inv, &mut self.planner)
        } else {
            (&mut self.fwd, &mut self.planner)
        };
        map.entry(n)
            .or_insert_with(|| {
                if inverse {
                    planner.plan_fft_inverse(n)
                } else {
                    planner.plan_fft_forward(n)
                }
            })
            .clone()
    }

    /// In-place 1D transform of `buf`; the inverse is unnormalised.
    pub fn fft1(&mut self, buf: &mut [C], inverse: bool) {
        let plan = self.plan(buf.len(), inverse);
        let need = plan.get_inplace_scratch_len();
        if self.scratch.len() < need {
            self.scratch.resize(need, C::new(0.0, 0.0));
        }
        plan.process_with_scratch(buf, &mut self.scratch[..need]);
    }

    /// In-place 2D transform of a row-major `ny x nx` array (`data[y*nx + x]`).
    /// The inverse is normalised by `1/(nx*ny)` so that `fft2 -> ifft2` is the identity.
    pub fn fft2(&mut self, data: &mut [C], ny: usize, nx: usize, inverse: bool) {
        assert_eq!(data.len(), ny * nx);
        for row in data.chunks_exact_mut(nx) {
            self.fft1(row, inverse);
        }
        let mut column = std::mem::take(&mut self.column);
        column.resize(ny, C::new(0.0, 0.0));
        for x in 0..nx {
            for y in 0..ny {
                column[y] = data[y * nx + x];
            }
            self.fft1(&mut column, inverse);
            for y in 0..ny {
                data[y * nx + x] = column[y];
            }
        }
        self.column = column;
        if inverse {
            let s = 1.0 / (nx * ny) as f64;
            data.iter_mut().for_each(|v| *v *= s);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn naive_dft(x: &[C], inverse: bool) -> Vec<C> {
        let n = x.len();
        let sign = if inverse { 1.0 } else { -1.0 };
        (0..n)
            .map(|k| {
                (0..n)
                    .map(|j| {
                        let ph =
                            sign * 2.0 * std::f64::consts::PI * ((k * j) % n) as f64 / n as f64;
                        x[j] * C::from_polar(1.0, ph)
                    })
                    .sum()
            })
            .collect()
    }

    fn signal(n: usize) -> Vec<C> {
        (0..n)
            .map(|i| {
                C::new(
                    ((i * 7 + 3) % 11) as f64 - 5.0,
                    ((i * 5 + 1) % 13) as f64 * 0.3,
                )
            })
            .collect()
    }

    #[test]
    fn matches_naive_dft_for_odd_even_and_prime_lengths() {
        let mut e = Engine::new();
        for n in [1, 2, 3, 4, 5, 7, 8, 9, 12, 13, 24, 31, 32, 72, 97] {
            let x = signal(n);
            for inverse in [false, true] {
                let mut y = x.clone();
                e.fft1(&mut y, inverse);
                let want = naive_dft(&x, inverse);
                for (a, b) in y.iter().zip(&want) {
                    assert_abs_diff_eq!(a.re, b.re, epsilon = 1e-9);
                    assert_abs_diff_eq!(a.im, b.im, epsilon = 1e-9);
                }
            }
        }
    }

    #[test]
    fn fft2_roundtrip_is_identity() {
        let mut e = Engine::new();
        for (ny, nx) in [(9, 8), (8, 8), (7, 12), (13, 5)] {
            let x = signal(ny * nx);
            let mut y = x.clone();
            e.fft2(&mut y, ny, nx, false);
            e.fft2(&mut y, ny, nx, true);
            for (a, b) in y.iter().zip(&x) {
                assert_abs_diff_eq!(a.re, b.re, epsilon = 1e-9);
                assert_abs_diff_eq!(a.im, b.im, epsilon = 1e-9);
            }
        }
    }
}
