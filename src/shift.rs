//! 1D local subvoxel-shift kernel (SPEC section 2).

use crate::fft::Engine;
use anyhow::{Result, bail};
use rustfft::num_complex::Complex64 as C;
use std::f64::consts::PI;

/// Kernel parameters: `nshifts` = M (shifts `s = -M..=M`, step `1/(2M)` voxel),
/// oscillation window `[min_w, max_w]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShiftParams {
    pub nshifts: usize,
    pub min_w: usize,
    pub max_w: usize,
}

impl Default for ShiftParams {
    fn default() -> Self {
        Self {
            nshifts: 20,
            min_w: 1,
            max_w: 3,
        }
    }
}

impl ShiftParams {
    pub fn validate(&self) -> Result<()> {
        if self.nshifts == 0 {
            bail!("nshifts must be at least 1");
        }
        if self.nshifts > 1000 || self.max_w > 100 {
            bail!("nshifts must be at most 1000 and maxw at most 100");
        }
        if self.min_w > self.max_w {
            bail!(
                "minw ({}) must not exceed maxw ({})",
                self.min_w,
                self.max_w
            );
        }
        Ok(())
    }

    /// Shortest line the kernel accepts (SPEC: `N > k2 + 2`).
    pub fn min_len(&self) -> usize {
        self.max_w + 3
    }
}

impl Engine {
    /// Remove Gibbs ringing from one real line: `out[x]` is the corrected value of `line[x]`.
    pub fn unring_line(&mut self, line: &[f64], out: &mut [f64], p: &ShiftParams) -> Result<()> {
        let n = line.len();
        assert_eq!(out.len(), n);
        p.validate()?;
        if n < p.min_len() {
            bail!(
                "line of length {n} is too short for maxw={} (need at least {})",
                p.max_w,
                p.min_len()
            );
        }
        let m = p.nshifts;
        let ns = 2 * m + 1;
        let mut shifted = std::mem::take(&mut self.shifted);
        let mut absdiff = std::mem::take(&mut self.absdiff);
        let mut spec = std::mem::take(&mut self.spec);
        let mut work = std::mem::take(&mut self.work);
        shifted.resize(ns * n, 0.0);
        absdiff.resize(ns * n, 0.0);
        spec.resize(n, C::new(0.0, 0.0));
        work.resize(n, C::new(0.0, 0.0));

        for (c, &v) in spec.iter_mut().zip(line) {
            *c = C::new(v, 0.0);
        }
        self.fft1(&mut spec, false);

        let inv_n = 1.0 / n as f64;
        for si in 0..ns {
            let s = si as isize - m as isize;
            let row = &mut shifted[si * n..(si + 1) * n];
            if s == 0 {
                row.copy_from_slice(line);
                continue;
            }
            let sigma = s as f64 / (2 * m) as f64;
            for (j, w) in work.iter_mut().enumerate() {
                if n.is_multiple_of(2) && j == n / 2 {
                    *w = C::new(0.0, 0.0); // Nyquist zeroed in shifted copies
                    continue;
                }
                let freq = if j <= (n - 1) / 2 {
                    j as f64
                } else {
                    j as f64 - n as f64
                };
                *w = spec[j] * C::from_polar(1.0, 2.0 * PI * freq * sigma / n as f64);
            }
            self.fft1(&mut work, true);
            for (r, w) in row.iter_mut().zip(work.iter()) {
                *r = w.re * inv_n;
            }
        }
        for si in 0..ns {
            let row = &shifted[si * n..(si + 1) * n];
            let d = &mut absdiff[si * n..(si + 1) * n];
            for j in 0..n {
                d[j] = (row[j] - row[(j + 1) % n]).abs();
            }
        }

        for x in 0..n {
            let mut best = f64::INFINITY;
            let mut best_si = m; // s = 0
            for si in 0..ns {
                let d = &absdiff[si * n..(si + 1) * n];
                let (mut dp, mut dm) = (0.0, 0.0);
                for k in p.min_w..=p.max_w {
                    dp += d[(x + k) % n];
                    dm += d[(x + n - k - 1) % n];
                }
                if dp < best {
                    best = dp;
                    best_si = si;
                }
                if dm < best {
                    best = dm;
                    best_si = si;
                }
            }
            let s = best_si as isize - m as isize;
            let row = &shifted[best_si * n..(best_si + 1) * n];
            out[x] = if s == 0 {
                row[x]
            } else {
                let a = s.unsigned_abs() as f64 / (2 * m) as f64;
                let nb = if s > 0 { (x + n - 1) % n } else { (x + 1) % n };
                (1.0 - a) * row[x] + a * row[nb]
            };
        }

        self.shifted = shifted;
        self.absdiff = absdiff;
        self.spec = spec;
        self.work = work;
        Ok(())
    }
}
