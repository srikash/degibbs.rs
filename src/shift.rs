//! 1D local subvoxel-shift kernel (SPEC section 2).

use crate::fft::Engine;
use anyhow::{Result, bail};
use rustfft::num_complex::Complex64 as C;
use std::{collections::HashMap, f64::consts::PI};

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

/// Working memory of the 1D kernel, kept between calls so lines of the same length reuse it.
#[derive(Default)]
pub(crate) struct KernelScratch {
    shifted: Vec<f64>,
    /// Imaginary parts of the shifted copies (complex input only).
    shifted_im: Vec<f64>,
    absdiff: Vec<f64>,
    spec: Vec<C>,
    work: Vec<C>,
    /// One absdiff row with `max_w + 1` wrapped samples on each side, so the window sums need no modulo.
    ext: Vec<f64>,
    best: Vec<f64>,
    best_si: Vec<usize>,
    /// `exp(2 pi i * freq * sigma_s / n)` for every shift and FFT index, per `(n, nshifts)`.
    phases: HashMap<(usize, usize), Vec<C>>,
}

fn phase_table(n: usize, m: usize) -> Vec<C> {
    let ns = 2 * m + 1;
    let mut t = Vec::with_capacity(ns * n);
    for si in 0..ns {
        let sigma = (si as isize - m as isize) as f64 / (2 * m) as f64;
        for j in 0..n {
            let freq = if j <= (n - 1) / 2 {
                j as f64
            } else {
                j as f64 - n as f64
            };
            t.push(C::from_polar(1.0, 2.0 * PI * freq * sigma / n as f64));
        }
    }
    t
}

impl Engine {
    /// Remove Gibbs ringing from one real line: `out[x]` is the corrected value of `line[x]`.
    pub fn unring_line(&mut self, line: &[f64], out: &mut [f64], p: &ShiftParams) -> Result<()> {
        self.unring_line_impl(line, None, out, None, p)
    }

    /// Remove Gibbs ringing from one complex line given as real and imaginary parts. One shift is chosen
    /// per sample from the oscillation of the complex values (the modulus of the complex difference, so the
    /// choice does not depend on the global phase), and both parts are interpolated with it.
    pub fn unring_line_complex(
        &mut self,
        re: &[f64],
        im: &[f64],
        out_re: &mut [f64],
        out_im: &mut [f64],
        p: &ShiftParams,
    ) -> Result<()> {
        assert_eq!(im.len(), re.len());
        assert_eq!(out_im.len(), re.len());
        self.unring_line_impl(re, Some(im), out_re, Some(out_im), p)
    }

    fn unring_line_impl(
        &mut self,
        line: &[f64],
        line_im: Option<&[f64]>,
        out: &mut [f64],
        mut out_im: Option<&mut [f64]>,
        p: &ShiftParams,
    ) -> Result<()> {
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
        let has_im = line_im.is_some();
        let m = p.nshifts;
        let ns = 2 * m + 1;
        let mut k = std::mem::take(&mut self.kernel);
        let table = k
            .phases
            .remove(&(n, m))
            .unwrap_or_else(|| phase_table(n, m));
        k.shifted.resize(ns * n, 0.0);
        if has_im {
            k.shifted_im.resize(ns * n, 0.0);
        }
        k.absdiff.resize(ns * n, 0.0);
        k.spec.resize(n, C::new(0.0, 0.0));
        k.work.resize(n, C::new(0.0, 0.0));
        let mut spec = std::mem::take(&mut k.spec);
        let mut work = std::mem::take(&mut k.work);

        for (j, c) in spec.iter_mut().enumerate() {
            *c = C::new(line[j], line_im.map_or(0.0, |v| v[j]));
        }
        self.fft1(&mut spec, false);

        let inv_n = 1.0 / n as f64;
        for si in 0..ns {
            let s = si as isize - m as isize;
            let row = &mut k.shifted[si * n..(si + 1) * n];
            if s == 0 {
                row.copy_from_slice(line);
                if let Some(v) = line_im {
                    k.shifted_im[si * n..(si + 1) * n].copy_from_slice(v);
                }
                continue;
            }
            let phases = &table[si * n..(si + 1) * n];
            for (j, w) in work.iter_mut().enumerate() {
                *w = spec[j] * phases[j];
            }
            if n.is_multiple_of(2) {
                work[n / 2] = C::new(0.0, 0.0); // Nyquist zeroed in shifted copies
            }
            self.fft1(&mut work, true);
            for (r, w) in row.iter_mut().zip(work.iter()) {
                *r = w.re * inv_n;
            }
            if has_im {
                for (r, w) in k.shifted_im[si * n..(si + 1) * n]
                    .iter_mut()
                    .zip(work.iter())
                {
                    *r = w.im * inv_n;
                }
            }
        }
        for si in 0..ns {
            let row = &k.shifted[si * n..(si + 1) * n];
            let d = &mut k.absdiff[si * n..(si + 1) * n];
            for j in 0..n - 1 {
                d[j] = (row[j] - row[j + 1]).abs();
            }
            d[n - 1] = (row[n - 1] - row[0]).abs();
            if has_im {
                let ri = &k.shifted_im[si * n..(si + 1) * n];
                for j in 0..n - 1 {
                    d[j] = d[j].hypot((ri[j] - ri[j + 1]).abs());
                }
                d[n - 1] = d[n - 1].hypot((ri[n - 1] - ri[0]).abs());
            }
        }

        // For every sample, the shift with the smallest oscillation on either side. Shifts are visited in
        // ascending order and the plus side before the minus side, with strict `<`, so the first minimum wins.
        let pad = p.max_w + 1;
        k.ext.resize(n + 2 * pad, 0.0);
        k.best.clear();
        k.best.resize(n, f64::INFINITY);
        k.best_si.clear();
        k.best_si.resize(n, m); // s = 0
        for si in 0..ns {
            let d = &k.absdiff[si * n..(si + 1) * n];
            k.ext[..pad].copy_from_slice(&d[n - pad..]);
            k.ext[pad..pad + n].copy_from_slice(d);
            k.ext[pad + n..].copy_from_slice(&d[..pad]);
            for x in 0..n {
                let (mut dp, mut dm) = (0.0, 0.0);
                for kk in p.min_w..=p.max_w {
                    dp += k.ext[pad + x + kk];
                    dm += k.ext[pad + x - kk - 1];
                }
                if dp < k.best[x] {
                    k.best[x] = dp;
                    k.best_si[x] = si;
                }
                if dm < k.best[x] {
                    k.best[x] = dm;
                    k.best_si[x] = si;
                }
            }
        }
        for x in 0..n {
            let best_si = k.best_si[x];
            let s = best_si as isize - m as isize;
            let (a, nb) = if s == 0 {
                (0.0, x)
            } else {
                (
                    s.unsigned_abs() as f64 / (2 * m) as f64,
                    if s > 0 { (x + n - 1) % n } else { (x + 1) % n },
                )
            };
            let row = &k.shifted[best_si * n..(best_si + 1) * n];
            out[x] = if s == 0 {
                row[x]
            } else {
                (1.0 - a) * row[x] + a * row[nb]
            };
            if let Some(o) = out_im.as_deref_mut() {
                let ri = &k.shifted_im[best_si * n..(best_si + 1) * n];
                o[x] = if s == 0 {
                    ri[x]
                } else {
                    (1.0 - a) * ri[x] + a * ri[nb]
                };
            }
        }

        k.spec = spec;
        k.work = work;
        k.phases.insert((n, m), table);
        self.kernel = k;
        Ok(())
    }
}
