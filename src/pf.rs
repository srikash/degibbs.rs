//! Partial-Fourier (RPG) pipelines (SPEC section 4). PE axis is y.

use crate::{
    fft::Engine,
    plane::{Axis, Plane},
    shift::ShiftParams,
};
use anyhow::{Result, bail};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pf {
    Full,
    P78,
    P68,
    P58,
}

impl Pf {
    /// Snap a PF fraction (e.g. a BIDS `PartialFourier` value) to a supported mode.
    pub fn from_fraction(f: f64) -> Result<Self> {
        const TOL: f64 = 0.01;
        for (v, pf) in [
            (1.0, Pf::Full),
            (7.0 / 8.0, Pf::P78),
            (6.0 / 8.0, Pf::P68),
            (5.0 / 8.0, Pf::P58),
        ] {
            if (f - v).abs() < TOL {
                return Ok(pf);
            }
        }
        bail!("unsupported partial Fourier fraction {f}: supported values are 1, 7/8, 6/8 and 5/8")
    }

    /// `(p, q)` of the up-sample / split pre-pass: `2*pf - 1 = p/q` (SPEC 4). `None` for the other modes.
    fn pre_pass(self) -> Option<(usize, usize)> {
        match self {
            Pf::P78 => Some((3, 4)),
            Pf::P58 => Some((1, 4)),
            _ => None,
        }
    }

    /// Every 1D line length the pipeline will feed to the kernel, for validation.
    fn line_lengths(self, nx: usize, ny: usize) -> Vec<usize> {
        let mut v = vec![nx, ny];
        match self {
            Pf::P68 => v.extend([ny.div_ceil(2), ny / 2]),
            Pf::P78 | Pf::P58 => {
                let (p, q) = self.pre_pass().unwrap();
                v.extend((0..q).map(|i| (p * ny).saturating_sub(i).div_ceil(q)));
            }
            Pf::Full => {}
        }
        v
    }
}

impl Engine {
    /// Correct one plane for ordinary and partial-Fourier-induced Gibbs ringing.
    pub fn rpg_plane(&mut self, plane: &mut Plane, pf: Pf, p: &ShiftParams) -> Result<()> {
        p.validate()?;
        let (nx, ny) = (plane.nx, plane.ny);
        if let Some(&bad) = pf.line_lengths(nx, ny).iter().find(|&&l| l < p.min_len()) {
            bail!(
                "plane {nx}x{ny} is too small for {pf:?}: a sub-line of length {bad} is shorter than {}",
                p.min_len()
            );
        }
        match pf {
            Pf::Full => self.sushi(plane, p),
            Pf::P68 => {
                // y branch: weighted Gy, unring y, then unring y on even/odd rows separately.
                let (mut jx, mut jy) = self.weighted_pair(plane);
                self.unring_axis(&mut jy, Axis::Y, p)?;
                self.unring_axis_sub(&mut jy, Axis::Y, 0, 2, p)?;
                self.unring_axis_sub(&mut jy, Axis::Y, 1, 2, p)?;
                // x branch: weighted Gx, unring x.
                self.unring_axis(&mut jx, Axis::X, p)?;
                plane.assign_sum(&jy, &jx);
                Ok(())
            }
            Pf::P78 | Pf::P58 => {
                let (up, q) = pf.pre_pass().unwrap();
                // 1. nearest-neighbour upsample along y by `up`.
                let uy = ny * up;
                let mut u = Plane::new(nx, uy);
                for y in 0..uy {
                    let src = (y / up) * nx;
                    u.data[y * nx..(y + 1) * nx].copy_from_slice(&plane.data[src..src + nx]);
                }
                // 2. split into q interleaved sub-images, SuShi each.
                for i in 0..q {
                    let rows: Vec<usize> = (i..uy).step_by(q).collect();
                    let mut sub = Plane::new(nx, rows.len());
                    for (r, &y) in rows.iter().enumerate() {
                        sub.data[r * nx..(r + 1) * nx]
                            .copy_from_slice(&u.data[y * nx..(y + 1) * nx]);
                    }
                    self.sushi(&mut sub, p)?;
                    for (r, &y) in rows.iter().enumerate() {
                        u.data[y * nx..(y + 1) * nx]
                            .copy_from_slice(&sub.data[r * nx..(r + 1) * nx]);
                    }
                }
                // 3. nearest-neighbour downsample (centre phase), 4. ordinary unring y.
                for y in 0..ny {
                    let src = (up * y + up / 2) * nx;
                    plane.data[y * nx..(y + 1) * nx].copy_from_slice(&u.data[src..src + nx]);
                }
                self.unring_axis(plane, Axis::Y, p)
            }
        }
    }
}
