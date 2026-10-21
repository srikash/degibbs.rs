//! Gibbs-ringing removal by local subvoxel shifts, with partial-Fourier (RPG) support.
//! The algorithm follows Kellner et al. (2016) and Lee et al. (2021); see `CITATION.bib`.

pub mod config;
pub mod fft;
pub mod io;
pub mod pf;
pub mod plane;
pub mod shift;
pub mod sidecar;

use anyhow::{Result, bail};
use ndarray::ArrayD;
use plane::{CPlane, Plane};
use rayon::prelude::*;
use rustfft::num_complex::Complex64 as C;

pub use pf::Pf;
pub use shift::ShiftParams;

/// What to correct and along which axes.
#[derive(Clone, Copy, Debug)]
pub struct Params {
    pub shift: ShiftParams,
    pub pf: Pf,
    /// Phase-encoding axis (the plane's y axis).
    pub pe_axis: usize,
    /// The other in-plane axis (the plane's x axis).
    pub other_axis: usize,
}

/// Upper bound on planes held (as f64) at once; chunks are also capped by `PLANE_BUDGET_BYTES`.
const MAX_PLANES_PER_CHUNK: usize = 512;
const PLANE_BUDGET_BYTES: usize = 256 << 20;

thread_local! {
    /// One FFT engine per worker thread, so plans and scratch buffers are reused across planes and chunks.
    static ENGINE: std::cell::RefCell<fft::Engine> = std::cell::RefCell::new(fft::Engine::new());
}

struct Layout {
    nx: usize,
    ny: usize,
    sx: usize,
    sy: usize,
    /// Base offsets of every plane: all combinations of the remaining axes.
    bases: Vec<usize>,
}

/// Check the axes and parameters, and make `vols` contiguous (all of equal shape).
fn prepare(vols: &mut [&mut ArrayD<f32>], p: &Params) -> Result<Layout> {
    let nd = vols[0].ndim();
    if nd < 2 {
        bail!("need at least a 2D image, got {nd}D");
    }
    if p.pe_axis >= nd || p.other_axis >= nd || p.pe_axis == p.other_axis {
        bail!(
            "invalid plane axes (pe {}, other {}) for a {nd}D image",
            p.pe_axis,
            p.other_axis
        );
    }
    p.shift.validate()?;
    // Every array is indexed with the strides of the first, so all must share one memory layout.
    let same_layout = vols.iter().all(|v| v.strides() == vols[0].strides());
    for v in vols.iter_mut() {
        if !same_layout
            || v.as_slice_memory_order().is_none()
            || v.strides().iter().any(|&s| s <= 0)
        {
            **v = v.as_standard_layout().into_owned();
        }
    }
    let shape = vols[0].shape().to_vec();
    let strides: Vec<usize> = vols[0].strides().iter().map(|&s| s as usize).collect();
    let rest: Vec<usize> = (0..nd)
        .filter(|&a| a != p.pe_axis && a != p.other_axis)
        .collect();
    let mut bases = vec![0usize];
    for &a in &rest {
        let (n, st) = (shape[a], strides[a]);
        bases = bases
            .iter()
            .flat_map(|&b| (0..n).map(move |i| b + i * st))
            .collect();
    }
    Ok(Layout {
        nx: shape[p.other_axis],
        ny: shape[p.pe_axis],
        sx: strides[p.other_axis],
        sy: strides[p.pe_axis],
        bases,
    })
}

/// Correct every `(other_axis, pe_axis)` plane of `vol` (2D or higher), in place.
/// All remaining axes are looped over. Planes are processed in parallel on the current rayon pool.
pub fn degibbs(vol: &mut ArrayD<f32>, p: &Params) -> Result<()> {
    let Layout {
        nx,
        ny,
        sx,
        sy,
        bases,
    } = prepare(&mut [vol], p)?;
    let data = vol
        .as_slice_memory_order_mut()
        .expect("made contiguous above");
    let per_chunk = (PLANE_BUDGET_BYTES / (nx * ny * 8))
        .clamp(rayon::current_num_threads(), MAX_PLANES_PER_CHUNK);
    for chunk in bases.chunks(per_chunk) {
        let mut planes: Vec<Plane> = chunk
            .iter()
            .map(|&b| {
                let mut pl = Plane::new(nx, ny);
                for y in 0..ny {
                    for x in 0..nx {
                        let v = data[b + y * sy + x * sx];
                        // NaN/Inf (e.g. outside a mask) are processed as 0 and restored below.
                        pl.data[y * nx + x] = if v.is_finite() { v as f64 } else { 0.0 };
                    }
                }
                pl
            })
            .collect();
        planes
            .par_iter_mut()
            .try_for_each(|pl| ENGINE.with(|e| e.borrow_mut().rpg_plane(pl, p.pf, &p.shift)))?;
        for (&b, pl) in chunk.iter().zip(&planes) {
            for y in 0..ny {
                for x in 0..nx {
                    let i = b + y * sy + x * sx;
                    if data[i].is_finite() {
                        data[i] = pl.data[y * nx + x] as f32;
                    }
                }
            }
        }
    }
    Ok(())
}

/// Complex counterpart of [`degibbs`] for an image given as real and imaginary parts of equal shape.
/// Only full-Fourier data is supported (`p.pf` must be [`Pf::Full`]): complex input is taken to be
/// POCS-reconstructed, so its k-space again ends symmetrically.
pub fn degibbs_complex(re: &mut ArrayD<f32>, im: &mut ArrayD<f32>, p: &Params) -> Result<()> {
    if re.shape() != im.shape() {
        bail!(
            "real and imaginary parts differ in shape ({:?} and {:?})",
            re.shape(),
            im.shape()
        );
    }
    if p.pf != Pf::Full {
        bail!("complex images are supported for full-Fourier data only (--pf 1)");
    }
    let Layout {
        nx,
        ny,
        sx,
        sy,
        bases,
    } = prepare(&mut [re, im], p)?;
    let dre = re
        .as_slice_memory_order_mut()
        .expect("made contiguous above");
    let dim = im
        .as_slice_memory_order_mut()
        .expect("made contiguous above");
    let per_chunk = (PLANE_BUDGET_BYTES / (nx * ny * 16))
        .clamp(rayon::current_num_threads(), MAX_PLANES_PER_CHUNK);
    for chunk in bases.chunks(per_chunk) {
        let mut planes: Vec<CPlane> = chunk
            .iter()
            .map(|&b| {
                let mut pl = CPlane::new(nx, ny);
                for y in 0..ny {
                    for x in 0..nx {
                        let i = b + y * sy + x * sx;
                        // A voxel with a NaN/Inf in either part is processed as 0 and restored below.
                        pl.data[y * nx + x] = if dre[i].is_finite() && dim[i].is_finite() {
                            C::new(dre[i] as f64, dim[i] as f64)
                        } else {
                            C::new(0.0, 0.0)
                        };
                    }
                }
                pl
            })
            .collect();
        planes
            .par_iter_mut()
            .try_for_each(|pl| ENGINE.with(|e| e.borrow_mut().rpg_plane_complex(pl, &p.shift)))?;
        for (&b, pl) in chunk.iter().zip(&planes) {
            for y in 0..ny {
                for x in 0..nx {
                    let i = b + y * sy + x * sx;
                    if dre[i].is_finite() && dim[i].is_finite() {
                        let v = pl.data[y * nx + x];
                        (dre[i], dim[i]) = (v.re as f32, v.im as f32);
                    }
                }
            }
        }
    }
    Ok(())
}
