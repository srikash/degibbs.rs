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
use plane::Plane;
use rayon::prelude::*;

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

/// Correct every `(other_axis, pe_axis)` plane of `vol` (2D or higher), in place.
/// All remaining axes are looped over. Planes are processed in parallel on the current rayon pool.
pub fn degibbs(vol: &mut ArrayD<f32>, p: &Params) -> Result<()> {
    let nd = vol.ndim();
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
    if vol.as_slice_memory_order().is_none() || vol.strides().iter().any(|&s| s <= 0) {
        *vol = vol.as_standard_layout().into_owned();
    }
    let shape = vol.shape().to_vec();
    let strides: Vec<usize> = vol.strides().iter().map(|&s| s as usize).collect();
    let (nx, ny) = (shape[p.other_axis], shape[p.pe_axis]);
    let (sx, sy) = (strides[p.other_axis], strides[p.pe_axis]);

    // Base offsets of every plane: all combinations of the remaining axes.
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
