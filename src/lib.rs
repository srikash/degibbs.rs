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

thread_local! {
    /// One FFT engine per worker thread, so plans and scratch buffers are reused across planes and chunks.
    static ENGINE: std::cell::RefCell<fft::Engine> = std::cell::RefCell::new(fft::Engine::new());
}

/// A raw view of an array that worker threads write through. Each plane touches only its own elements
/// (the base offsets of two planes never produce the same index), so the accesses never overlap.
#[derive(Clone, Copy)]
struct Shared(*mut f32);
unsafe impl Send for Shared {}
unsafe impl Sync for Shared {}
impl Shared {
    /// SAFETY: `i` must be inside the array and not written by another thread at the same time.
    unsafe fn get(self, i: usize) -> f32 {
        unsafe { *self.0.add(i) }
    }
    /// SAFETY: as for `get`.
    unsafe fn set(self, i: usize, v: f32) {
        unsafe { *self.0.add(i) = v }
    }
}

/// Run `f` for every plane base: planes in parallel when there are at least as many as threads, otherwise one
/// plane at a time with the lines of each plane spread over the threads. The result does not depend on which.
fn for_each_plane<F>(bases: &[usize], f: F) -> Result<()>
where
    F: Fn(&mut fft::Engine, usize) -> Result<()> + Sync,
{
    if bases.len() < rayon::current_num_threads() {
        let mut engine = fft::Engine::new();
        engine.par_lines = true;
        bases.iter().try_for_each(|&b| f(&mut engine, b))
    } else {
        bases
            .par_iter()
            .try_for_each(|&b| ENGINE.with(|e| f(&mut e.borrow_mut(), b)))
    }
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
    if vols[0].shape()[p.pe_axis] == 0 || vols[0].shape()[p.other_axis] == 0 {
        bail!(
            "the plane axes must not be empty (shape {:?})",
            vols[0].shape()
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
    let data = Shared(
        vol.as_slice_memory_order_mut()
            .expect("made contiguous above")
            .as_mut_ptr(),
    );
    for_each_plane(&bases, |engine, b| {
        let mut pl = Plane::new(nx, ny);
        for y in 0..ny {
            for x in 0..nx {
                // SAFETY: this plane's own element, see `Shared`.
                let v = unsafe { data.get(b + y * sy + x * sx) };
                // NaN/Inf (e.g. outside a mask) are processed as 0 and restored below.
                pl.data[y * nx + x] = if v.is_finite() { v as f64 } else { 0.0 };
            }
        }
        engine.rpg_plane(&mut pl, p.pf, &p.shift)?;
        for y in 0..ny {
            for x in 0..nx {
                let i = b + y * sy + x * sx;
                // SAFETY: as above.
                unsafe {
                    if data.get(i).is_finite() {
                        data.set(i, pl.data[y * nx + x] as f32);
                    }
                }
            }
        }
        Ok(())
    })
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
    let dre = Shared(
        re.as_slice_memory_order_mut()
            .expect("made contiguous above")
            .as_mut_ptr(),
    );
    let dim = Shared(
        im.as_slice_memory_order_mut()
            .expect("made contiguous above")
            .as_mut_ptr(),
    );
    for_each_plane(&bases, |engine, b| {
        let mut pl = CPlane::new(nx, ny);
        for y in 0..ny {
            for x in 0..nx {
                let i = b + y * sy + x * sx;
                // SAFETY: this plane's own elements, see `Shared`.
                let (r, m) = unsafe { (dre.get(i), dim.get(i)) };
                // A voxel with a NaN/Inf in either part is processed as 0 and restored below.
                pl.data[y * nx + x] = if r.is_finite() && m.is_finite() {
                    C::new(r as f64, m as f64)
                } else {
                    C::new(0.0, 0.0)
                };
            }
        }
        engine.rpg_plane_complex(&mut pl, &p.shift)?;
        for y in 0..ny {
            for x in 0..nx {
                let i = b + y * sy + x * sx;
                // SAFETY: as above.
                unsafe {
                    if dre.get(i).is_finite() && dim.get(i).is_finite() {
                        let v = pl.data[y * nx + x];
                        dre.set(i, v.re as f32);
                        dim.set(i, v.im as f32);
                    }
                }
            }
        }
        Ok(())
    })
}
