//! Thin adapter over the `nifti` crate: read to `f32` (complex files as real and imaginary parts), write back
//! reusing the input header.

use anyhow::{Context, Result};
use ndarray::ArrayD;
use nifti::{
    IntoNdArray, NiftiHeader, NiftiObject, ReaderOptions, extension::ExtensionSequence,
    writer::WriterOptions,
};
use rustfft::num_complex::{Complex32, Complex64};
use std::path::Path;

pub struct Image {
    /// Real data, or the real part of a complex image.
    pub data: ArrayD<f32>,
    /// Imaginary part of a complex image.
    pub imag: Option<ArrayD<f32>>,
    pub header: NiftiHeader,
    /// Header extensions of the source file, if it had any (carried over on write).
    pub extensions: Option<ExtensionSequence>,
}

/// Read a `.nii`/`.nii.gz`; `scl_slope`/`scl_inter` are applied.
pub fn read(path: &Path) -> Result<Image> {
    let obj = ReaderOptions::new()
        .read_file(path)
        .with_context(|| format!("reading {}", path.display()))?;
    let header = obj.header().clone();
    // NIfTI codes: 32 complex64, 128 RGB24, 1792 complex128, 2048 complex256, 2304 RGBA32.
    if matches!(header.datatype, 128 | 2048 | 2304) {
        anyhow::bail!(
            "{}: RGB and 256-bit complex images are not supported",
            path.display()
        );
    }
    let is_complex = matches!(header.datatype, 32 | 1792);
    if is_complex && header.scl_slope != 0.0 && (header.scl_slope != 1.0 || header.scl_inter != 0.0)
    {
        anyhow::bail!(
            "{}: complex images with scl_slope/scl_inter scaling are not supported",
            path.display()
        );
    }
    let extensions = Some(obj.extensions().clone()).filter(|e| !e.is_empty());
    let decode = |e: nifti::NiftiError| anyhow::anyhow!("decoding {}: {e}", path.display());
    let (data, imag) = match header.datatype {
        32 => {
            let c = obj
                .into_volume()
                .into_ndarray::<Complex32>()
                .map_err(decode)?;
            (c.mapv(|v| v.re), Some(c.mapv(|v| v.im)))
        }
        1792 => {
            let c = obj
                .into_volume()
                .into_ndarray::<Complex64>()
                .map_err(decode)?;
            (c.mapv(|v| v.re as f32), Some(c.mapv(|v| v.im as f32)))
        }
        _ => (
            obj.into_volume().into_ndarray::<f32>().map_err(decode)?,
            None,
        ),
    };
    Ok(Image {
        data,
        imag,
        header,
        extensions,
    })
}

/// Write float32 data (complex64 if `imag` is set) with the geometry of `header`; `.gz` extension compresses.
pub fn write(path: &Path, img: &Image) -> Result<()> {
    let name = path.to_string_lossy();
    if !(name.ends_with(".nii") || name.ends_with(".nii.gz")) {
        anyhow::bail!(
            "output must end in .nii or .nii.gz (got {})",
            path.display()
        );
    }
    // The input's display range no longer describes the corrected data.
    let header = NiftiHeader {
        cal_min: 0.0,
        cal_max: 0.0,
        ..img.header.clone()
    };
    let mut opts = WriterOptions::new(path);
    opts = opts.reference_header(&header);
    if let Some(ext) = &img.extensions {
        opts = opts.with_extensions(ext.clone());
    }
    match &img.imag {
        None => opts.write_nifti(&img.data),
        Some(im) => {
            let z = ndarray::Zip::from(&img.data)
                .and(im)
                .map_collect(|&r, &i| Complex32::new(r, i));
            opts.write_nifti(&z)
        }
    }
    .with_context(|| format!("writing {}", path.display()))
}
