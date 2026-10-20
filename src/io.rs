//! Thin adapter over the `nifti` crate: read to `f32`, write back reusing the input header.

use anyhow::{Context, Result};
use ndarray::ArrayD;
use nifti::{
    IntoNdArray, NiftiHeader, NiftiObject, ReaderOptions, extension::ExtensionSequence,
    writer::WriterOptions,
};
use std::path::Path;

pub struct Image {
    pub data: ArrayD<f32>,
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
    if matches!(header.datatype, 32 | 128 | 1792 | 2048 | 2304) {
        anyhow::bail!(
            "{}: complex and RGB images are not supported (magnitude data only)",
            path.display()
        );
    }
    let extensions = Some(obj.extensions().clone()).filter(|e| !e.is_empty());
    let data = obj
        .into_volume()
        .into_ndarray::<f32>()
        .with_context(|| format!("decoding {}", path.display()))?;
    Ok(Image {
        data,
        header,
        extensions,
    })
}

/// Write float32 data with the geometry of `header`; `.gz` extension compresses.
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
    opts.write_nifti(&img.data)
        .with_context(|| format!("writing {}", path.display()))
}
