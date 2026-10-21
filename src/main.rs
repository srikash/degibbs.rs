use anyhow::{Context, Result};
use clap::Parser;
use degibbs::{
    config::{Requested, resolve},
    io,
    shift::ShiftParams,
    sidecar::{self, Sidecar},
};
use std::path::PathBuf;

/// Gibbs-ringing removal by local subvoxel shifts (Kellner 2016) with partial-Fourier support (Lee 2021).
///
/// The input is a magnitude NIfTI image, or complex data (a complex NIfTI, or a magnitude image plus --phase).
/// The phase-encoding axis and partial-Fourier factor can be read from a BIDS JSON sidecar
/// (PhaseEncodingDirection, PartialFourier, SliceEncodingDirection); flags override it.
///
/// Complex input is assumed to be POCS-reconstructed full-Fourier data and is processed as --pf 1, ignoring any
/// PartialFourier in the sidecar. Zero-filled partial-Fourier data should be given as magnitude.
#[derive(Parser, Debug)]
#[command(version)]
struct Cli {
    /// Input image (.nii or .nii.gz), 2D, 3D or 4D: magnitude, or complex64/complex128
    input: PathBuf,
    /// Output image (.nii or .nii.gz): float32, or complex64 for a complex input
    output: PathBuf,
    /// Phase image in radians, same shape as the input magnitude. Makes the input complex (assumed to be
    /// POCS-reconstructed, processed as --pf 1); the corrected magnitude goes to OUTPUT and the phase to --phase-out
    #[arg(long, requires = "phase_out")]
    phase: Option<PathBuf>,
    /// Where to write the corrected phase (radians, float32) when --phase is given
    #[arg(long, requires = "phase")]
    phase_out: Option<PathBuf>,
    /// Phase-encoding axis in the file's voxel axes: i|j|k or 0|1|2 (a +/- sign is accepted and ignored)
    #[arg(long)]
    pe_axis: Option<String>,
    /// In-plane axis pair containing the PE axis [default: 0,1, or derived from SliceEncodingDirection]
    #[arg(long)]
    axes: Option<String>,
    /// Partial Fourier factor: 1, 7/8, 6/8 or 5/8 (decimals like 0.75 accepted)
    #[arg(long)]
    pf: Option<String>,
    /// BIDS JSON sidecar [default: <input>.json next to the input, if present]
    #[arg(long, conflicts_with = "no_json")]
    json: Option<PathBuf>,
    /// Do not look for a sidecar
    #[arg(long)]
    no_json: bool,
    /// Number of subvoxel shifts M (shifts of s/(2M) voxel, s = -M..M)
    #[arg(long, default_value_t = 20)]
    nshifts: usize,
    /// Start of the oscillation window
    #[arg(long, default_value_t = 1)]
    minw: usize,
    /// End of the oscillation window
    #[arg(long, default_value_t = 3)]
    maxw: usize,
    /// Worker threads [default: all cores]
    #[arg(long)]
    threads: Option<usize>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let sidecar = if cli.no_json {
        None
    } else if let Some(path) = &cli.json {
        Some(Sidecar::from_path(path)?)
    } else {
        let path = sidecar::default_path(&cli.input);
        // An auto-discovered sidecar must not block runs that the flags fully specify.
        if path.exists() {
            Sidecar::from_path(&path)
                .map_err(|e| eprintln!("degibbs: warning: ignoring unreadable sidecar ({e:#})"))
                .ok()
        } else {
            None
        }
    };

    let mut img = io::read(&cli.input)?;
    let mut phase_img = None;
    if let Some(path) = &cli.phase {
        if img.imag.is_some() {
            anyhow::bail!(
                "--phase needs a magnitude input, but {} is already complex",
                cli.input.display()
            );
        }
        let ph = io::read(path)?;
        img.imag = Some(combine_magnitude_and_phase(&mut img.data, &ph, path)?);
        phase_img = Some(ph);
    }
    let req = Requested {
        pe_axis: cli.pe_axis,
        axes: cli.axes,
        pf: cli.pf,
        shift: ShiftParams {
            nshifts: cli.nshifts,
            min_w: cli.minw,
            max_w: cli.maxw,
        },
        complex: img.imag.is_some(),
    };
    let resolved = resolve(&req, sidecar.as_ref(), img.data.ndim())?;
    for n in &resolved.notices {
        eprintln!("degibbs: note: {n}");
    }
    let p = resolved.params;
    eprintln!(
        "degibbs: shape {:?}, PE axis {}, other plane axis {}, PF {:?}",
        img.data.shape(),
        p.pe_axis,
        p.other_axis,
        p.pf
    );

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads.unwrap_or(0))
        .build()
        .context("building thread pool")?;
    pool.install(|| match img.imag.as_mut() {
        Some(im) => degibbs::degibbs_complex(&mut img.data, im, &p),
        None => degibbs::degibbs(&mut img.data, &p),
    })?;
    match (phase_img, &cli.phase_out) {
        (Some(ph), Some(phase_out)) => {
            // Back to magnitude and phase: the magnitude is written like any other output.
            let im = img.imag.take().expect("set when --phase is given");
            let phase = ndarray::Zip::from(&im)
                .and(&img.data)
                .map_collect(|&i, &r| i.atan2(r));
            img.data = ndarray::Zip::from(&img.data)
                .and(&im)
                .map_collect(|&r, &i| r.hypot(i));
            io::write(
                phase_out,
                &io::Image {
                    data: phase,
                    imag: None,
                    header: ph.header,
                    extensions: ph.extensions,
                },
            )?;
            io::write(&cli.output, &img)
        }
        _ => io::write(&cli.output, &img),
    }
}

/// Build the complex image from a magnitude (replaced by its real part) and a phase image in radians.
/// Returns the imaginary part.
fn combine_magnitude_and_phase(
    mag: &mut ndarray::ArrayD<f32>,
    phase: &io::Image,
    path: &std::path::Path,
) -> Result<ndarray::ArrayD<f32>> {
    if phase.imag.is_some() {
        anyhow::bail!("{}: the phase image must not be complex", path.display());
    }
    if phase.data.shape() != mag.shape() {
        anyhow::bail!(
            "{}: phase image shape {:?} differs from the magnitude shape {:?}",
            path.display(),
            phase.data.shape(),
            mag.shape()
        );
    }
    let (lo, hi) = phase
        .data
        .iter()
        .filter(|v| v.is_finite())
        .fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), &v| {
            (lo.min(v), hi.max(v))
        });
    // Only the cosine and sine of the phase are used, so wrapped phase in -pi..pi or 0..2pi is accepted.
    const TOL: f32 = 0.0005;
    let pi = std::f32::consts::PI;
    if lo < -pi - TOL || hi > 2.0 * pi + TOL || hi - lo > 2.0 * pi + TOL {
        anyhow::bail!(
            "{}: phase values range from {lo} to {hi}, which is not wrapped phase in radians (-pi..pi or 0..2pi). \
             Convert to radians first; Siemens phase stored as integers from -4096 to 4095 becomes radians when \
             multiplied by pi/4096",
            path.display()
        );
    }
    let im = ndarray::Zip::from(&*mag)
        .and(&phase.data)
        .map_collect(|&m, &p| (m as f64 * (p as f64).sin()) as f32);
    mag.zip_mut_with(&phase.data, |m, &p| {
        *m = (*m as f64 * (p as f64).cos()) as f32
    });
    Ok(im)
}
