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
/// Operates on magnitude NIfTI images. The phase-encoding axis and partial-Fourier factor can be read from a
/// BIDS JSON sidecar (PhaseEncodingDirection, PartialFourier, SliceEncodingDirection); flags override it.
#[derive(Parser, Debug)]
#[command(version)]
struct Cli {
    /// Input image (.nii or .nii.gz), 2D, 3D or 4D
    input: PathBuf,
    /// Output image (.nii or .nii.gz), float32
    output: PathBuf,
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
    /// Number of subvoxel shifts M (shifts of s/(2M) voxel, s = -M..M-1)
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
    let req = Requested {
        pe_axis: cli.pe_axis,
        axes: cli.axes,
        pf: cli.pf,
        shift: ShiftParams {
            nshifts: cli.nshifts,
            min_w: cli.minw,
            max_w: cli.maxw,
        },
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
    pool.install(|| degibbs::degibbs(&mut img.data, &p))?;
    io::write(&cli.output, &img)
}
