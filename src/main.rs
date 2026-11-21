use anyhow::{Context, Result, bail};
use clap::{CommandFactory, Parser};
use degibbs::{
    Pf,
    config::{Requested, parse_axes, parse_pf, resolve},
    io,
    shift::ShiftParams,
    sidecar::{self, Sidecar, parse_axis},
};
use serde_json::{Value, json};
use std::path::{Path, PathBuf};
use std::time::Instant;

const SHORT_HELP_FOOTER: &str =
    "Expert options (slice plane, search settings, threads): degibbs --full-help";

/// Removes Gibbs ringing: the faint ripples next to sharp edges in MRI images (for example at the edge of the brain
/// or the ventricles). Method: Kellner 2016, with partial-Fourier support after Lee 2021.
///
/// The phase-encoding direction and partial Fourier factor are read from the image's BIDS JSON sidecar when it has
/// them; the options below override or replace it. A magnitude image must have a known partial Fourier factor:
/// degibbs never guesses it. Complex input is processed as full Fourier.
#[derive(Parser, Debug)]
#[command(version, disable_help_flag = true)]
struct Cli {
    /// Image to correct (.nii or .nii.gz), 2D, 3D or 4D: a magnitude image, or a complex image (complex64 or
    /// complex128)
    input: PathBuf,
    /// Where to write the corrected image (.nii or .nii.gz): float32, or complex64 for a complex input. A record of
    /// the settings used is written beside it as <OUTPUT>.json
    output: PathBuf,

    /// Phase-encoding direction: the image axis along which the scanner stepped through k-space line by line. One
    /// of i, j, k (first, second, third image axis) or 0, 1, 2; a trailing sign (j-) is accepted and ignored.
    /// [default: PhaseEncodingDirection from the sidecar]
    #[arg(long, value_name = "AXIS", help_heading = "Acquisition")]
    pe_axis: Option<String>,
    /// Partial Fourier factor: the fraction of k-space the scanner acquired along the phase-encoding direction. One
    /// of 1 (all of it), 7/8, 6/8 or 5/8, or the same as a decimal (0.75). Needed for magnitude images; taken from
    /// PartialFourier in the sidecar if not given, and degibbs stops if it cannot be found
    #[arg(long, value_name = "FACTOR", help_heading = "Acquisition")]
    pf: Option<String>,
    /// The scanner already filled in the missing k-space with POCS reconstruction, so treat the image as full
    /// Fourier. Any PartialFourier in the sidecar is ignored
    #[arg(long, help_heading = "Acquisition")]
    pocs: bool,

    /// JSON sidecar to read the acquisition details from. Use it when the input is a processed copy and only the
    /// original scan has a sidecar. [default: the input's name with .json in place of .nii or .nii.gz, if present]
    #[arg(
        long,
        value_name = "FILE",
        conflicts_with = "no_json",
        help_heading = "Files"
    )]
    json: Option<PathBuf>,
    /// Do not read any sidecar: give --pe-axis and --pf (or --pocs) yourself
    #[arg(long, help_heading = "Files")]
    no_json: bool,
    /// Phase image in radians, same shape as a magnitude input. The magnitude and phase are then corrected together,
    /// as for a complex image (full Fourier data only). Needs --phase-out
    #[arg(
        long,
        value_name = "FILE",
        requires = "phase_out",
        help_heading = "Files"
    )]
    phase: Option<PathBuf>,
    /// Where to write the corrected phase image (radians). Needs --phase
    #[arg(long, value_name = "FILE", requires = "phase", help_heading = "Files")]
    phase_out: Option<PathBuf>,
    /// Replace output files that already exist. Without it degibbs stops rather than overwrite anything. The input
    /// is never overwritten
    #[arg(long, help_heading = "Files")]
    overwrite: bool,

    /// Check the inputs and print the settings that would be used, then stop without correcting anything
    #[arg(long, help_heading = "Other")]
    dry_run: bool,
    /// Show this help
    #[arg(short = 'h', long = "help", help_heading = "Other")]
    help: bool,
    /// Show all options, including the expert ones
    #[arg(long, help_heading = "Other")]
    full_help: bool,

    /// The two image axes that form each slice, for example 0,1. Ringing is removed slice by slice in this plane,
    /// and the plane must contain the phase-encoding axis. [default: the two axes other than
    /// SliceEncodingDirection in the sidecar, otherwise 0,1]
    #[arg(long, value_name = "A,B", help_heading = "Expert")]
    axes: Option<String>,
    /// How finely to search for the best sub-voxel shift: 2M+1 shifts evenly spaced from -1/2 to +1/2 voxel are
    /// tried at each voxel. Larger is finer but slower. Maximum 1000
    #[arg(long, value_name = "M", default_value_t = 20, help_heading = "Expert")]
    nshifts: usize,
    /// Nearest distance, in voxels, used to judge ringing. For each trial shift, degibbs adds up the voxel-to-voxel
    /// changes from --minw to --maxw voxels away on each side and keeps the shift with the least change. 0 is
    /// allowed; must not exceed --maxw
    #[arg(long, value_name = "N", default_value_t = 1, help_heading = "Expert")]
    minw: usize,
    /// Farthest distance, in voxels, used to judge ringing. Maximum 100
    #[arg(long, value_name = "N", default_value_t = 3, help_heading = "Expert")]
    maxw: usize,
    /// Number of CPU threads. The result is identical for any number. [default: all CPUs this process may use
    /// (respects container and Slurm limits, and RAYON_NUM_THREADS)]
    #[arg(long, value_name = "N", help_heading = "Expert")]
    threads: Option<usize>,
}

const EXPERT: [&str; 5] = ["axes", "nshifts", "minw", "maxw", "threads"];

/// `-h`/`--help` show the short help, `--full-help` everything. Handled before parsing so the required inputs
/// are not demanded.
fn print_help_and_exit() {
    let mut full = false;
    let mut short = false;
    for a in std::env::args().skip(1) {
        match a.as_str() {
            "--" => break,
            "--full-help" => full = true,
            "-h" | "--help" => short = true,
            _ => {}
        }
    }
    if !(full || short) {
        return;
    }
    let mut cmd = Cli::command();
    if full {
        let _ = cmd.print_long_help();
    } else {
        for id in EXPERT {
            cmd = cmd.mut_arg(id, |a| a.hide(true));
        }
        cmd = cmd.after_help(SHORT_HELP_FOOTER);
        let _ = cmd.print_long_help();
    }
    std::process::exit(0);
}

fn main() {
    print_help_and_exit();
    if let Err(e) = run(Cli::parse()) {
        eprintln!("Error: {e:#}");
        std::process::exit(1);
    }
}

/// The file `path` resolved as far as it exists, so `./a.nii` and `a.nii` compare equal before either is written.
fn canonical(path: &Path) -> PathBuf {
    if let Ok(p) = path.canonicalize() {
        return p;
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    match (parent.canonicalize(), path.file_name()) {
        (Ok(d), Some(n)) => d.join(n),
        _ => path.to_path_buf(),
    }
}

/// Everything that can be checked from the options alone, before any file is read.
fn check_options(cli: &Cli, shift: &ShiftParams) -> Result<Vec<(PathBuf, PathBuf)>> {
    if let Some(s) = &cli.pe_axis {
        parse_axis(s)?;
    }
    if let Some(s) = &cli.axes {
        parse_axes(s)?;
    }
    let pf = cli.pf.as_deref().map(parse_pf).transpose()?;
    if cli.pocs && pf.is_some_and(|p| p != Pf::Full) {
        bail!(
            "--pocs says the scanner already filled in the missing k-space (full Fourier), which contradicts --pf {}: \
             drop one of them",
            cli.pf.as_deref().unwrap_or_default()
        );
    }
    shift.validate()?;

    // (image, its output sidecar) for every file that will be written
    let mut outs = vec![(cli.output.clone(), sidecar::default_path(&cli.output))];
    if let Some(p) = &cli.phase_out {
        outs.push((p.clone(), sidecar::default_path(p)));
    }
    let mut inputs = vec![&cli.input];
    inputs.extend(cli.phase.as_ref());
    for (img, side) in &outs {
        let name = img.to_string_lossy();
        if !(name.ends_with(".nii") || name.ends_with(".nii.gz")) {
            bail!("output must end in .nii or .nii.gz (got {})", img.display());
        }
        for input in &inputs {
            if canonical(img) == canonical(input) {
                bail!(
                    "the output {} is the input image: degibbs never overwrites its input",
                    img.display()
                );
            }
        }
        if !cli.no_json {
            let given = cli.json.as_ref().cloned();
            let found = given.unwrap_or_else(|| sidecar::default_path(&cli.input));
            if canonical(side) == canonical(&found) {
                bail!(
                    "the output sidecar {} is the input sidecar: degibbs never overwrites its input",
                    side.display()
                );
            }
        }
    }
    if outs.len() == 2 && canonical(&outs[0].0) == canonical(&outs[1].0) {
        bail!("OUTPUT and --phase-out are the same file");
    }
    if !cli.overwrite {
        for p in outs.iter().flat_map(|(a, b)| [a, b]) {
            if p.exists() {
                bail!(
                    "{} already exists: use --overwrite to replace it",
                    p.display()
                );
            }
        }
    }
    Ok(outs)
}

fn pf_name(pf: Pf, pocs: bool) -> &'static str {
    match (pf, pocs) {
        (_, true) => "full Fourier (POCS)",
        (Pf::Full, _) => "1 (full Fourier)",
        (Pf::P78, _) => "7/8",
        (Pf::P68, _) => "6/8",
        (Pf::P58, _) => "5/8",
    }
}

fn run(cli: Cli) -> Result<()> {
    let started = Instant::now();
    let shift = ShiftParams {
        nshifts: cli.nshifts,
        min_w: cli.minw,
        max_w: cli.maxw,
    };
    let outs = check_options(&cli, &shift)?;

    let (sidecar, raw) = if cli.no_json {
        (None, None)
    } else if let Some(path) = &cli.json {
        let (s, r) = Sidecar::from_path_raw(path)?;
        (Some(s), Some(r))
    } else {
        let path = sidecar::default_path(&cli.input);
        // An auto-discovered sidecar must not block runs that the flags fully specify.
        if path.exists() {
            match Sidecar::from_path_raw(&path) {
                Ok((s, r)) => (Some(s), Some(r)),
                Err(e) => {
                    eprintln!("degibbs: warning: ignoring unreadable sidecar ({e:#})");
                    (None, None)
                }
            }
        } else {
            (None, None)
        }
    };

    let (ndim, file_is_complex) = io::peek(&cli.input)?;
    if cli.phase.is_some() && file_is_complex {
        bail!(
            "--phase needs a magnitude input, but {} is already complex",
            cli.input.display()
        );
    }
    let req = Requested {
        pe_axis: cli.pe_axis.clone(),
        axes: cli.axes.clone(),
        pf: cli.pf.clone(),
        pocs: cli.pocs,
        shift,
        complex: file_is_complex || cli.phase.is_some(),
    };
    let resolved = resolve(&req, sidecar.as_ref(), ndim)?;
    for n in &resolved.notices {
        eprintln!("degibbs: note: {n}");
    }
    let p = resolved.params;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads.unwrap_or(0))
        .build()
        .context("building thread pool")?;
    let threads = pool.current_num_threads();
    let pf_text = pf_name(p.pf, cli.pocs);
    eprintln!(
        "degibbs: {}D image, phase-encoding axis {}, slice plane axes {},{}, partial Fourier {pf_text}, \
         shifts {} (M), window {}..{}, {threads} threads",
        ndim,
        p.pe_axis,
        p.pe_axis.min(p.other_axis),
        p.pe_axis.max(p.other_axis),
        cli.nshifts,
        cli.minw,
        cli.maxw
    );
    if cli.dry_run {
        eprintln!("degibbs: dry run: would write {}", cli.output.display());
        return Ok(());
    }

    let mut img = io::read(&cli.input)?;
    let mut phase_img = None;
    if let Some(path) = &cli.phase {
        let ph = io::read(path)?;
        img.imag = Some(combine_magnitude_and_phase(&mut img.data, &ph, path)?);
        phase_img = Some(ph);
    }
    if img.data.ndim() != ndim {
        bail!("{}: the image changed while reading", cli.input.display());
    }
    pool.install(|| match img.imag.as_mut() {
        Some(im) => degibbs::degibbs_complex(&mut img.data, im, &p),
        None => degibbs::degibbs(&mut img.data, &p),
    })?;

    let record = json!({
        "Version": env!("CARGO_PKG_VERSION"),
        "Command": std::env::args().collect::<Vec<_>>().join(" "),
        "Source": cli.input.file_name().map(|n| n.to_string_lossy()).unwrap_or_default(),
        "PhaseEncodingAxis": p.pe_axis,
        "SlicePlaneAxes": [p.pe_axis.min(p.other_axis), p.pe_axis.max(p.other_axis)],
        "PartialFourier": pf_text,
        "Pocs": cli.pocs,
        "Nshifts": cli.nshifts,
        "Minw": cli.minw,
        "Maxw": cli.maxw,
        "Threads": threads,
    });

    if let (Some(mut ph), Some(phase_out)) = (phase_img, &cli.phase_out) {
        {
            // Back to magnitude and phase, in place: the phase image's buffer receives the new phase.
            let im = img.imag.take().expect("set when --phase is given");
            ndarray::Zip::from(&mut ph.data)
                .and(&im)
                .and(&img.data)
                .for_each(|p, &i, &r| *p = i.atan2(r));
            ndarray::Zip::from(&mut img.data)
                .and(&im)
                .for_each(|r, &i| *r = r.hypot(i));
            drop(im);
            io::write(phase_out, &ph)?;
            write_sidecar(&outs[1].1, raw.as_ref(), &record)?;
        }
    }
    io::write(&cli.output, &img)?;
    write_sidecar(&outs[0].1, raw.as_ref(), &record)?;
    eprintln!(
        "degibbs: done in {:.1} s: wrote {}",
        started.elapsed().as_secs_f64(),
        cli.output.display()
    );
    Ok(())
}

/// The input sidecar's keys carried through unchanged, plus one `DeGibbs` record.
fn write_sidecar(path: &Path, input: Option<&Value>, record: &Value) -> Result<()> {
    let mut obj = match input {
        Some(Value::Object(m)) => m.clone(),
        _ => Default::default(),
    };
    if obj.contains_key("DeGibbs") {
        eprintln!("degibbs: note: replacing the DeGibbs record already in the input sidecar");
    }
    obj.insert("DeGibbs".into(), record.clone());
    let text = serde_json::to_string_pretty(&Value::Object(obj))?;
    std::fs::write(path, text + "\n").with_context(|| format!("writing {}", path.display()))
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
