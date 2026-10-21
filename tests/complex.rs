//! Complex input: library behaviour and the CLI rules (complex NIfTI, magnitude plus --phase).
mod common;

use assert_cmd::Command;
use common::corrupted_and_truth_sl;
use degibbs::{Params, Pf, ShiftParams, degibbs, degibbs_complex, io};
use ndarray::{ArrayD, IxDyn, ShapeBuilder};
use nifti::NiftiHeader;
use std::path::Path;

const N: usize = 48;

fn params() -> Params {
    Params {
        shift: ShiftParams::default(),
        pf: Pf::Full,
        pe_axis: 1,
        other_axis: 0,
    }
}

/// A Gibbs-corrupted Shepp-Logan magnitude (with a smooth phase) as `[nx, ny, 2]`, plus a smooth phase map.
fn phantom() -> (ArrayD<f32>, ArrayD<f32>) {
    let (plane, _) = corrupted_and_truth_sl(N, N, 8, 1.0, -1, true);
    let mut mag = ArrayD::<f32>::zeros(IxDyn(&[N, N, 2]));
    let mut ph = ArrayD::<f32>::zeros(IxDyn(&[N, N, 2]));
    for x in 0..N {
        for y in 0..N {
            for z in 0..2 {
                mag[[x, y, z]] = (plane.data[y * N + x] * (1.0 + 0.2 * z as f64)) as f32;
                ph[[x, y, z]] =
                    (1.2 * (x as f32 / 9.0).sin() + (y as f32 / 13.0).cos() - 0.5) * 0.8;
            }
        }
    }
    (mag, ph)
}

fn to_complex(mag: &ArrayD<f32>, ph: &ArrayD<f32>, rot: f32) -> (ArrayD<f32>, ArrayD<f32>) {
    let re = ndarray::Zip::from(mag)
        .and(ph)
        .map_collect(|&m, &p| m * (p + rot).cos());
    let im = ndarray::Zip::from(mag)
        .and(ph)
        .map_collect(|&m, &p| m * (p + rot).sin());
    (re, im)
}

fn max_abs(a: &ArrayD<f32>) -> f32 {
    a.iter().fold(0.0, |m, v| m.max(v.abs()))
}

#[test]
fn zero_imaginary_part_reproduces_the_magnitude_path() {
    let (mag, _) = phantom();
    let mut real = mag.clone();
    degibbs(&mut real, &params()).unwrap();
    let (mut re, mut im) = (mag.clone(), ArrayD::<f32>::zeros(mag.raw_dim()));
    degibbs_complex(&mut re, &mut im, &params()).unwrap();
    let worst = ndarray::Zip::from(&real)
        .and(&re)
        .fold(0.0f32, |m, &a, &b| m.max((a - b).abs()));
    let peak = max_abs(&mag);
    eprintln!(
        "zero imaginary part: max difference {worst:.2e} ({:.1e} of the peak)",
        worst / peak
    );
    assert!(worst < 1e-4 * peak);
    assert!(max_abs(&im) < 1e-4 * peak);
}

#[test]
fn a_global_phase_rotation_does_not_change_the_result() {
    let (mag, ph) = phantom();
    let (re0, im0) = to_complex(&mag, &ph, 0.0);
    let (mut a_re, mut a_im) = (re0.clone(), im0.clone());
    degibbs_complex(&mut a_re, &mut a_im, &params()).unwrap();
    for deg in [30.0f32, 45.0, 90.0] {
        let theta = deg.to_radians();
        let (mut re, mut im) = to_complex(&mag, &ph, theta);
        degibbs_complex(&mut re, &mut im, &params()).unwrap();
        // rotate back by -theta
        let (c, s) = (theta.cos(), theta.sin());
        let (mut sum, mut cnt, mut worst) = (0.0f64, 0usize, 0.0f32);
        for ((&r, &i), (&r0, &i0)) in re.iter().zip(im.iter()).zip(a_re.iter().zip(a_im.iter())) {
            let (rb, ib) = (r * c + i * s, -r * s + i * c);
            let d = ((rb - r0).powi(2) + (ib - i0).powi(2)).sqrt();
            sum += (d as f64).powi(2);
            cnt += 1;
            worst = worst.max(d);
        }
        let rms = (sum / cnt as f64).sqrt() as f32;
        let peak = max_abs(&mag);
        eprintln!(
            "rotation {deg} deg: rms difference {:.2e}, max {:.2e} of the peak",
            rms / peak,
            worst / peak
        );
        assert!(rms < 1e-4 * peak, "{deg} deg: {}", rms / peak);
        assert!(worst < 1e-3 * peak, "{deg} deg: max {}", worst / peak);
    }
}

#[test]
fn mismatched_parts_and_partial_fourier_are_rejected() {
    let (mag, _) = phantom();
    let (mut re, mut im) = (mag.clone(), ArrayD::<f32>::zeros(IxDyn(&[N, N, 3])));
    assert!(degibbs_complex(&mut re, &mut im, &params()).is_err());
    let mut im = ArrayD::<f32>::zeros(mag.raw_dim());
    let mut p = params();
    p.pf = Pf::P68;
    assert!(degibbs_complex(&mut re, &mut im, &p).is_err());
}

#[test]
fn non_finite_voxels_in_either_part_stay_put() {
    let (mag, ph) = phantom();
    let (mut re, mut im) = to_complex(&mag, &ph, 0.0);
    re[[3, 4, 0]] = f32::NAN;
    im[[10, 11, 1]] = f32::INFINITY;
    degibbs_complex(&mut re, &mut im, &params()).unwrap();
    assert!(re[[3, 4, 0]].is_nan());
    assert!(im[[10, 11, 1]].is_infinite());
    assert!(re.iter().filter(|v| !v.is_finite()).count() == 1);
    assert!(im.iter().filter(|v| !v.is_finite()).count() == 1);
}

// ---- command line -------------------------------------------------------------------------------------------

fn header() -> NiftiHeader {
    NiftiHeader {
        pixdim: [1.0, 2.0, 2.5, 3.0, 1.0, 0.0, 0.0, 0.0],
        sform_code: 1,
        srow_x: [2.0, 0.0, 0.0, -10.0],
        srow_y: [0.0, 2.5, 0.0, -20.0],
        srow_z: [0.0, 0.0, 3.0, -30.0],
        ..Default::default()
    }
}

fn save(path: &Path, data: ArrayD<f32>, imag: Option<ArrayD<f32>>) {
    let img = io::Image {
        data,
        imag,
        header: header(),
        extensions: None,
    };
    io::write(path, &img).unwrap();
}

fn cli() -> Command {
    Command::cargo_bin("degibbs").unwrap()
}

fn stderr(a: &assert_cmd::assert::Assert) -> String {
    String::from_utf8_lossy(&a.get_output().stderr).into_owned()
}

#[test]
fn complex_nifti_round_trips_as_complex64_and_runs_as_full_fourier() {
    let dir = tempfile::tempdir().unwrap();
    let (mag, ph) = phantom();
    let (re, im) = to_complex(&mag, &ph, 0.0);
    let input = dir.path().join("c.nii.gz");
    save(&input, re.clone(), Some(im.clone()));
    // the sidecar says 0.75 but complex input is processed as --pf 1, with a notice
    std::fs::write(
        dir.path().join("c.json"),
        r#"{"PhaseEncodingDirection":"j","PartialFourier":0.75}"#,
    )
    .unwrap();
    let output = dir.path().join("o.nii.gz");
    let run = cli().arg(&input).arg(&output).assert().success();
    let err = stderr(&run);
    assert!(err.contains("POCS"), "{err}");
    assert!(err.contains("PF Full"), "{err}");

    let got = io::read(&output).unwrap();
    assert_eq!(got.header.datatype, 32);
    assert_eq!(got.header.pixdim, header().pixdim);
    let (mut want_re, mut want_im) = (re, im);
    degibbs_complex(&mut want_re, &mut want_im, &params()).unwrap();
    let got_im = got.imag.expect("complex output");
    for (a, b) in got.data.iter().zip(want_re.iter()) {
        assert!((a - b).abs() < 1e-4);
    }
    for (a, b) in got_im.iter().zip(want_im.iter()) {
        assert!((a - b).abs() < 1e-4);
    }
}

#[test]
fn explicit_partial_fourier_with_complex_input_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let (mag, ph) = phantom();
    let (re, im) = to_complex(&mag, &ph, 0.0);
    let input = dir.path().join("c.nii");
    save(&input, re, Some(im));
    let out = dir.path().join("o.nii");
    let run = cli()
        .args(["--pe-axis", "j", "--pf", "6/8"])
        .arg(&input)
        .arg(&out)
        .assert()
        .failure();
    assert!(stderr(&run).contains("full-Fourier"));
    cli()
        .args(["--pe-axis", "j", "--pf", "1"])
        .arg(&input)
        .arg(&out)
        .assert()
        .success();
}

#[test]
fn magnitude_plus_phase_matches_the_complex_route() {
    let dir = tempfile::tempdir().unwrap();
    let (mag, ph) = phantom();
    let (mp, pp) = (dir.path().join("m.nii.gz"), dir.path().join("p.nii.gz"));
    save(&mp, mag.clone(), None);
    save(&pp, ph.clone(), None);
    std::fs::write(
        dir.path().join("m.json"),
        r#"{"PhaseEncodingDirection":"j","PartialFourier":0.75}"#,
    )
    .unwrap();
    let (om, op) = (dir.path().join("om.nii.gz"), dir.path().join("op.nii.gz"));
    let run = cli()
        .arg(&mp)
        .arg(&om)
        .arg("--phase")
        .arg(&pp)
        .arg("--phase-out")
        .arg(&op)
        .assert()
        .success();
    assert!(stderr(&run).contains("POCS"));

    let (mut re, mut im) = to_complex(&mag, &ph, 0.0);
    degibbs_complex(&mut re, &mut im, &params()).unwrap();
    let got_m = io::read(&om).unwrap();
    let got_p = io::read(&op).unwrap();
    assert!(got_m.imag.is_none() && got_p.imag.is_none());
    let peak = max_abs(&mag);
    for ((&m, &p), (&r, &i)) in got_m
        .data
        .iter()
        .zip(got_p.data.iter())
        .zip(re.iter().zip(im.iter()))
    {
        assert!((m - r.hypot(i)).abs() < 1e-3 * peak);
        // compare phases through the unit vector, which is safe at the +-pi wrap
        if r.hypot(i) > 0.05 * peak {
            assert!(
                ((p.cos() - i.atan2(r).cos()).abs()) < 1e-3,
                "{p} vs {}",
                i.atan2(r)
            );
        }
    }
}

#[test]
fn phase_in_the_wrong_units_gets_a_conversion_hint() {
    let dir = tempfile::tempdir().unwrap();
    let (mag, ph) = phantom();
    let (mp, pp) = (dir.path().join("m.nii"), dir.path().join("p.nii"));
    save(&mp, mag, None);
    save(&pp, ph.mapv(|v| v * 1000.0), None); // like integer-stored Siemens phase
    let run = cli()
        .args(["--pe-axis", "j"])
        .arg(&mp)
        .arg(dir.path().join("o.nii"))
        .arg("--phase")
        .arg(&pp)
        .arg("--phase-out")
        .arg(dir.path().join("po.nii"))
        .assert()
        .failure();
    let err = stderr(&run);
    assert!(err.contains("pi/4096"), "{err}");
}

#[test]
fn phase_flags_come_as_a_pair_and_shapes_must_match() {
    let dir = tempfile::tempdir().unwrap();
    let (mag, ph) = phantom();
    let (mp, pp) = (dir.path().join("m.nii"), dir.path().join("p.nii"));
    save(&mp, mag.clone(), None);
    save(&pp, ph.clone(), None);
    let o = dir.path().join("o.nii");
    cli()
        .args(["--pe-axis", "j"])
        .arg(&mp)
        .arg(&o)
        .arg("--phase")
        .arg(&pp)
        .assert()
        .failure();
    // a phase image of a different shape
    let small = ph.slice(ndarray::s![..N - 2, .., ..]).to_owned().into_dyn();
    let sp = dir.path().join("s.nii");
    save(&sp, small, None);
    let run = cli()
        .args(["--pe-axis", "j"])
        .arg(&mp)
        .arg(&o)
        .arg("--phase")
        .arg(&sp)
        .arg("--phase-out")
        .arg(dir.path().join("po.nii"))
        .assert()
        .failure();
    assert!(stderr(&run).contains("shape"));
    // --phase with an input that is already complex
    let (re, im) = to_complex(&mag, &ph, 0.0);
    let cp = dir.path().join("c.nii");
    save(&cp, re, Some(im));
    cli()
        .args(["--pe-axis", "j"])
        .arg(&cp)
        .arg(&o)
        .arg("--phase")
        .arg(&pp)
        .arg("--phase-out")
        .arg(dir.path().join("po.nii"))
        .assert()
        .failure();
}

#[test]
fn parts_with_different_memory_layouts_are_handled() {
    let (mag, ph) = phantom();
    let (re, im) = to_complex(&mag, &ph, 0.0);
    let (mut c_re, mut c_im) = (re.clone(), im.clone());
    degibbs_complex(&mut c_re, &mut c_im, &params()).unwrap();
    // same values, but the real part in Fortran order and the imaginary part in C order
    let mut f_re = ArrayD::<f32>::zeros(ndarray::IxDyn(re.shape()).f());
    f_re.assign(&re);
    let mut x_im = im.clone();
    degibbs_complex(&mut f_re, &mut x_im, &params()).unwrap();
    for (a, b) in c_re.iter().zip(f_re.iter()) {
        assert!((a - b).abs() < 1e-4);
    }
    for (a, b) in c_im.iter().zip(x_im.iter()) {
        assert!((a - b).abs() < 1e-4);
    }
}

#[test]
fn phase_in_zero_to_two_pi_is_accepted() {
    let dir = tempfile::tempdir().unwrap();
    let (mag, ph) = phantom();
    let (mp, pa, pb) = (
        dir.path().join("m.nii"),
        dir.path().join("a.nii"),
        dir.path().join("b.nii"),
    );
    save(&mp, mag, None);
    save(&pa, ph.clone(), None);
    // the same phase wrapped into 0..2pi
    let two_pi = 2.0 * std::f32::consts::PI;
    save(&pb, ph.mapv(|v| v.rem_euclid(two_pi)), None);
    let run = |p: &Path, out: &str| {
        cli()
            .args(["--pe-axis", "j"])
            .arg(&mp)
            .arg(dir.path().join(out))
            .arg("--phase")
            .arg(p)
            .arg("--phase-out")
            .arg(dir.path().join(format!("p_{out}")))
            .assert()
            .success();
        io::read(&dir.path().join(out)).unwrap().data
    };
    let (a, b) = (run(&pa, "oa.nii"), run(&pb, "ob.nii"));
    for (x, y) in a.iter().zip(b.iter()) {
        assert!((x - y).abs() < 1e-2, "{x} vs {y}");
    }
}
