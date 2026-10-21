mod common;

use approx::assert_abs_diff_eq;
use assert_cmd::Command;
use common::corrupted_and_truth;
use degibbs::{Params, Pf, ShiftParams, degibbs, fft::Engine, io, plane::Plane};
use ndarray::{ArrayD, IxDyn, ShapeBuilder};
use nifti::NiftiHeader;
use std::path::Path;

fn params(pf: Pf, pe: usize, other: usize) -> Params {
    Params {
        shift: ShiftParams::default(),
        pf,
        pe_axis: pe,
        other_axis: other,
    }
}

/// F-ordered `[nx, ny, nz, nt]` volume whose planes are scaled copies of one corrupted phantom.
fn volume(nx: usize, ny: usize, nz: usize, nt: usize, nu: f64) -> (ArrayD<f32>, Plane) {
    let (plane, _) = corrupted_and_truth(nx, ny, 8, nu);
    let mut v = ArrayD::<f32>::zeros(IxDyn(&[nx, ny, nz, nt]).f());
    for ((x, y, z, t), _) in ndarray::Array4::<f32>::zeros((nx, ny, nz, nt)).indexed_iter() {
        v[[x, y, z, t]] = (plane.data[y * nx + x] * (1.0 + 0.1 * z as f64 + 0.3 * t as f64)) as f32;
    }
    (v, plane)
}

fn write_header() -> NiftiHeader {
    NiftiHeader {
        pixdim: [1.0, 2.0, 2.5, 3.0, 1.0, 0.0, 0.0, 0.0],
        sform_code: 1,
        srow_x: [2.0, 0.0, 0.0, -10.0],
        srow_y: [0.0, 2.5, 0.0, -20.0],
        srow_z: [0.0, 0.0, 3.0, -30.0],
        ..Default::default()
    }
}

fn save(path: &Path, data: ArrayD<f32>) {
    let img = io::Image {
        data,
        imag: None,
        header: write_header(),
        extensions: None,
    };
    io::write(path, &img).unwrap();
}

#[test]
fn nifti_roundtrip_keeps_data_and_geometry_for_nii_and_gz() {
    let dir = tempfile::tempdir().unwrap();
    let (v, _) = volume(24, 20, 3, 2, 1.0);
    for name in ["a.nii", "a.nii.gz"] {
        let path = dir.path().join(name);
        save(&path, v.clone());
        let back = io::read(&path).unwrap();
        assert_eq!(back.data.shape(), v.shape());
        assert_eq!(back.data, v);
        assert_eq!(back.header.pixdim, write_header().pixdim);
        assert_eq!(back.header.srow_x, write_header().srow_x);
        assert_eq!(back.header.datatype, 16); // float32
    }
}

#[test]
fn every_plane_of_a_4d_volume_matches_the_plane_pipeline() {
    let (mut v, plane) = volume(40, 36, 3, 2, 0.75);
    let orig = v.clone();
    degibbs(&mut v, &params(Pf::P68, 1, 0)).unwrap();
    let mut e = Engine::new();
    for z in 0..3 {
        for t in 0..2 {
            let mut p = Plane::new(40, 36);
            for y in 0..36 {
                for x in 0..40 {
                    p.data[y * 40 + x] = orig[[x, y, z, t]] as f64;
                }
            }
            e.rpg_plane(&mut p, Pf::P68, &ShiftParams::default())
                .unwrap();
            for y in 0..36 {
                for x in 0..40 {
                    assert_abs_diff_eq!(v[[x, y, z, t]] as f64, p.data[y * 40 + x], epsilon = 1e-4);
                }
            }
        }
    }
    let _ = plane;
}

#[test]
fn pe_on_axis_0_equals_pe_on_axis_1_of_the_permuted_volume() {
    let (v, _) = volume(36, 40, 2, 1, 0.75);
    let mut a = v.clone();
    degibbs(&mut a, &params(Pf::P78, 1, 0)).unwrap();
    let mut b = v.clone().permuted_axes(vec![1, 0, 2, 3]);
    degibbs(&mut b, &params(Pf::P78, 0, 1)).unwrap();
    let b = b.permuted_axes(vec![1, 0, 2, 3]);
    for (x, y) in a.iter().zip(b.iter()) {
        assert_abs_diff_eq!(x, y, epsilon = 1e-5);
    }
}

#[test]
fn pe_along_the_third_axis_works() {
    let (v, _) = volume(36, 40, 2, 1, 1.0);
    // move y to axis 2: shape [nx, nz, ny, nt]
    let mut a = v.clone().permuted_axes(vec![0, 2, 1, 3]);
    degibbs(&mut a, &params(Pf::Full, 2, 0)).unwrap();
    let mut b = v.clone();
    degibbs(&mut b, &params(Pf::Full, 1, 0)).unwrap();
    let a = a.permuted_axes(vec![0, 2, 1, 3]);
    for (x, y) in a.iter().zip(b.iter()) {
        assert_abs_diff_eq!(x, y, epsilon = 1e-5);
    }
}

#[test]
fn invalid_axes_are_rejected() {
    let (mut v, _) = volume(24, 24, 1, 1, 1.0);
    assert!(degibbs(&mut v, &params(Pf::Full, 1, 1)).is_err());
    assert!(degibbs(&mut v, &params(Pf::Full, 7, 0)).is_err());
}

fn cli() -> Command {
    Command::cargo_bin("degibbs").unwrap()
}

#[test]
fn cli_reads_pe_and_pf_from_the_sidecar() {
    let dir = tempfile::tempdir().unwrap();
    let (v, _) = volume(40, 36, 2, 2, 0.75);
    let input = dir.path().join("dwi.nii.gz");
    save(&input, v.clone());
    std::fs::write(
        dir.path().join("dwi.json"),
        r#"{"PhaseEncodingDirection":"j-","PartialFourier":0.75}"#,
    )
    .unwrap();
    let output = dir.path().join("out.nii.gz");
    cli().arg(&input).arg(&output).assert().success();

    let mut want = v;
    degibbs(&mut want, &params(Pf::P68, 1, 0)).unwrap();
    let got = io::read(&output).unwrap();
    assert_eq!(got.header.pixdim, write_header().pixdim);
    for (a, b) in got.data.iter().zip(want.iter()) {
        assert_abs_diff_eq!(a, b, epsilon = 1e-5);
    }
}

#[test]
fn cli_flags_override_the_sidecar_and_missing_pf_prints_a_notice() {
    let dir = tempfile::tempdir().unwrap();
    let (v, _) = volume(36, 40, 1, 1, 1.0);
    let input = dir.path().join("a.nii");
    save(&input, v.clone());
    std::fs::write(
        dir.path().join("a.json"),
        r#"{"PhaseEncodingDirection":"i"}"#,
    )
    .unwrap();
    let output = dir.path().join("o.nii");
    let assert = cli()
        .arg(&input)
        .arg(&output)
        .args(["--pe-axis", "j"])
        .assert()
        .success();
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
    assert!(stderr.contains("assuming full Fourier"), "{stderr}");
    let mut want = v;
    degibbs(&mut want, &params(Pf::Full, 1, 0)).unwrap();
    for (a, b) in io::read(&output).unwrap().data.iter().zip(want.iter()) {
        assert_abs_diff_eq!(a, b, epsilon = 1e-5);
    }
}

#[test]
fn cli_without_any_pe_information_fails_clearly() {
    let dir = tempfile::tempdir().unwrap();
    let (v, _) = volume(24, 24, 1, 1, 1.0);
    let input = dir.path().join("a.nii");
    save(&input, v);
    let out = cli()
        .arg(&input)
        .arg(dir.path().join("o.nii"))
        .assert()
        .failure();
    assert!(String::from_utf8_lossy(&out.get_output().stderr).contains("phase-encoding axis"));
}

#[test]
fn cli_rejects_unsupported_pf() {
    let dir = tempfile::tempdir().unwrap();
    let (v, _) = volume(24, 24, 1, 1, 1.0);
    let input = dir.path().join("a.nii");
    save(&input, v);
    cli()
        .arg(&input)
        .arg(dir.path().join("o.nii"))
        .args(["--pe-axis", "j", "--pf", "0.5"])
        .assert()
        .failure();
}

#[test]
fn non_finite_voxels_stay_put_and_do_not_contaminate_the_plane() {
    let (mut v, _) = volume(40, 36, 1, 1, 1.0);
    v[[0, 0, 0, 0]] = f32::NAN;
    v[[5, 3, 0, 0]] = f32::INFINITY;
    degibbs(&mut v, &params(Pf::Full, 1, 0)).unwrap();
    assert!(v[[0, 0, 0, 0]].is_nan());
    assert!(v[[5, 3, 0, 0]].is_infinite());
    assert_eq!(v.iter().filter(|x| !x.is_finite()).count(), 2);
}

#[test]
fn cli_ignores_a_malformed_auto_sidecar_when_flags_suffice_and_rejects_bad_output_names() {
    let dir = tempfile::tempdir().unwrap();
    let (v, _) = volume(24, 24, 1, 1, 1.0);
    let input = dir.path().join("a.nii");
    save(&input, v);
    std::fs::write(dir.path().join("a.json"), r#"{"PartialFourier":"6/8"}"#).unwrap();
    cli()
        .arg(&input)
        .arg(dir.path().join("o.nii"))
        .args(["--pe-axis", "j"])
        .assert()
        .success();
    cli()
        .arg(&input)
        .arg(dir.path().join("o.img"))
        .args(["--pe-axis", "j"])
        .assert()
        .failure();
}

#[test]
fn absurd_parameters_are_rejected() {
    let (mut v, _) = volume(24, 24, 1, 1, 1.0);
    let mut p = params(Pf::Full, 1, 0);
    p.shift.max_w = usize::MAX;
    assert!(degibbs(&mut v, &p).is_err());
    p.shift = ShiftParams {
        nshifts: 100_000_000,
        ..ShiftParams::default()
    };
    assert!(degibbs(&mut v, &p).is_err());
}
