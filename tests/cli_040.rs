//! The 0.4 command line: required partial Fourier factor, --pocs, early option checks, --dry-run,
//! overwrite rules, the output sidecar and the help tiers.
mod common;

use common::{cli, save};
use ndarray::{ArrayD, IxDyn, ShapeBuilder};
use std::path::Path;

fn image(path: &Path) {
    let mut v = ArrayD::<f32>::zeros(IxDyn(&[24, 24, 2]).f());
    for ((x, y, z), p) in ndarray::Array3::<f32>::zeros((24, 24, 2))
        .indexed_iter()
        .map(|(i, v)| (i, *v))
    {
        let _ = p;
        v[[x, y, z]] = if (6..18).contains(&x) && (6..18).contains(&y) {
            100.0
        } else {
            10.0
        } + z as f32;
    }
    save(path, v, None);
}

fn stderr(a: &assert_cmd::assert::Assert) -> String {
    String::from_utf8_lossy(&a.get_output().stderr).into_owned()
}

fn setup() -> (tempfile::TempDir, std::path::PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("a.nii");
    image(&input);
    (dir, input)
}

#[test]
fn a_missing_partial_fourier_factor_is_an_error_naming_the_fixes() {
    let (dir, input) = setup();
    // a sidecar without the key does not help either
    std::fs::write(
        dir.path().join("a.json"),
        r#"{"PhaseEncodingDirection":"j"}"#,
    )
    .unwrap();
    let out = dir.path().join("o.nii");
    let a = cli().arg(&input).arg(&out).assert().failure();
    let e = stderr(&a);
    assert!(e.contains("--pf 1") && e.contains("--pocs"), "{e}");
    assert!(!out.exists());
}

#[test]
fn pocs_ignores_the_sidecar_factor_and_needs_no_pf() {
    let (dir, input) = setup();
    std::fs::write(
        dir.path().join("a.json"),
        r#"{"PhaseEncodingDirection":"j","PartialFourier":0.75}"#,
    )
    .unwrap();
    let out = dir.path().join("o.nii");
    let a = cli().arg(&input).arg(&out).arg("--pocs").assert().success();
    assert!(stderr(&a).contains("ignoring PartialFourier"));
    assert!(out.exists());
    // contradicting --pf is an error, agreeing --pf 1 is fine
    let out2 = dir.path().join("o2.nii");
    cli()
        .arg(&input)
        .arg(&out2)
        .args(["--pocs", "--pf", "6/8"])
        .assert()
        .failure();
    cli()
        .arg(&input)
        .arg(&out2)
        .args(["--pocs", "--pf", "1"])
        .assert()
        .success();
}

#[test]
fn bad_options_fail_before_the_image_is_read() {
    let dir = tempfile::tempdir().unwrap();
    let missing = dir.path().join("nope.nii");
    let out = dir.path().join("o.nii");
    for args in [
        vec!["--pe-axis", "x"],
        vec!["--pe-axis", "j", "--pf", "0.5"],
        vec!["--pe-axis", "j", "--pf", "1", "--nshifts", "5000"],
        vec!["--pe-axis", "j", "--pf", "1", "--minw", "5", "--maxw", "2"],
        vec!["--pe-axis", "j", "--pf", "1", "--axes", "1,1"],
    ] {
        let a = cli().arg(&missing).arg(&out).args(&args).assert().failure();
        // the option error is reported, not "reading nope.nii"
        assert!(!stderr(&a).contains("nope.nii"), "{args:?}: {}", stderr(&a));
    }
}

#[test]
fn dry_run_prints_the_settings_and_writes_nothing() {
    let (dir, input) = setup();
    let out = dir.path().join("o.nii");
    let a = cli()
        .arg(&input)
        .arg(&out)
        .args(["--pe-axis", "j", "--pf", "6/8", "--dry-run"])
        .assert()
        .success();
    let e = stderr(&a);
    assert!(
        e.contains("partial Fourier 6/8") && e.contains("dry run"),
        "{e}"
    );
    assert!(!out.exists() && !dir.path().join("o.json").exists());
    // and it fails where a real run would
    cli()
        .arg(&input)
        .arg(&out)
        .args(["--pe-axis", "j", "--dry-run"])
        .assert()
        .failure();
}

#[test]
fn existing_outputs_need_overwrite_and_the_input_is_never_replaced() {
    let (dir, input) = setup();
    let out = dir.path().join("o.nii");
    let base = ["--pe-axis", "j", "--pf", "1"];
    cli().arg(&input).arg(&out).args(base).assert().success();
    let first = std::fs::read(&out).unwrap();
    std::fs::write(&out, b"keep me").unwrap();
    cli().arg(&input).arg(&out).args(base).assert().failure();
    assert_eq!(std::fs::read(&out).unwrap(), b"keep me");
    cli()
        .arg(&input)
        .arg(&out)
        .args(base)
        .arg("--overwrite")
        .assert()
        .success();
    assert_eq!(std::fs::read(&out).unwrap(), first);
    // a leftover sidecar alone also blocks
    std::fs::remove_file(&out).unwrap();
    cli().arg(&input).arg(&out).args(base).assert().failure();
    // output = input is refused even with --overwrite, also spelled differently
    let again = dir.path().join(".").join("a.nii");
    cli()
        .arg(&input)
        .arg(&again)
        .args(base)
        .arg("--overwrite")
        .assert()
        .failure();
}

#[test]
fn the_output_sidecar_carries_the_input_keys_and_a_degibbs_record() {
    let (dir, input) = setup();
    std::fs::write(
        dir.path().join("a.json"),
        r#"{"PhaseEncodingDirection":"j-","TotalReadoutTime":0.05,"PartialFourier":1}"#,
    )
    .unwrap();
    let out = dir.path().join("o.nii.gz");
    cli()
        .arg(&input)
        .arg(&out)
        .args(["--threads", "2"])
        .assert()
        .success();
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.path().join("o.json")).unwrap()).unwrap();
    assert_eq!(v["PhaseEncodingDirection"], "j-");
    assert_eq!(v["TotalReadoutTime"], 0.05);
    let d = &v["Degibbs"];
    assert_eq!(d["Version"], env!("CARGO_PKG_VERSION"));
    assert_eq!(d["PhaseEncodingAxis"], 1);
    assert_eq!(d["Nshifts"], 20);
    assert_eq!(d["Threads"], 2);
    assert_eq!(d["Source"], "a.nii");
    // the input sidecar is untouched
    assert!(
        !std::fs::read_to_string(dir.path().join("a.json"))
            .unwrap()
            .contains("Degibbs")
    );
}

#[test]
fn without_an_input_sidecar_the_output_sidecar_has_only_the_record() {
    let (dir, input) = setup();
    let out = dir.path().join("o.nii");
    cli()
        .arg(&input)
        .arg(&out)
        .args(["--pe-axis", "j", "--pocs"])
        .assert()
        .success();
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.path().join("o.json")).unwrap()).unwrap();
    assert_eq!(v.as_object().unwrap().len(), 1);
    assert_eq!(v["Degibbs"]["Pocs"], true);
}

#[test]
fn the_output_does_not_depend_on_the_thread_count() {
    let (dir, input) = setup();
    let mut outs = Vec::new();
    for t in ["1", "4"] {
        let out = dir.path().join(format!("o{t}.nii"));
        cli()
            .arg(&input)
            .arg(&out)
            .args(["--pe-axis", "j", "--pf", "1", "--threads", t])
            .assert()
            .success();
        outs.push(std::fs::read(&out).unwrap());
    }
    assert_eq!(outs[0], outs[1]);
}

#[test]
fn short_help_hides_the_expert_options_and_full_help_shows_them() {
    let short = cli().arg("-h").assert().success();
    let s = String::from_utf8_lossy(&short.get_output().stdout).into_owned();
    assert!(s.contains("--pocs") && s.contains("--dry-run") && s.contains("--overwrite"));
    for expert in ["--nshifts", "--minw", "--maxw", "--axes", "--threads"] {
        assert!(!s.contains(expert), "{expert} in short help");
    }
    assert!(s.contains("--full-help"));
    let long = cli().arg("--help").assert().success();
    assert_eq!(long.get_output().stdout, short.get_output().stdout);
    let full = cli().arg("--full-help").assert().success();
    let f = String::from_utf8_lossy(&full.get_output().stdout).into_owned();
    for expert in ["--nshifts", "--minw", "--maxw", "--axes", "--threads"] {
        assert!(f.contains(expert), "{expert} missing from full help");
    }
}

#[test]
fn a_single_2d_slice_gives_the_same_output_at_any_thread_count() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("s.nii");
    let mut v = ArrayD::<f32>::zeros(IxDyn(&[40, 36]).f());
    for ((x, y), p) in v.indexed_iter_mut().map(|(i, p)| ((i[0], i[1]), p)) {
        *p = if (10..30).contains(&x) && (9..27).contains(&y) {
            100.0
        } else {
            5.0
        };
    }
    save(&input, v, None);
    let mut outs = Vec::new();
    for t in ["1", "8"] {
        let out = dir.path().join(format!("o{t}.nii"));
        cli()
            .arg(&input)
            .arg(&out)
            .args(["--pe-axis", "j", "--pf", "6/8", "--threads", t])
            .assert()
            .success();
        outs.push(std::fs::read(&out).unwrap());
    }
    assert_eq!(outs[0], outs[1]);
}
