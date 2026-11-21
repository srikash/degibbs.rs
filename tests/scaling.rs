//! Integer input with scl_slope/scl_inter is scaled once: on reading, not again on writing.
mod common;
use common::{cli, header};
use ndarray::{ArrayD, IxDyn, ShapeBuilder};
use nifti::writer::WriterOptions;

#[test]
fn scaled_integer_input_is_scaled_once_on_output() {
    let dir = tempfile::tempdir().unwrap();
    let mut v = ArrayD::<i16>::zeros(IxDyn(&[24, 24, 2]).f());
    for ((x, y, z), p) in v.indexed_iter_mut().map(|(i, p)| ((i[0], i[1], i[2]), p)) {
        *p = if (6..18).contains(&x) && (6..18).contains(&y) {
            100
        } else {
            10
        } + z as i16;
    }
    let mut h = header();
    h.scl_slope = 2.0;
    h.scl_inter = 5.0;
    let input = dir.path().join("a.nii");
    WriterOptions::new(&input)
        .reference_header(&h)
        .write_nifti(&v)
        .unwrap();
    let mut bytes = std::fs::read(&input).unwrap();
    bytes[112..116].copy_from_slice(&2.0f32.to_le_bytes());
    bytes[116..120].copy_from_slice(&5.0f32.to_le_bytes());
    std::fs::write(&input, bytes).unwrap();
    let rin = degibbs::io::read(&input).unwrap();
    assert_eq!(rin.data[[0, 0, 0]], 25.0); // 10 * 2 + 5
    let out = dir.path().join("o.nii");
    cli()
        .arg(&input)
        .arg(&out)
        .args(["--pe-axis", "j", "--pf", "1"])
        .assert()
        .success();
    let r = degibbs::io::read(&out).unwrap();
    // the slope was applied on reading and is not applied again: the output header carries none
    assert!(r.header.scl_slope == 1.0 && r.header.scl_inter == 0.0);
    assert!(
        (r.data[[0, 0, 0]] - 25.0).abs() < 1.0,
        "{}",
        r.data[[0, 0, 0]]
    );
}
