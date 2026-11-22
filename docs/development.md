# Development

Needs Rust 1.88 or newer.

## Tests

```
cargo test --release
```

The tests cover the algorithm (against reference implementations and phantoms), the command line (option checks, `--pocs`, overwrite rules, the output sidecar, the help tiers) and determinism: the output is byte-identical for any `--threads`.

Before sending a change, also run `cargo fmt --check` and `cargo clippy --all-targets -- -D warnings`. CI runs the same.

## Output must not change within a minor version

Within a minor version (for example 0.4.x) the output is bit-identical for the same input, settings and CPU. This is about the image values: `.nii.gz` files compressed with `pigz` differ in file bytes (`--pigz off` removes the difference). A change that moves any output value, however slightly, waits for the next minor version. To check a change against an earlier build:

```
git worktree add ../degibbs-ref v0.4.0 && (cd ../degibbs-ref && cargo build --release)
DEGIBBS_REF_BIN=../degibbs-ref/target/release/degibbs \
    cargo test --release --test reference_build -- --ignored --nocapture
```

This prints, for every partial Fourier setting, how many voxels differ and by how much.

## Phantom validation

`validation/phantom.py` and `validation/compare_reference.py` produce the numbers and figures in [validation.md](validation.md). They need Python with the packages in `validation/requirements.txt`. The commands are at the top of `validation.md`.

## Benchmark

`bench/run.sh` times a run at several thread counts and records wall time and peak memory in a CSV. It is meant for a workstation and is not part of CI.

```
bench/run.sh -t "1 2 4 8" -o results.csv big_dwi.nii.gz
```
