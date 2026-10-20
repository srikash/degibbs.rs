# degibbs.rs

Gibbs-ringing removal for magnitude MRI in pure Rust: local subvoxel shifts (Kellner et al., MRM 2016) with the
partial-Fourier extension "RPG" (Lee et al., MRM 2021). A single self-contained binary: NIfTI in, NIfTI out, so it fits in
MRtrix3, FSL, ANTs or any other pipeline. The algorithm is specified in [`SPEC.md`](SPEC.md); the papers are the
only reference. Numeric results may differ slightly from other implementations of the same papers, and the
choices the papers leave open are listed in `SPEC.md`.

## Usage

```
degibbs in.nii.gz out.nii.gz --pe-axis j --pf 6/8
degibbs dwi.nii.gz dwi_degibbs.nii.gz            # reads dwi.json (BIDS sidecar)
```

| flag | meaning |
|---|---|
| `--pe-axis i\|j\|k` | phase-encoding axis in the file's voxel axes (a `+`/`-` sign is accepted and ignored) |
| `--pf 1\|7/8\|6/8\|5/8` | partial Fourier factor (decimals such as `0.75` work) |
| `--axes 0,1` | in-plane axis pair containing the PE axis (default `0,1`, or from `SliceEncodingDirection`) |
| `--json FILE` / `--no-json` | sidecar to read; default is `<input>.json` next to the input |
| `--nshifts 20 --minw 1 --maxw 3` | subvoxel shifts per voxel and the oscillation window |
| `--threads N` | worker threads (default: all cores) |

Input may be 2D, 3D or 4D; every plane spanned by `--axes` is corrected, and all other axes are looped over.
The output is float32 with the input's header and geometry.

### BIDS sidecar

`PhaseEncodingDirection` (or `PhaseEncodingAxis`) gives `--pe-axis`, `PartialFourier` gives `--pf`, and `SliceEncodingDirection` gives the
default `--axes`. Explicit flags win. If no partial Fourier factor is available the tool assumes full Fourier and
says so on stderr; that is the main way to get a wrong result, so check it for partial-Fourier data. Sidecar axes
describe the NIfTI they were written with and become stale if the image is reoriented afterwards (for example
with `fslreorient2std`); the tool cannot detect this.

### MRtrix3 / FSL

```
mrconvert dwi.mif dwi.nii.gz -json_export dwi.json -export_grad_fsl dwi.bvec dwi.bval
degibbs dwi.nii.gz dwi_degibbs.nii.gz
mrconvert dwi_degibbs.nii.gz out.mif -fslgrad dwi.bvec dwi.bval
```

## Scope and limits

- Magnitude (real-valued) images only. Complex data is not supported.
- Supported partial Fourier factors: 1, 7/8, 6/8, 5/8. As the RPG paper reports, 5/8 keeps residual ringing and blurs;
  use 6/8 where there is a choice.
- Each plane needs lines longer than `maxw + 2` samples, including the sub-lines the partial-Fourier scheme creates.

## Build and test

```
cargo build --release     # target/release/degibbs
cargo test --release
```

Tests need no external tools. Correctness is checked against analytic ground truth (phantoms truncated in k-space,
compared with their block-averaged high-resolution versions), a naive reference implementation of the 1D kernel,
and metamorphic checks (transposed planes, permuted axes).

## Licence

The NYU Non-Commercial Research License, copied unchanged from [NYU-DiffusionMRI/rpg-degibbs](https://github.com/NYU-DiffusionMRI/rpg-degibbs) (see `LICENSE`): non-commercial research only, not for clinical care, no patent rights granted. The licence requires publications to cite Lee et al., MRM 2021 (doi:10.1002/mrm.28830). Please also cite Kellner et al., MRM 2016 (doi:10.1002/mrm.26054), whose method it builds on.
