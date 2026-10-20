# degibbs.rs

Gibbs-ringing removal for magnitude MRI. It uses local subvoxel shifts (Kellner et al., 2016) and handles partial-Fourier acquisitions with the RPG extension (Lee et al., 2021).

It is one self-contained binary in pure Rust. NIfTI in, NIfTI out, so it works with MRtrix3, FSL, ANTs or any other pipeline.

## Install

```
cargo install --git https://github.com/srikash/degibbs.rs
```

This puts `degibbs` on your `PATH`. You need Rust 1.85 or newer and nothing else.

## Usage

```
degibbs in.nii.gz out.nii.gz --pe-axis j --pf 6/8
degibbs dwi.nii.gz dwi_degibbs.nii.gz            # reads dwi.json (BIDS sidecar)
```

| flag | meaning |
|---|---|
| `--pe-axis i\|j\|k` | phase-encoding axis in the file's voxel axes (a `+`/`-` sign is accepted and ignored) |
| `--pf 1\|7/8\|6/8\|5/8` | partial-Fourier factor (decimals such as `0.75` work) |
| `--axes 0,1` | in-plane axis pair containing the PE axis (default `0,1`, or taken from `SliceEncodingDirection`) |
| `--json FILE` / `--no-json` | sidecar to read; the default is `<input>.json` next to the input |
| `--nshifts 20 --minw 1 --maxw 3` | subvoxel shifts per voxel and the oscillation window |
| `--threads N` | worker threads (default: all cores) |

Input may be 2D, 3D or 4D. Every plane spanned by `--axes` is corrected, and all other axes are looped over. The output is float32 with the input's header and geometry.

### BIDS sidecar

`PhaseEncodingDirection` (or `PhaseEncodingAxis`) sets `--pe-axis`, `PartialFourier` sets `--pf`, and `SliceEncodingDirection` sets the default `--axes`. Explicit flags win.

If no partial-Fourier factor is available, the tool assumes full Fourier and says so on stderr. This is the most likely way to get a wrong result, so check it for partial-Fourier data. Sidecar axes describe the NIfTI they were written with and become stale if the image is reoriented afterwards (for example with `fslreorient2std`). The tool cannot detect this.

### MRtrix3 and FSL

```
mrconvert dwi.mif dwi.nii.gz -json_export dwi.json -export_grad_fsl dwi.bvec dwi.bval
degibbs dwi.nii.gz dwi_degibbs.nii.gz
mrconvert dwi_degibbs.nii.gz out.mif -fslgrad dwi.bvec dwi.bval
```

## Choosing `--pf`

Pass the partial-Fourier factor the scan was acquired with. The sidecar's `PartialFourier` sets it automatically.

| acquired | `--pf` | expected result |
|---|---|---|
| full Fourier | `1` | the standard Kellner correction |
| 7/8 | `7/8` | removes more ringing than `1` |
| 6/8 | `6/8` | removes more ringing than `1` |
| 5/8 | `5/8` | little or no benefit; ringing remains and the image blurs |

A mismatched value gives a smaller gain, and `5/8` on any other data makes the result worse. The RPG paper recommends 6/8 over 5/8 when you have a choice at acquisition.

## Validation on a Gibbs phantom

The phantom is a 96 × 96 slice of a high-resolution object (ellipses and a rectangle, with a smooth phase), sampled 8× finer. Its k-space is truncated to 96 lines, and for partial Fourier the missing lines on one side of k-space are zeroed. The magnitude is then taken, as for zero-filled scanner data. The truth is the high-resolution object averaged into voxels. Errors are RMSE against that truth, measured more than 4 voxels from any edge, so they show ringing rather than edge blur.

Far-field RMSE (lower is better; the bold value is the best correction for each row):

| acquired | none | `--pf 1` | `--pf 7/8` | `--pf 6/8` | `--pf 5/8` |
|---|---|---|---|---|---|
| full Fourier | 0.649 | **0.392** | 0.405 | **0.392** | 1.228 |
| 7/8 | 0.782 | 0.668 | **0.639** | 0.677 | 1.269 |
| 6/8 | 1.023 | 0.929 | 0.915 | **0.852** | 1.381 |
| 5/8 | 1.878 | 1.835 | 1.841 | **1.826** | 1.891 |

- Full-Fourier data improves by about 40%.
- For 7/8 and 6/8 data, the matching setting removes 17–18% of the error and does better than treating the data as full Fourier.
- The `5/8` setting increases the error on every other data set. On 5/8 data nothing helps much: the best result is a 3% gain, and the matching `5/8` setting is slightly worse than doing nothing.
- Edge error is 14–15.5 in every case. It comes from the resolution lost to k-space truncation, which no correction recovers.
- Zeroing the other side of k-space gives the same picture. With noise at SNR 30 (five seeds) the error floor is about 1.7, the best setting gains about 10% (5% on 5/8 data), the matching setting is still best for 7/8 and 6/8 data, and `5/8` still hurts the other three.

![Phantom images and error maps for each acquired partial-Fourier factor](docs/qc_images.png)

*Left: truth, uncorrected, corrected with the matching factor and corrected as full Fourier. Right: absolute error maps for the same three images, with the far-field and near-edge RMSE in each title.*

![Phase-encode profiles and RMSE heatmaps](docs/qc_profiles.png)

*Phase-encode profiles through the phantom, and heatmaps of far-field RMSE, absolute and relative to uncorrected, for every acquired and correction factor.*

This is one phantom. Real data will differ, so read the table for the ordering, not for the size of the gains.

## Limits

- Magnitude (real-valued) images only. Complex images are not supported.
- Supported partial-Fourier factors are 1, 7/8, 6/8 and 5/8.
- Each plane needs lines longer than `maxw + 2` samples, including the sub-lines the partial-Fourier scheme creates.

## Licence

The NYU Non-Commercial Research License, copied unchanged from [NYU-DiffusionMRI/rpg-degibbs](https://github.com/NYU-DiffusionMRI/rpg-degibbs) (see `LICENSE`): non-commercial research only, not for clinical care, no patent rights granted.

## Citing

The licence (clause 6) requires any publication of research that involves this software to cite:

- **Lee, Novikov, Fieremans.** Removal of partial Fourier-induced Gibbs (RPG) ringing artifacts in MRI. *Magn Reson Med* 86(5):2733-2750, 2021. [doi:10.1002/mrm.28830](https://doi.org/10.1002/mrm.28830). This is the partial-Fourier method implemented here (`--pf 7/8`, `6/8`, `5/8`).

Please also cite the method that RPG extends, which is the full-Fourier part of this tool (`--pf 1`) and the building block of the partial-Fourier pipelines:

- **Kellner, Dhital, Kiselev, Reisert.** Gibbs-ringing artifact removal based on local subvoxel-shifts. *Magn Reson Med* 76(5):1574-1581, 2016. [doi:10.1002/mrm.26054](https://doi.org/10.1002/mrm.26054).

Both references are in [`CITATION.bib`](CITATION.bib) for BibTeX and BibLaTeX.
