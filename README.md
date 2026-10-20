# degibbs.rs

Gibbs-ringing removal for magnitude MRI, in pure Rust.

- **Method:** local subvoxel shifts (Kellner et al., 2016), with the RPG extension for partial-Fourier acquisitions (Lee et al., 2021).
- **Input and output:** NIfTI in, NIfTI out. Works with MRtrix3, FSL, ANTs or any other pipeline.
- **Partial Fourier:** full Fourier, 7/8, 6/8 and 5/8.
- **Binary:** one self-contained executable, no other dependencies.
- **Licence:** non-commercial research only, not for clinical care.

## Contents

- [Quick start](#quick-start)
- [Usage](#usage)
  - [Options](#options)
  - [BIDS sidecar](#bids-sidecar)
  - [MRtrix3 and FSL](#mrtrix3-and-fsl)
- [Containers (Docker and Apptainer)](#containers-docker-and-apptainer)
- [Choosing the partial-Fourier factor](#choosing-the-partial-fourier-factor)
- [Validation on a Gibbs phantom](#validation-on-a-gibbs-phantom)
  - [Setup](#setup)
  - [Results](#results)
  - [Findings](#findings)
  - [Figures](#figures)
- [Limits](#limits)
- [Licence and citing](#licence-and-citing)

## Quick start

Download the static Linux binary (x86-64) from the [v0.1.5 release](https://github.com/srikash/degibbs.rs/releases/tag/v0.1.5). It runs on Rocky 9 and 10, Ubuntu 22.04 and 24.04, and Fedora 44. Use `gh release download v0.1.5 -R srikash/degibbs.rs` while the repository is private.

```
tar xzf degibbs-v0.1.5-x86_64-linux-musl.tar.gz
```

Or build from source (needs Rust 1.85 or newer; this puts `degibbs` on your `PATH`):

```
cargo install --git https://github.com/srikash/degibbs.rs
```

Run it:

```
degibbs in.nii.gz out.nii.gz --pe-axis j --pf 6/8
degibbs dwi.nii.gz dwi_degibbs.nii.gz            # reads dwi.json (BIDS sidecar)
```

## Usage

```
degibbs <input> <output> [options]
```

Input may be 2D, 3D or 4D. Every plane spanned by `--axes` is corrected, and all other axes are looped over. The output is float32 with the input's header and geometry.

### Options

| flag | meaning | default |
|---|---|---|
| `--pe-axis i\|j\|k` (or `0\|1\|2`) | phase-encoding axis in the file's voxel axes (a `+`/`-` sign is accepted and ignored) | from the sidecar; required otherwise |
| `--pf 1\|7/8\|6/8\|5/8` | partial-Fourier factor (decimals such as `0.75` work) | from the sidecar, else `1` |
| `--axes 0,1` | in-plane axis pair containing the PE axis | `0,1`, or from `SliceEncodingDirection` |
| `--json FILE` | BIDS sidecar to read | `<input>.json` next to the input |
| `--no-json` | do not look for a sidecar | off |
| `--nshifts N` | subvoxel shifts per voxel | `20` |
| `--minw N`, `--maxw N` | start and end of the oscillation window | `1`, `3` |
| `--threads N` | worker threads | all cores |

### BIDS sidecar

Keys the tool reads:

| key | sets |
|---|---|
| `PhaseEncodingDirection` (or `PhaseEncodingAxis`) | `--pe-axis` |
| `PartialFourier` | `--pf` |
| `SliceEncodingDirection` | the default `--axes` |

Explicit flags override the sidecar.

> [!WARNING]
> If no partial-Fourier factor is available, the tool assumes full Fourier and says so on stderr. This is the most likely way to get a wrong result. Check that message for partial-Fourier data.

> [!NOTE]
> Sidecar axes describe the NIfTI they were written with. They become stale if the image is reoriented afterwards (for example with `fslreorient2std`), and the tool cannot detect this.

### MRtrix3 and FSL

```
mrconvert dwi.mif dwi.nii.gz -json_export dwi.json -export_grad_fsl dwi.bvec dwi.bval
degibbs dwi.nii.gz dwi_degibbs.nii.gz
mrconvert dwi_degibbs.nii.gz out.mif -fslgrad dwi.bvec dwi.bval
```

## Containers (Docker and Apptainer)

A multi-architecture image (linux/amd64 and linux/arm64, about 18 MB) is published to the GitHub Container Registry as `ghcr.io/srikash/degibbs.rs`. Tags: `latest` (tracks `main`), `main`, `0.1.5` (the current release), and `sha-<commit>`.

The package is private while the repository is private. Log in first with a GitHub token that has the `read:packages` scope.

**Docker**

```
echo "$GITHUB_TOKEN" | docker login ghcr.io -u <github-user> --password-stdin
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data ghcr.io/srikash/degibbs.rs \
    in.nii.gz out.nii.gz --pe-axis j --pf 6/8
```

The container's working directory is `/data`, so relative paths resolve inside the mounted folder. Passing `-u` keeps the output files owned by you.

**Apptainer (HPC)**

```
export APPTAINER_DOCKER_USERNAME=<github-user>
export APPTAINER_DOCKER_PASSWORD=<token>
apptainer pull degibbs.sif docker://ghcr.io/srikash/degibbs.rs:latest

apptainer run degibbs.sif in.nii.gz out.nii.gz --pe-axis j --pf 6/8
```

Apptainer binds your current directory and home by default. Use `--bind` for other paths such as a project or scratch filesystem. The tool uses the CPUs your job is allowed to use, so set `--cpus-per-task` in Slurm, or pass `--threads N`.

## Choosing the partial-Fourier factor

Pass the factor the scan was acquired with. The sidecar's `PartialFourier` sets it automatically.

| acquired | `--pf` | expected result |
|---|---|---|
| full Fourier | `1` | the standard Kellner correction |
| 7/8 | `7/8` | removes more ringing than `1`; edges lose some accuracy |
| 6/8 | `6/8` | removes more ringing than `1`; edges lose some accuracy |
| 5/8 | `5/8` | modest benefit (about 10% in the phantom); strong ringing remains and edges blur |

A mismatched value gives a smaller gain, and `5/8` on any other data makes the result worse. The RPG paper recommends 6/8 over 5/8 when you have a choice at acquisition.

## Validation on a Gibbs phantom

### Setup

- The phantom is the modified Shepp-Logan head (Toft, 1996), scaled so the skull is 100 and the brain 20. It is multiplied by a smooth phase map, as in real data.
- The RPG paper reports the same behaviour on its 90 × 90 Shepp-Logan phantom: robust removal for 7/8 and 6/8, and blurring with remaining ringing for 5/8.
- A 96 × 96 slice is sampled from the object at 8× the resolution, so k-space truncation produces true Gibbs ringing.
- For partial Fourier, the missing lines on one side of k-space are zeroed.
- The magnitude is then taken, as for zero-filled scanner data.
- The truth is the high-resolution object averaged into voxels.
- Far-field error is the RMSE against that truth, measured more than 4 voxels from any edge. It shows ringing rather than edge blur. Near-edge error is measured within 2 voxels of an edge.

### Results

Far-field RMSE (lower is better; bold marks the best correction in each row, and both when they are within 2%):

| acquired | none | `--pf 1` | `--pf 7/8` | `--pf 6/8` | `--pf 5/8` |
|---|---|---|---|---|---|
| full Fourier | 0.779 | **0.477** | 0.486 | **0.472** | 0.841 |
| 7/8 | 0.771 | 0.566 | **0.524** | 0.561 | 0.841 |
| 6/8 | 0.852 | 0.687 | 0.673 | **0.622** | 0.886 |
| 5/8 | 1.368 | 1.279 | 1.279 | 1.261 | **1.236** |

### Findings

- Full-Fourier data improves by 39%.
- For 7/8 and 6/8 data, the matching setting removes 32% and 27% of the far-field error. That is 7% and 9% better than treating the data as full Fourier.
- The `5/8` setting increases the far-field error on every other data set, by 4% to 9%.
- On 5/8 data the matching setting gives a 10% gain and the best result. A large error remains (1.236, against 0.47 to 0.62 for the other data sets).
- Partial-Fourier corrections cost some edge accuracy. Near-edge error rises from 4.0 to 4.3 (7/8), 4.2 to 4.7 (6/8) and 6.4 to 8.4 (5/8), while full-Fourier data improves from 4.2 to 3.5.
- Over the whole image, edges included, the partial-Fourier corrections do not lower the RMSE. It rises from 2.78 to 2.96 (7/8), 2.96 to 3.24 (6/8) and 4.52 to 5.76 (5/8). Full-Fourier data improves from 2.91 to 2.39, and `--pf 1` on 7/8 data gives 2.64.
- Zeroing the other side of k-space gives the same picture.
- With noise at SNR 30 (five seeds), the matching setting is still best for 7/8, 6/8 and 5/8 data. The gains are 26%, 24% and 10%. The `5/8` setting still increases the error on every other data set.
- The size of the gain depends on the phase. With zero phase (a purely real object), the gains shrink to 8% to 15%.

### Figures

![Phantom images and error maps for each acquired partial-Fourier factor](docs/qc_images.png)

*Left: truth, uncorrected, corrected with the matching factor and corrected as full Fourier. Right: absolute error maps for the same three images, with the far-field and near-edge RMSE in each title.*

![Phase-encode profiles and RMSE heatmaps](docs/qc_profiles.png)

*Phase-encode profiles through the brain (column 58), and heatmaps of far-field RMSE, absolute and relative to uncorrected, for every acquired and correction factor.*

This is one phantom. Real data will differ, so read the table for the ordering, not for the size of the gains.

## Limits

- Magnitude (real-valued) images only. Complex images are not supported.
- Supported partial-Fourier factors are 1, 7/8, 6/8 and 5/8.
- For `5/8`, the original RPG code applies only the phase-encode unringing to magnitude data. This tool runs the pipeline described in the paper, so its `--pf 5/8` output differs from the original there. The other settings agree with it to within 0.02% (relative RMS) on the phantom.
- The partial-Fourier pipelines model zero-filled reconstruction. The phantom tests do not cover data reconstructed with POCS, homodyne or similar methods.
- Each plane needs lines longer than `maxw + 2` samples, including the sub-lines the partial-Fourier scheme creates.

## Licence and citing

**Licence.** The NYU Non-Commercial Research License, copied unchanged from [NYU-DiffusionMRI/rpg-degibbs](https://github.com/NYU-DiffusionMRI/rpg-degibbs) (see `LICENSE`). Non-commercial research only, not for clinical care, no patent rights granted.

**Required citation.** Clause 6 of the licence requires any publication of research that involves this software to cite:

- Lee, Novikov, Fieremans. Removal of partial Fourier-induced Gibbs (RPG) ringing artifacts in MRI. *Magn Reson Med* 86(5):2733-2750, 2021. [doi:10.1002/mrm.28830](https://doi.org/10.1002/mrm.28830). This is the partial-Fourier method implemented here (`--pf 7/8`, `6/8`, `5/8`).

**Also cite.** Please also cite the method RPG extends. It is the full-Fourier part of this tool (`--pf 1`) and the building block of the partial-Fourier pipelines:

- Kellner, Dhital, Kiselev, Reisert. Gibbs-ringing artifact removal based on local subvoxel-shifts. *Magn Reson Med* 76(5):1574-1581, 2016. [doi:10.1002/mrm.26054](https://doi.org/10.1002/mrm.26054).

Both references are in [`CITATION.bib`](CITATION.bib) for BibTeX and BibLaTeX.
