# degibbs.rs

Gibbs-ringing removal for magnitude MRI, in pure Rust.

- **Method:** local subvoxel shifts (Kellner et al., 2016), with the RPG extension for partial-Fourier acquisitions (Lee et al., 2021).
- **Input and output:** NIfTI in, NIfTI out. Works with MRtrix3, FSL, ANTs or any other pipeline.
- **Partial Fourier:** full Fourier, 7/8, 6/8 and 5/8.
- **Binary:** one self-contained executable, no other dependencies. Also available as a Docker, Podman and Apptainer image.
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

Download the static Linux binary (x86-64) from the [v0.1.5 release](https://github.com/srikash/degibbs.rs/releases/tag/v0.1.5). It runs on Rocky 9 and 10, Ubuntu 22.04 and 24.04, and Fedora 44.

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

The image is published to the GitHub Container Registry as `ghcr.io/srikash/degibbs.rs`. It is about 18 MB and runs on Intel/AMD and ARM (including Apple silicon) machines.

- `latest` and `main` follow the newest commit.
- `0.1.5` is the current release. Use a version tag in analysis scripts so the software stays fixed.
- `sha-<commit>` pins an exact commit.

The container's entry point is `degibbs`, so everything after the image name is a normal `degibbs` command. Its working directory is `/data`.

### Docker (laptop or workstation)

```
docker run --rm ghcr.io/srikash/degibbs.rs:0.1.5 --version

docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data ghcr.io/srikash/degibbs.rs:0.1.5 \
    dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

- `-v "$PWD":/data` makes the current folder visible inside the container, so file names are relative to it.
- `-u` keeps the output files owned by you rather than root.
- Podman works with the same image. Add `--userns=keep-id` instead of `-u`, and `:Z` after `/data` on SELinux systems such as Fedora and Rocky.

### Apptainer (HPC)

Pull the image once into a single `.sif` file, on a login or transfer node:

```
apptainer pull degibbs_0.1.5.sif docker://ghcr.io/srikash/degibbs.rs:0.1.5
```

Then run it like any other program. Apptainer shares your user ID, so output files are yours.

```
apptainer run degibbs_0.1.5.sif dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

- Apptainer makes your home and current folder visible by default. Add `--bind /path/to/study:/data` for a project or scratch filesystem, and give absolute paths inside the container.
- A Slurm job, with the thread count matched to the CPUs you asked for:

```
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

apptainer run --bind /scratch/$USER/study:/data degibbs_0.1.5.sif \
    /data/sub-01/dwi.nii.gz /data/sub-01/dwi_degibbs.nii.gz --threads "$SLURM_CPUS_PER_TASK"
```

- Older clusters call the program `singularity`. The commands are the same.

## Choosing the partial-Fourier factor

Pass the factor the scan was acquired with. The sidecar's `PartialFourier` sets it automatically.

| acquired | `--pf` | expected result |
|---|---|---|
| full Fourier | `1` | the standard Kellner correction |
| 7/8 | `7/8` | removes more ringing than `1` |
| 6/8 | `6/8` | removes more ringing than `1` |
| 5/8 | `5/8` | large gain on the phantom when the phase varies slowly (17% when it varies fast); the RPG paper reports blurring and leftover ringing, so check the result |

A mismatched value gives a smaller gain, and `5/8` on any other data is the worst setting. Thin bright structures can get worse with the partial-Fourier settings (see the skull ring below). The RPG paper recommends 6/8 over 5/8 when you have a choice at acquisition.

## Validation on a Gibbs phantom

The script is [`validation/phantom.py`](validation/phantom.py). It regenerates every table and figure below:

```
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --figures docs   # tables and figures
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --phase none     # zero phase
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --phase fast     # phase turns over in 10 voxels
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --side lo        # other end of k-space unsampled
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --snr 20 --offsets 4 --seeds 3
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --snr 8 --offsets 4 --seeds 3
```

### Setup

- The phantom is the modified Shepp-Logan head (Toft, 1996), scaled so the skull is 100 and the brain 20, with a slowly varying phase map as in real data.
- A 96 × 96 slice is sampled from the object at 8× the resolution, so k-space truncation produces true Gibbs ringing.
- For partial Fourier, the missing lines at one end of k-space are zeroed.
- The magnitude is then taken, as for zero-filled scanner data.
- The truth is the high-resolution object averaged into voxels.
- Errors are RMSE against that truth, averaged over 8 sub-voxel offsets of the object.
- Far-field error is measured more than 4 voxels from any edge, so it shows ringing rather than blur. It is reported for brain tissue (550 pixels, inside the skull with ventricles excluded) and for the background outside the head (3096 pixels).
- Near-edge error is measured within 2 voxels of an edge inside the brain (1952 pixels). The skull ring and its 2-voxel margin (1532 pixels) are reported separately.
- Noise is quoted as brain intensity divided by the standard deviation of each complex channel. It is added on the sampled k-space lines only, so zero-filled partial-Fourier images carry less noise than full-Fourier ones.

### Results

Far-field RMSE in the brain (lower is better; bold marks the best correction in each row):

| acquired | none | `--pf 1` | `--pf 7/8` | `--pf 6/8` | `--pf 5/8` |
|---|---|---|---|---|---|
| full Fourier | 1.266 | **0.182** | 0.284 | 0.187 | 0.657 |
| 7/8 | 1.229 | 0.559 | **0.251** | 0.576 | 0.714 |
| 6/8 | 1.354 | 0.875 | 0.783 | **0.471** | 0.976 |
| 5/8 | 1.840 | 1.631 | 1.651 | 1.557 | **0.780** |

Far-field RMSE in the background (truth is zero):

| acquired | none | `--pf 1` | `--pf 7/8` | `--pf 6/8` | `--pf 5/8` |
|---|---|---|---|---|---|
| full Fourier | 0.731 | 0.717 | 0.712 | **0.706** | 0.955 |
| 7/8 | 0.753 | 0.740 | 0.737 | **0.731** | 0.936 |
| 6/8 | 0.812 | 0.801 | 0.798 | **0.790** | 0.946 |
| 5/8 | 1.338 | 1.310 | 1.314 | **1.306** | 1.360 |

### Findings

- In the brain, the matching setting removes 86% (full Fourier), 80% (7/8), 65% (6/8) and 58% (5/8) of the far-field error. It is the best setting in every row with this slowly varying phase.
- On 7/8 and 6/8 data the matching setting has 55% and 46% less error than `--pf 1`. This holds at every offset, between 50% and 61% less for 7/8 and between 45% and 47% less for 6/8.
- On full-Fourier data `--pf 1` beats `--pf 6/8` by 2% to 4% at every offset.
- The `5/8` setting is the worst correction on every other data set, but it still beats no correction inside the brain.
- Near edges inside the brain, the matching setting also lowers the error: 1.37 to 0.92 (full Fourier), 1.27 to 1.05 (7/8), 1.49 to 1.24 (6/8) and 2.24 to 2.03 (5/8).
- In the background, no setting removes more than 4% of the error, because the magnitude of the ringing there is not an oscillation about a plateau.
- Around the skull ring, a thin structure at five times the brain intensity, the partial-Fourier corrections raise the error: 6.06 to 6.95 (7/8), 6.46 to 7.53 (6/8) and 9.98 to 13.18 (5/8). Full-Fourier data improves there, from 6.16 to 5.25. Other thin, bright structures may behave the same way. The whole-image RMSE rises accordingly, for example from 2.76 to 3.01 for 7/8 data.
- Zero phase, or zeroing the other end of k-space, gives the same ordering and gains within 5 percentage points.
- With a fast phase map that turns over in about 10 voxels, the uncorrected error rises sharply for 6/8 and 5/8 data (1.77 and 4.00, against 1.35 and 1.84). The matching setting removes 33% (full Fourier), 36% (7/8), 23% (6/8) and 17% (5/8), and it is no longer always the best: `--pf 7/8` beats `--pf 1` on full-Fourier data (0.62 against 0.69), `--pf 7/8` edges out `--pf 6/8` on 6/8 data (1.34 against 1.36), and `--pf 6/8` beats `--pf 5/8` on 5/8 data (3.23 against 3.34).
- At a noise level of 20, the matching setting removes 51% (full Fourier), 51% (7/8), 48% (6/8) and 52% (5/8) of the error. At 8 the gains are 29%, 34%, 32% and 41%. The truth has no noise, so part of this drop is noise smoothing rather than ringing removal. For example, `--pf 6/8` beats `--pf 1` on full-Fourier data at both noise levels.
- The RPG paper reports robust removal for 7/8 and 6/8 on its own 90 × 90 Shepp-Logan phantom with added phase, and blurring with leftover ringing for 5/8. It attributes the 5/8 result to phase varying within the long ringing interval. The slow phase used here does not test that. The fast phase map does, and there the 5/8 gain falls to 17%, in line with the paper's warning.

### Figures

![Phantom images and error maps for each acquired partial-Fourier factor](docs/qc_images.png)

*One offset, no noise. Left: truth, uncorrected, corrected with the matching factor and corrected as full Fourier. Right: absolute error maps for the same three images, with the far-field RMSE in the brain and in the background in each title.*

![Phase-encode profiles and RMSE heatmaps](docs/qc_profiles.png)

*Phase-encode profiles through the brain (column 58, one offset), and far-field RMSE heatmaps for the brain and the background (average over 8 offsets).*

This is one phantom. Real data will differ, and the ordering depends on how fast the phase varies. Read the tables for the ordering, not for the size of the gains.

## Limits

- Magnitude (real-valued) images only. Complex images are not supported.
- Supported partial-Fourier factors are 1, 7/8, 6/8 and 5/8.
- For `5/8` magnitude input, the C++ module in DESIGNER-v2 discards the sub-image step and only unrings along the phase-encode axis. The NYU MATLAB code and this tool run the pipeline described in the paper. For the other settings, this tool agrees with the DESIGNER-v2 module to within 0.002% (relative RMS) on the phantom.
- The partial-Fourier pipelines model zero-filled reconstruction. The phantom tests do not cover data reconstructed with POCS, homodyne or similar methods.
- Each plane needs lines longer than `maxw + 2` samples, including the sub-lines the partial-Fourier scheme creates.

## Licence and citing

**Licence.** The NYU Non-Commercial Research License, copied unchanged from [NYU-DiffusionMRI/rpg-degibbs](https://github.com/NYU-DiffusionMRI/rpg-degibbs) (see `LICENSE`). Non-commercial research only, not for clinical care, no patent rights granted.

**Required citation.** Clause 6 of the licence requires any publication of research that involves this software to cite:

- Lee, Novikov, Fieremans. Removal of partial Fourier-induced Gibbs (RPG) ringing artifacts in MRI. *Magn Reson Med* 86(5):2733-2750, 2021. [doi:10.1002/mrm.28830](https://doi.org/10.1002/mrm.28830). This is the partial-Fourier method implemented here (`--pf 7/8`, `6/8`, `5/8`).

**Also cite.** Please also cite the method RPG extends. It is the full-Fourier part of this tool (`--pf 1`) and the building block of the partial-Fourier pipelines:

- Kellner, Dhital, Kiselev, Reisert. Gibbs-ringing artifact removal based on local subvoxel-shifts. *Magn Reson Med* 76(5):1574-1581, 2016. [doi:10.1002/mrm.26054](https://doi.org/10.1002/mrm.26054).

Both references are in [`CITATION.bib`](CITATION.bib) for BibTeX and BibLaTeX.
