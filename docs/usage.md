# Usage

```
degibbs <input> <output> [options]
```

Input may be 2D, 3D or 4D. Every plane spanned by `--axes` is corrected, and all other axes are looped over. The output is float32 with the input's header and geometry.

## Options

| flag | meaning | default |
|---|---|---|
| `--pe-axis i\|j\|k` (or `0\|1\|2`) | phase-encoding axis in the file's voxel axes (a `+`/`-` sign is accepted and ignored) | from the sidecar; required otherwise |
| `--pf 1\|7/8\|6/8\|5/8` | partial-Fourier factor (decimals such as `0.75` work) | from the sidecar, else `1` |
| `--axes 0,1` | in-plane axis pair containing the PE axis | `0,1`, or from `SliceEncodingDirection` |
| `--json FILE` | BIDS sidecar to read | `<input>.json` next to the input |
| `--no-json` | do not look for a sidecar | off |
| `--phase FILE` | phase image in radians; makes the input complex, see [Complex images](#complex-images) | none |
| `--phase-out FILE` | where to write the corrected phase (required with `--phase`) | none |
| `--nshifts N` | subvoxel shifts per voxel | `20` |
| `--minw N`, `--maxw N` | start and end of the oscillation window | `1`, `3` |
| `--threads N` | worker threads | all cores |

## BIDS sidecar

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

## MRtrix3 and FSL

```
mrconvert dwi.mif dwi.nii.gz -json_export dwi.json -export_grad_fsl dwi.bvec dwi.bval
degibbs dwi.nii.gz dwi_degibbs.nii.gz
mrconvert dwi_degibbs.nii.gz out.mif -fslgrad dwi.bvec dwi.bval
```

## Containers (Docker and Apptainer)

The image is published to the GitHub Container Registry as `ghcr.io/srikash/degibbs.rs`. It is about 18 MB and runs on Intel/AMD and ARM (including Apple silicon) machines.

- `latest` and `main` follow the newest commit.
- `0.1.6` is the current release. Use a version tag in analysis scripts so the software stays fixed.
- `sha-<commit>` pins an exact commit.

The container's entry point is `degibbs`, so everything after the image name is a normal `degibbs` command. Its working directory is `/data`.

## Docker (laptop or workstation)

```
docker run --rm ghcr.io/srikash/degibbs.rs:0.1.6 --version

docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data ghcr.io/srikash/degibbs.rs:0.1.6 \
    dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

- `-v "$PWD":/data` makes the current folder visible inside the container, so file names are relative to it.
- `-u` keeps the output files owned by you rather than root.
- Podman works with the same image. Add `--userns=keep-id` instead of `-u`, and `:Z` after `/data` on SELinux systems such as Fedora and Rocky.

## Apptainer (HPC)

Pull the image once into a single `.sif` file, on a login or transfer node:

```
apptainer pull degibbs_0.1.6.sif docker://ghcr.io/srikash/degibbs.rs:0.1.6
```

Then run it like any other program. Apptainer shares your user ID, so output files are yours.

```
apptainer run degibbs_0.1.6.sif dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

- Apptainer makes your home and current folder visible by default. Add `--bind /path/to/study:/data` for a project or scratch filesystem, and give absolute paths inside the container.
- A Slurm job, with the thread count matched to the CPUs you asked for:

```
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

apptainer run --bind /scratch/$USER/study:/data degibbs_0.1.6.sif \
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
| 5/8 | `5/8` | large gain on the phantom with slowly varying phase, much less with fast phase; the RPG paper reports blurring and leftover ringing, so check the result |
| partial Fourier, then POCS reconstruction | `1` | treated as full Fourier (see below) |

A mismatched value gives a smaller gain, and `5/8` on any other data is the worst setting. Partial-Fourier settings can make thin bright structures such as the skull ring in the phantom worse. The RPG paper recommends 6/8 over 5/8 when you have a choice at acquisition. The evidence is in the [phantom validation](validation.md).

## Complex images

The tool accepts complex data in two forms, both treated as full-Fourier (`--pf 1`) data:

```
degibbs complex.nii.gz out.nii.gz --pe-axis j                                   # complex64 or complex128 NIfTI
degibbs mag.nii.gz out_mag.nii.gz --phase phase.nii.gz --phase-out out_phase.nii.gz --pe-axis j
```

- Complex input is assumed to be POCS-reconstructed. A POCS reconstruction fills in the missing k-space, so the data again ends symmetrically and ordinary unringing applies.
- The sidecar's `PartialFourier` is ignored for complex input, with a notice. An explicit `--pf` below 1 is an error.
- On the full-Fourier phantom the complex path matches the magnitude path inside the brain and removes more error in the background and around the skull ring; see the [validation](validation.md#complex-input).
- The phase image must be wrapped phase in radians (-pi to pi, or 0 to 2 pi). Siemens phase stored as integers from -4096 to 4095 becomes radians when multiplied by pi/4096.
- If the data were zero-filled instead of reconstructed, pass the magnitude image and the matching `--pf`. Complex zero-filled data processed as `--pf 1` keeps its wider partial-Fourier ringing.
- POCS needs a good phase estimate. On the phantom, POCS followed by the complex path beats zero-filled magnitude with the matching `--pf` when the phase varies slowly, and loses in the brain when it varies fast; see the [comparison](validation.md#pocs-and-complex-input-for-partial-fourier-data).

## Limits

- Supported partial-Fourier factors are 1, 7/8, 6/8 and 5/8. Complex input is supported for full-Fourier data only.
- The partial-Fourier settings model zero-filled reconstruction.
- A scanner k-space filter, or a smooth join between measured and reconstructed data (homodyne and some POCS versions), reshapes the ringing, so the correction is only approximate.
- Each plane needs lines longer than `maxw + 2` samples, including the sub-lines the partial-Fourier scheme creates.
- `5/8` magnitude input follows the RPG paper. The DESIGNER-v2 C++ module differs there (see the [validation notes](validation.md#agreement-with-the-designer-v2-rpg-module)).

[Back to the README](../README.md)

