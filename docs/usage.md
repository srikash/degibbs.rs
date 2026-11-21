# Usage

```
degibbs [options] <input> <output>
```

Input may be 2D, 3D or 4D. Every plane spanned by `--axes` is corrected, and all other axes are looped over. The output is float32 (complex64 for a complex NIfTI input) with the input's header and geometry. A record of the settings used is written beside it (see [Output sidecar](#output-sidecar)).

`degibbs -h` shows the options you normally need. `degibbs --full-help` also shows the expert options. `degibbs --version` prints the version.

**Stability.** Within a minor version (for example every 0.4.x) the option names, defaults, output format and sidecar keys do not change, and the output is bit-identical for the same input, settings and CPU. Pin `0.4` in pipelines.

## Required

| argument | meaning |
|---|---|
| `<input>` | image to correct (`.nii` or `.nii.gz`): a magnitude image, or a complex image |
| `<output>` | where to write the corrected image (`.nii` or `.nii.gz`) |

## Acquisition

Read from the input's [sidecar](#bids-sidecar) when not given.

| option | meaning | default |
|---|---|---|
| `--pe-axis i\|j\|k` (or `0\|1\|2`) | phase-encoding direction: the image axis along which the scanner stepped through k-space line by line (a `+`/`-` sign is accepted and ignored) | `PhaseEncodingDirection` from the sidecar; required otherwise |
| `--pf 1\|7/8\|6/8\|5/8` | partial Fourier factor: the fraction of k-space the scanner acquired along the phase-encoding direction (decimals such as `0.75` work) | `PartialFourier` from the sidecar; **required for magnitude images** if the sidecar has none |
| `--pocs` | the scanner already filled in the missing k-space with POCS reconstruction: treat the image as full Fourier and ignore any `PartialFourier` in the sidecar. Needs no `--pf` | off |

## Files

| option | meaning | default |
|---|---|---|
| `--json FILE` | sidecar to read the acquisition details from; use it when only the original scan has one | the input's name with `.json`, if present |
| `--no-json` | do not read any sidecar | off |
| `--phase FILE` | phase image in radians; makes the input complex, see [Complex images](#complex-images) | none |
| `--phase-out FILE` | where to write the corrected phase (required with `--phase`) | none |
| `--overwrite` | replace output files that already exist | off: DeGibbs stops instead |

## Other

| option | meaning |
|---|---|
| `--dry-run` | check everything and print the settings that would be used, then stop without correcting anything |
| `-h`, `--help` | show the options above |
| `--full-help` | show all options, including the expert ones |
| `-V`, `--version` | show the version |

## Expert options

The defaults suit almost all data. The names follow MRtrix `mrdegibbs`.

| option | meaning | default |
|---|---|---|
| `--axes a,b` | the two image axes that form each slice; the plane must contain the phase-encoding axis | `0,1`, or the two axes other than `SliceEncodingDirection` |
| `--nshifts M` | how finely to search for the best sub-voxel shift: 2M+1 shifts evenly spaced from -1/2 to +1/2 voxel are tried at each voxel. Larger is finer but slower. At most 1000 | `20` (41 shifts) |
| `--minw N` | nearest distance, in voxels, used to judge ringing. For each trial shift, DeGibbs adds up the voxel-to-voxel changes from `--minw` to `--maxw` voxels away on each side and keeps the shift with the least change. May be 0; must not exceed `--maxw` | `1` |
| `--maxw N` | farthest distance, in voxels, used to judge ringing. At most 100 | `3` |
| `--threads N` | number of CPU threads | all CPUs the process may use |

By default DeGibbs uses every CPU the process is allowed to use, so it respects container limits and Slurm allocations; the `RAYON_NUM_THREADS` environment variable also works. The result does not depend on the thread count. It can differ in the last bits between CPUs with different SIMD support.

## Safety

- DeGibbs stops if an output file or output sidecar already exists, unless `--overwrite` is given.
- DeGibbs never writes over its input, with or without `--overwrite`.
- All options are checked before the image is read, so a typo fails at once.

## Output sidecar

Beside every output image DeGibbs writes a JSON file with the same name (`out.nii.gz` gives `out.json`). It holds every key of the input sidecar unchanged (so `PhaseEncodingDirection`, `TotalReadoutTime` and the like stay available to later steps), plus one `DeGibbs` object: the version, the full command, the source file name, the resolved phase-encoding axis, slice plane, partial Fourier factor, `Pocs`, `Nshifts`, `Minw`, `Maxw` and `Threads`. With `--phase-out`, the corrected phase gets its own sidecar the same way. If the input sidecar already has a `DeGibbs` object (DeGibbs run twice), it is replaced.

## BIDS sidecar

Keys the tool reads:

| key | sets |
|---|---|
| `PhaseEncodingDirection` (or `PhaseEncodingAxis`) | `--pe-axis` |
| `PartialFourier` | `--pf` |
| `SliceEncodingDirection` | the default `--axes` |

Explicit options override the sidecar.

> [!IMPORTANT]
> DeGibbs never guesses the partial Fourier factor of a magnitude image. Some DICOM converters leave `PartialFourier` out depending on the scanner, and some write no JSON at all, and a processed copy of a scan often has no sidecar. If the factor is missing, DeGibbs stops and tells you to pass `--pf 1` (full Fourier), `--pf 6/8` and so on (partial Fourier), or `--pocs` (the scanner already filled in the missing k-space). If only the original scan has a sidecar, pass it with `--json`.

> [!NOTE]
> Sidecar axes describe the NIfTI they were written with. They become stale if the image is reoriented afterwards (for example with `fslreorient2std`), and the tool cannot detect this.

## MRtrix3 and FSL

```
mrconvert dwi.mif dwi.nii.gz -json_export dwi.json -export_grad_fsl dwi.bvec dwi.bval
degibbs dwi.nii.gz dwi_degibbs.nii.gz --pf 6/8   # drop --pf if dwi.json has PartialFourier; use --pocs or --pf 1 for full Fourier
mrconvert dwi_degibbs.nii.gz out.mif -fslgrad dwi.bvec dwi.bval
```

## Choosing the partial-Fourier factor

Pass the factor the scan was acquired with. The sidecar's `PartialFourier` sets it automatically, and DeGibbs stops if it cannot be found.

| acquired | `--pf` | expected result |
|---|---|---|
| full Fourier | `1` | the standard Kellner correction |
| 7/8 | `7/8` | removes more ringing than `1` |
| 6/8 | `6/8` | removes more ringing than `1` |
| 5/8 | `5/8` | large gain on the phantom with slowly varying phase, much less with fast phase; the RPG paper reports blurring and leftover ringing, so check the result |
| partial Fourier, then POCS reconstruction on the scanner | `--pocs` (same as `--pf 1`) | treated as full Fourier (see below) |

A mismatched value gives a smaller gain, and `5/8` on any other data is the worst setting. Partial-Fourier settings can make thin bright structures such as the skull ring in the phantom worse. The RPG paper recommends 6/8 over 5/8 when you have a choice at acquisition. The evidence is in the [phantom validation](validation.md).

## Complex images

The tool accepts complex data in two forms, both treated as full-Fourier (`--pf 1`) data:

```
degibbs complex.nii.gz out.nii.gz --pe-axis j                                   # complex64 or complex128 NIfTI
degibbs mag.nii.gz out_mag.nii.gz --phase phase.nii.gz --phase-out out_phase.nii.gz --pe-axis j
```

- Complex input is assumed to be POCS-reconstructed. A POCS reconstruction fills in the missing k-space, so the data again ends symmetrically and ordinary unringing applies.
- The sidecar's `PartialFourier` is ignored for complex input, with a notice. An explicit `--pf` below 1 is an error. `--pocs` is accepted and changes nothing.
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

