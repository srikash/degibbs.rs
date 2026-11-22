# Changelog

## 0.4.1

No change to the program: the output is identical to 0.4.0.

- The "Cite this repository" box on GitHub now cites the software. Lee 2021 and Kellner 2016 are listed as references in `CITATION.cff`.
- The container image is built from the crate published on crates.io, so it always matches `cargo install degibbs`.
- `latest` is the newest release. The `main` and `sha-<commit>` image tags are no longer published.
- 0.4.0 was released on GitHub only; 0.4.1 is the first version on crates.io.

## 0.4.0

The first release with a stability promise: within 0.4.x, option names, defaults, output format and sidecar keys do not change, and the output is bit-identical for the same input, settings and CPU. With `--pigz` the image values are the same but the `.nii.gz` file bytes differ.

### Changes from 0.2 that can break scripts

- **The partial Fourier factor of a magnitude image is now required.** It comes from `--pf`, from `PartialFourier` in the sidecar, or from the new `--pocs`. Before, DeGibbs assumed full Fourier with a notice. If your scan was fully sampled, add `--pf 1`; if the scanner already filled in the missing k-space (POCS), add `--pocs`; if only the original scan has a sidecar, add `--json`.
- **Existing output files are no longer replaced** unless you pass `--overwrite`. DeGibbs never writes over its input.
- **A `.json` file is written beside every output** with the input sidecar's keys plus a `DeGibbs` record of the version and the settings used. In the record, `PartialFourier` is a number, axes are `i`, `j`, `k`, and `Command` is an array of arguments.
- `-h` and `--help` show a short help; the expert options (`--axes`, `--nshifts`, `--minw`, `--maxw`, `--threads`) are in `--full-help`.

### Fixed

- **An output name with extra dots no longer writes over the input.** A name such as `a.fixed.nii` was written as `a.nii` by the NIfTI writer, so `degibbs a.nii a.fixed.nii` replaced the input and reported a file it had not written. Outputs are now written to a temporary name and renamed, so a half-written output is never left behind either.

### New

- `--pocs`: the scanner already filled in the missing k-space; treat the image as full Fourier.
- `--dry-run`: check the inputs and print the resolved settings without correcting anything.
- `--overwrite`.
- `--pigz` (expert).
- Every option is checked before the image is read.
- Plain-language definitions for every option.

### Lighter and tidier

- Two shifted copies share one inverse FFT in the shift search.
- Slices are copied in and out in parallel, and a single 2D slice or a few slices use all threads.
- On a 4D test image (126 MB compressed) peak memory fell from 320 to 257 MB and the correction itself ran about 4% faster.
- The `--phase` route holds two full volumes in memory instead of five.
- `.nii.gz` files use the `zlib-rs` backend. `--pigz` (expert) compresses a `.nii.gz` output with the multi-threaded `pigz` instead: smaller files and quicker writes, with the same image data. The container image includes pigz.
