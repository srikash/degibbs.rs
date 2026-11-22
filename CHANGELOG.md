# Changelog

## 0.4.0

The first release with a stability promise: within 0.4.x, option names, defaults, output format and sidecar keys do not change, and the output is bit-identical for the same input, settings and CPU. This is about the image values: `.nii.gz` files compressed with `pigz` differ in file bytes (`--pigz off` removes the difference).

### Changes from 0.2 that can break scripts

- **The partial Fourier factor of a magnitude image is now required.** It comes from `--pf`, from `PartialFourier` in the sidecar, or from the new `--pocs`. Before, DeGibbs assumed full Fourier with a notice. If your scan was fully sampled, add `--pf 1`; if the scanner already filled in the missing k-space (POCS), add `--pocs`; if only the original scan has a sidecar, add `--json`.
- **Existing output files are no longer replaced** unless you pass `--overwrite`. DeGibbs never writes over its input.
- **A `.json` file is written beside every output** with the input sidecar's keys plus a `DeGibbs` record of the version and the settings used. In the record, `PartialFourier` is a number, axes are `i`, `j`, `k`, and `Command` is an array of arguments.
- `-h` and `--help` show a short help; the expert options (`--axes`, `--nshifts`, `--minw`, `--maxw`, `--threads`) are in `--full-help`.

### New

- `--pocs`: the scanner already filled in the missing k-space; treat the image as full Fourier.
- `--dry-run`: check the inputs and print the resolved settings without correcting anything.
- `--overwrite`.
- `--pigz auto|on|off` (expert).
- Every option is checked before the image is read.
- Plain-language definitions for every option.

### Faster and lighter

- Two shifted copies share one inverse FFT in the shift search.
- Slices are copied in and out in parallel, and a single 2D slice or a few slices use all threads.
- The `--phase` route holds two full volumes in memory instead of five.
- `.nii.gz` files use the `zlib-rs` backend, and a `.nii.gz` output is compressed with the multi-threaded `pigz` when it is installed (always in the container image); see `--pigz`. Pigz output is smaller and quicker to write; the image data are the same, but the compressed bytes differ from the built-in compressor's. Use `--pigz off` for byte-for-byte reproducible `.nii.gz` files.
