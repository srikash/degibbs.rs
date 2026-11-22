# degibbs.rs

Gibbs-ringing removal for MRI in pure Rust (assisted in no small part by Claude 5.5).

- **Method:** local subvoxel shifts (Kellner et al., 2016), with the RPG extension for partial-Fourier acquisitions (Lee et al., 2021).
- **Input and output:** NIfTI in, NIfTI out. Works with MRtrix3, FSL, ANTs or any other pipeline.
- **Partial Fourier:** full Fourier, 7/8, 6/8 and 5/8 for magnitude images; complex images as full Fourier (after POCS).
- **Validation:** [phantom results and figures](docs/validation.md).
- **Binary:** one self-contained executable, no other dependencies. Also available as a Docker, Podman and Apptainer image.
- **Licence:** non-commercial research only, not for clinical care.

## Install

Download the static Linux binary (x86-64) from the [v0.2.1 release](https://github.com/srikash/degibbs.rs/releases/tag/v0.2.1). It runs on Rocky 9 and 10, Ubuntu 22.04 and 24.04, and Fedora 44.

```
tar xzf degibbs-v0.2.1-x86_64-linux-musl.tar.gz
```

Or build from source (needs Rust 1.88 or newer; this puts `degibbs` on your `PATH`):

```
cargo install --git https://github.com/srikash/degibbs.rs
```

## Run

```
degibbs in.nii.gz out.nii.gz --pe-axis j --pf 6/8
degibbs dwi.nii.gz dwi_degibbs.nii.gz            # reads dwi.json (BIDS sidecar)
degibbs in.nii.gz out.nii.gz --pe-axis j --pocs  # the scanner already filled in the missing k-space
degibbs -h                                       # the options you normally need
degibbs --full-help                              # everything, including expert options
```

> [!IMPORTANT]
> DeGibbs never guesses the partial Fourier factor. If the sidecar has no `PartialFourier` (some converters leave it out, and a processed copy often has no sidecar), pass `--pf` (`1` for full Fourier), `--pocs`, or `--json` pointing at the original scan's sidecar. Without one of these DeGibbs stops with an error.

Existing outputs are never replaced unless you pass `--overwrite`. Every output gets a `.json` beside it with the settings used.

**Stability.** Within a minor version (for example 0.4.x), option names, defaults, output format and sidecar keys do not change, and the output is bit-identical for the same input, settings and CPU. Pin `0.4` in pipelines.

## Documentation

- [Usage](docs/usage.md): options, BIDS sidecar, MRtrix3 and FSL, choosing the partial-Fourier factor, complex images, limits.
- [Containers](docs/containers.md): Docker, Podman and Apptainer, including a Slurm job.
- [Validation](docs/validation.md): phantom results and figures.
- [Development](docs/development.md): tests, checking a change against an earlier build, the phantom validation and the benchmark.

## Licence and citing

**Licence.** The NYU Non-Commercial Research License, copied unchanged from [NYU-DiffusionMRI/rpg-degibbs](https://github.com/NYU-DiffusionMRI/rpg-degibbs) (see `LICENCE`). Non-commercial research only, not for clinical care, no patent rights granted.

**Required citation.** Clause 6 of the licence requires any publication of research that involves this software to cite:

- Lee, Novikov, Fieremans. Removal of partial Fourier-induced Gibbs (RPG) ringing artifacts in MRI. *Magn Reson Med* 86(5):2733-2750, 2021. [doi:10.1002/mrm.28830](https://doi.org/10.1002/mrm.28830). This is the partial-Fourier method implemented here (`--pf 7/8`, `6/8`, `5/8`).

**Also cite.** Please also cite the method RPG extends. It is the full-Fourier part of this tool (`--pf 1`) and the building block of the partial-Fourier pipelines:

- Kellner, Dhital, Kiselev, Reisert. Gibbs-ringing artifact removal based on local subvoxel-shifts. *Magn Reson Med* 76(5):1574-1581, 2016. [doi:10.1002/mrm.26054](https://doi.org/10.1002/mrm.26054).

Both references are in [`CITATION.bib`](CITATION.bib) for BibTeX and BibLaTeX, and [`CITATION.cff`](CITATION.cff) drives the "Cite this repository" button on GitHub.
