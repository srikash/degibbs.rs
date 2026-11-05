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

Or build from source (needs Rust 1.85 or newer; this puts `degibbs` on your `PATH`):

```
cargo install --git https://github.com/srikash/degibbs.rs
```

## Run

```
degibbs in.nii.gz out.nii.gz --pe-axis j --pf 6/8
degibbs dwi.nii.gz dwi_degibbs.nii.gz            # reads dwi.json (BIDS sidecar)
```

> [!WARNING]
> If no partial-Fourier factor is available, the tool assumes full Fourier and says so on stderr. This is the most likely way to get a wrong result. Check that message for partial-Fourier data.

## Documentation

- [Usage](docs/usage.md): options, BIDS sidecar, MRtrix3 and FSL, choosing the partial-Fourier factor, complex images, limits.
- [Containers](docs/containers.md): Docker, Podman and Apptainer, including a Slurm job.
- [Validation](docs/validation.md): phantom results and figures.

## Licence and citing

**Licence.** The NYU Non-Commercial Research License, copied unchanged from [NYU-DiffusionMRI/rpg-degibbs](https://github.com/NYU-DiffusionMRI/rpg-degibbs) (see `LICENCE`). Non-commercial research only, not for clinical care, no patent rights granted.

**Required citation.** Clause 6 of the licence requires any publication of research that involves this software to cite:

- Lee, Novikov, Fieremans. Removal of partial Fourier-induced Gibbs (RPG) ringing artifacts in MRI. *Magn Reson Med* 86(5):2733-2750, 2021. [doi:10.1002/mrm.28830](https://doi.org/10.1002/mrm.28830). This is the partial-Fourier method implemented here (`--pf 7/8`, `6/8`, `5/8`).

**Also cite.** Please also cite the method RPG extends. It is the full-Fourier part of this tool (`--pf 1`) and the building block of the partial-Fourier pipelines:

- Kellner, Dhital, Kiselev, Reisert. Gibbs-ringing artifact removal based on local subvoxel-shifts. *Magn Reson Med* 76(5):1574-1581, 2016. [doi:10.1002/mrm.26054](https://doi.org/10.1002/mrm.26054).

Both references are in [`CITATION.bib`](CITATION.bib) for BibTeX and BibLaTeX.
