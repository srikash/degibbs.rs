# Containers (Docker and Apptainer)

The image is published to the GitHub Container Registry as `ghcr.io/srikash/degibbs.rs`. It is about 18 MB and runs on Intel/AMD and ARM (including Apple silicon) machines.

- `latest` and `main` follow the newest commit.
- `0.2.1` is the current release. Use a version tag in analysis scripts so the software stays fixed.
- `sha-<commit>` pins an exact commit.

The container's entry point is `degibbs`, so everything after the image name is a normal `degibbs` command. Its working directory is `/data`.

## Docker (laptop or workstation)

```
docker run --rm ghcr.io/srikash/degibbs.rs:0.2.1 --version

docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data ghcr.io/srikash/degibbs.rs:0.2.1 \
    dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

- `-v "$PWD":/data` makes the current folder visible inside the container, so file names are relative to it.
- `-u` keeps the output files owned by you rather than root.
- Podman works with the same image. Add `--userns=keep-id` instead of `-u`, and `:Z` after `/data` on SELinux systems such as Fedora and Rocky.

## Apptainer (HPC)

Pull the image once into a single `.sif` file, on a login or transfer node:

```
apptainer pull degibbs_0.2.1.sif docker://ghcr.io/srikash/degibbs.rs:0.2.1
```

Then run it like any other program. Apptainer shares your user ID, so output files are yours.

```
apptainer run degibbs_0.2.1.sif dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

- Apptainer makes your home and current folder visible by default. Add `--bind /path/to/study:/data` for a project or scratch filesystem, and give absolute paths inside the container.
- A Slurm job, with the thread count matched to the CPUs you asked for:

```
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

apptainer run --bind /scratch/$USER/study:/data degibbs_0.2.1.sif \
    /data/sub-01/dwi.nii.gz /data/sub-01/dwi_degibbs.nii.gz --threads "$SLURM_CPUS_PER_TASK"
```

- Older clusters call the program `singularity`. The commands are the same.

[Back to the README](../README.md)

