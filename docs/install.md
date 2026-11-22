# Install

Pick one. Each ends with a check that prints the version.

## 1. Native

**Download the binary** (Linux, x86-64, static; runs on Rocky 9 and 10, Ubuntu 22.04 and 24.04, Fedora 44):

```
curl -LO https://github.com/srikash/degibbs.rs/releases/download/v0.4.1/degibbs-v0.4.1-x86_64-linux-musl.tar.gz
tar xzf degibbs-v0.4.1-x86_64-linux-musl.tar.gz
./degibbs --version
```

Move `degibbs` to a folder on your `PATH`, for example `~/.local/bin`.

**Or install with Cargo** (any platform; needs Rust 1.88 or newer):

```
cargo install degibbs
degibbs --version
```

The optional `--pigz` option needs the [pigz](https://zlib.net/pigz/) program on your `PATH`.

## 2. Docker

```
docker run --rm ghcr.io/srikash/degibbs.rs:0.4.1 --version
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data ghcr.io/srikash/degibbs.rs:0.4.1 \
    dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

`-v "$PWD":/data` shows the current folder inside the container, and `-u` keeps the outputs owned by you. The image includes pigz.

## 3. Apptainer (also called Singularity)

No `.sif` file is published. Pull the image once, on a login or transfer node, and Apptainer builds it for you:

```
apptainer pull degibbs_0.4.1.sif docker://ghcr.io/srikash/degibbs.rs:0.4.1
apptainer run degibbs_0.4.1.sif --version
apptainer run degibbs_0.4.1.sif dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

Add `--bind /path/to/study:/data` for a project or scratch filesystem. A Slurm job example is in [Containers](containers.md).

## 4. Podman

The same image works. Use `--userns=keep-id` in place of `-u`, and `:Z` after `/data` on SELinux systems such as Fedora and Rocky:

```
podman run --rm ghcr.io/srikash/degibbs.rs:0.4.1 --version
podman run --rm --userns=keep-id -v "$PWD":/data:Z ghcr.io/srikash/degibbs.rs:0.4.1 \
    dwi.nii.gz dwi_degibbs.nii.gz --pe-axis j --pf 6/8
```

Use a version tag in analysis scripts so the software stays fixed. [Usage](usage.md) lists every option.

[Back to the README](../README.md)
