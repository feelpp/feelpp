# Docker Image Metadata

`packaging/docker/` is the repository-owned home for source-build OCI image
metadata used by `fpp-pkg image ...`.

Current scope:

- Ubuntu and Debian environment metadata under `config/`
- shared environment assets under `assets/`
- base environment templates plus component/full Dockerfile templates under
  `templates/`

Design rules:

- keep only source-of-truth inputs here
- do not commit generated Docker contexts here
- do not mirror the full legacy external Docker repository here
- keep Ubuntu and Debian environment images on the `apt` path
- keep Spack image generation under `packaging/spack/` plus `ops/`

Current supported variant:

- `feelpp-env`
  - base development image used by CI to build Feel++ from Ubuntu, Debian, and
    Spack environments

Current planned environment group:

- `env-current`
  - `ubuntu:noble` / `ubuntu-24.04`
  - `ubuntu:resolute` / `ubuntu-26.04`
  - `debian:trixie` / `debian-13`
  - `spack:openmpi` / `spack-openmpi`

Later planned rows:

- Fedora after a `dnf`/RPM backend exists
- `debian:testing`
- `debian:sid`

Generated bake contexts belong under the packaging job root and are produced by
`fpp-pkg image bake --target ...`.

## preCICE dependency environments

Noble, Trixie, and Resolute environment recipes install preCICE **3.4.1**.
Distribution metadata pins the official amd64 package URL and SHA-256. The arm64
path builds the same verified source release against the distribution's MPI,
Boost, Eigen, PETSc and Python dependencies. Arm64 runtime validation remains a
publication gate; generating a Dockerfile does not verify that build.

The generated context copies `assets/install-precice.sh` only for these targets.
The installer checks package identity/version/architecture, verifies checksums,
runs `ldconfig` and `precice-version`, and records its inputs in
`/usr/local/share/feelpp/image/precice.tsv`. Rebuild environment images before
building the corresponding Feel++ core images with preCICE enabled.
