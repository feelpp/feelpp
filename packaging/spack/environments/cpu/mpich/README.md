# CPU MPICH 3.4.3 environment

This environment supports the Feel++ core component image `spack:mpich`.
It uses MPICH 3.4.3 with the CH3 device. The repository's OCI bake pins
Spack `v1.0.0`; the local validation installation is Spack `1.2.0.dev0`
(`7e864787bd15a516314d86002ce1cbabde7cbbe9`). The builtin recipe
accepts `device=ch3 netmod=tcp`, but emits
`--with-device=ch3:nemesis:tcp`. The small `packages/mpich` overlay in the
repository's existing `feelpp` package repo changes just that argument to
`--with-device=ch3` for this spec. No other MPICH configure arguments change.

The manifest sets `mpi` to `mpich@3.4.3 device=ch3 netmod=tcp`, uses
`concretizer.unify: true`, and installs with 15 build jobs and three concurrent
packages. Its view is `$SPACK_USER_CACHE_PATH/views/feelpp/cpu-mpich`.

Local concretization on Ubuntu 24.04/zen2 selected
`mpich@3.4.3 device=ch3 netmod=tcp /s6qwr5qy`, `petsc@3.25.0 /cymaisjq`,
and `slepc@3.25.0 /ndhxt3q3`. Boost, ARPACK, FFTW, mpi4py, and PETSc all
depend directly on that same MPICH hash. SLEPc and slepc4py reach it through
PETSc. The local `spack.lock` is generated state and is excluded from image
contexts; reconcretization inside the OCI image may produce different hashes.

The OCI bake with Spack `v1.0.0` (`73eaea13f381e3495299284856fd02a64e1d154c`)
concretized MPICH `3.4.3 /4zscrufl`, PETSc `3.25.5 /7yfidbqb`, and SLEPc
`3.25.2 /4drpfgt3`. The image lockfile has one MPI provider, and direct MPI
dependencies including Boost, ARPACK, FFTW, HDF5, Hypre, MUMPS, ScaLAPACK,
mpi4py, and PETSc all point to that MPICH node. `mpichversion` in the image
reports `MPICH Device: ch3:nemesis` and the exact `--with-device=ch3`
configure argument. `mpicc -show` links `-lmpi` from that MPICH installation.

```bash
export FEELPP_REPO_ROOT="$PWD"
spack -e packaging/spack/environments/cpu/mpich concretize -f
spack -e packaging/spack/environments/cpu/mpich spec -Il mpich@3.4.3
spack -e packaging/spack/environments/cpu/mpich python -c 'import spack.environment as ev; s=next(x for x in ev.active_environment().all_specs() if x.name=="mpich"); print(s); print([a for a in s.package.configure_args() if a.startswith("--with-device=")])'
```

For LUMI, use its documented `singularity-bindings` host Cray MPICH binding
and launch with `srun`: https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/container-jobs/#using-the-host-mpi.
LUMI's MPICH 3.4.3/CH3 example is at
https://docs.lumi-supercomputer.eu/software/containers/singularity/#building-containers-on-local-hardware.
A container MPICH run only tests the image's own MPI and is not evidence that
the host Cray MPI or Slingshot provider was loaded.

The 2026-09-26 core images are `ghcr.io/feelpp/feelpp-env:spack-mpich` and
`ghcr.io/feelpp/feelpp:spack-mpich`; the core SIF is available as
`oras://ghcr.io/feelpp/feelpp:spack-mpich-sif` (SIF SHA-256
`27ed1b37c822d5c0d6a445e4335300201eed775c05da31a38b900aeb2efa5a35`).
On LUMI, host-MPI probes completed on one, two, and three distinct nodes
(jobs `22362050`, `22362051`, `22362052`). Each rank loaded Cray MPICH
`8.1.29.7`, host libfabric `1.22.0`, and `libcxi.so.1`. The available
`singularity-bindings/24.03` module contains stale libfabric, PALS, ROCm,
and `libjson-c` paths after LUMI's system update; the probe substituted the
current host paths and omitted unavailable binds. This confirms a small MPI
communicator and host library selection, not a Feel++ solver workload.
