# Feel++ CPU/OpenMPI macOS Environment

This environment is the macOS-oriented OpenMPI baseline for local Feel++
development.

It differs from `cpu/openmpi/` in a few deliberate ways:

- OpenMPI uses `fabrics=none`; the Linux `cma,ucx` fabric selection pulls UCX,
  and UCX currently fails on macOS while probing `librt`.
- BLAS/LAPACK use the Homebrew OpenBLAS external at
  `/opt/homebrew/opt/openblas`. This avoids the parallel PETSc/MUMPS solve
  crash observed with the `veclibfort`/Accelerate stack, while also avoiding
  the Spack OpenBLAS test-link failure against `libgfortran` with AppleClang.
  Homebrew OpenBLAS is OpenMP-enabled; use `OMP_NUM_THREADS=1` for MPI test
  runs unless you intentionally want nested threading.
- PETSc is built with both MUMPS and SuperLU_DIST. MUMPS remains the default
  direct solver in Feel++. On macOS, Feel++ defaults MUMPS `ICNTL(20)` to `0`
  so MUMPS uses centralized dense right-hand sides; this avoids the Apple
  Silicon parallel distributed-RHS crash at the cost of more root memory for
  large or many RHS solves. SuperLU_DIST provides a second distributed LU
  implementation for macOS diagnostics and fallback runs:
  `--pc-factor-mat-solver-package-type=superlu-dist`.
- Gmsh remains `+opencascade`; the repository overlay passes Tcl's library
  directory explicitly so OpenCASCADE can configure on Darwin.
- Gmsh is built `~cgns~med`; Feel++ does not use CGNS here, and disabling it
  avoids a macOS/Homebrew CGNS header collision during the Gmsh build.
- The standalone Spack `py-gmsh` package is not included because the recipe is
  pinned to a Linux wheel. The `gmsh` package itself installs `gmsh.py`.

Activate and install it through the Feel++ helper:

```sh
source .venv/bin/activate
fpp-spack init
fpp-spack install openmpi-macosx -- -j4 -p2 --fail-fast --show-log-on-error
```

For an interactive shell, activate the same environment with:

```sh
eval "$(fpp-spack activate openmpi-macosx)"
```

After activation, plain Spack commands also use the repository environment:

```sh
spack spec
spack install -j4 -p2 --fail-fast --show-log-on-error
```
