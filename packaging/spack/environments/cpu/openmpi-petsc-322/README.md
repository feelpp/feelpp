# Feel++ Shared CPU/OpenMPI PETSc 3.22 Environment

This environment mirrors `cpu/openmpi/` but pins the solver stack to the
PETSc/SLEPc 3.22 line. Its purpose is to validate Feel++ against the solver
generation used by Debian Trixie-like deployments while keeping the same
OpenMPI/UCX CPU workflow as the main Spack environment.

Pinned solver stack:

- `petsc@3.22.5`
- `py-petsc4py@3.22.5`
- `slepc@3.22.2`
- `py-slepc4py@3.22.2`

Typical use:

```bash
spacktivate ./packaging/spack/environments/cpu/openmpi-petsc-322
spack concretize -f
spack install
cmake --preset release-clang-spack
cmake --build build/release-clang-spack -j
ctest --preset release-clang-spack -R toolbox_heat-testsuite_test2d_1-np-1
```

Use the same `SPACK_USER_CONFIG_PATH` and `SPACK_USER_CACHE_PATH` conventions
as `cpu/openmpi/` so the generated view and repo cache are not mixed with an
unrelated user-local Spack setup.
