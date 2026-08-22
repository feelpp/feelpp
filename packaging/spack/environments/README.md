# Shared Spack Environments

This directory contains shared, repository-owned Spack environments.

These environments are intended to describe supported Feel++ development and CI
stacks. They are not the place for:

- one developer's local compiler externals
- private mirrors or credentials
- user-local cache state

Environment status:

- `cpu/openmpi/`
  - supported shared baseline
  - focuses on the immediate CPU/OpenMPI PETSc/HPDDM workflow
- `cpu/openmpi5/`
  - portable OpenMPI 5 validation environment for Gaya and LUMI-C
  - enables UCX and OFI while keeping CXI/site externals out of the manifest
- `cpu/openmpi-petsc-322/`
  - supported CPU/OpenMPI validation environment
  - mirrors `cpu/openmpi/` while pinning the PETSc/SLEPc stack to 3.22
- `cpu/openmpi-macosx/`
  - supported macOS CPU/OpenMPI environment
  - avoids UCX/CMA and uses OpenBLAS for BLAS/LAPACK
- `openmpi4/`
  - legacy imported bootstrap manifest from the old hidden `.spack/` path
  - visible for comparison and migration
  - not the preferred baseline for new shared support
