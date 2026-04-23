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
- `openmpi4/`
  - legacy imported bootstrap manifest from the old hidden `.spack/` path
  - visible for comparison and migration
  - not the preferred baseline for new shared support
