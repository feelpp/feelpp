# Feel++ Shared CPU/OpenMPI Environment

This is the first supported shared Spack environment for Feel++.

Design goals:

- provide a visible repository-owned baseline
- focus on the immediate CPU/OpenMPI workflow
- encode PETSc with explicit `+hpddm` support
- pin the shared solver stack to the PETSc/SLEPc `3.25.0` line
- stay smaller than the imported `openmpi4/` bootstrap manifest
- avoid embedding one developer's local compiler and mirror assumptions

Scope:

- CPU only
- OpenMPI-based MPI stack
- PETSc/SLEPc solver stack
- PETSc/SLEPc HPDDM runtime stack
- Python support needed by Feel++ bindings and validation
- Python checker/runtime support for SymPy-backed quickstart scripts
- core packages needed for the current Feel++-with-Spack development workflow

Deliberately not in scope yet:

- CUDA or ROCm environments
- Kokkos-specific accelerator variants
- site-specific compiler externals
- private mirrors or credentials
- broad optional developer package sets unrelated to the immediate PETSc/HPDDM
  workflow

Notes:

- this environment now uses the repository-owned local overlay under
  `packaging/spack/repo/` for narrowly scoped package fixes needed by the
  shared env
- the generated environment view is intentionally kept outside the repository
  checkout via Spack's user cache path so CMake exports do not capture
  source-tree-prefixed include directories
- use the same `SPACK_USER_CONFIG_PATH` consistently when working with this
  env; if Spack falls back to `~/.spack`, it may reuse an older builtin repo
  cache that only exposes PETSc/SLEPc `3.23.x`
- set `SPACK_USER_CACHE_PATH` before activating or installing the env if you
  want that generated view and related Spack caches on fast local storage
- the current overlay fixes `cln` by selecting the upstream CMake build system
  for `@1.3:` instead of the failing autoreconf-based path
- the current overlay also fixes `ginac` fetches by using a reachable source
  mirror while keeping the CMake-based package model
- the current overlay keeps `petsc+hpddm` on the full PETSc install target
  rather than `install-lib`, which is still the safer path for HPDDM-enabled
  shared environments
- `PCHPDDM` runtime support also requires `slepc+hpddm`; upstream PETSc 3.25
  expects `libhpddm_petsc` to come from the SLEPc side when PETSc itself was
  not configured with `--download-slepc`
- the shared manifest pins `petsc`, `slepc`, `py-petsc4py`, and
  `py-slepc4py` to `3.25.0` so the C/C++ and Python solver stack stays on one
  coherent release line
- Feel++ CTest runs should propagate `PETSC_DIR` and `SLEPC_DIR` from the
  configured CMake paths so PETSc can locate `libhpddm_petsc` at runtime
- this environment still references the external `numpex/spack.numpex`
  repository as a transitional package source for `feelpp`
- for Spack `v1.x`, that repository must be configured as a named Git-based
  repo, not as a plain list entry
- the community `spack/spack-packages` repository should not be listed here as
  a raw repo URL because Spack already provides it as the builtin package repo
- that repository reference is intended to be temporary until Feel++ either
  has an in-repo incubation overlay or is upstreamed to `spack/spack-packages`
- site-local compiler or mirror preferences should be expressed through local
  Spack configuration or example files under `packaging/spack/includes/site/`
