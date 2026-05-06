# Temporary Spack Overlay Repository

This directory is the repository-owned local Spack overlay used while package
recipes are being stabilized.

Current scope:

- carry minimal overrides needed to unblock the shared Feel++ Spack
  environments
- keep those overrides versioned with the repository instead of relying on
  user-local edits inside a Spack clone
- prefer small, targeted package overrides over a broad fork of builtin
  packages
- follow the Spack package API v2 layout under `spack_repo/<namespace>/`

Current overlays:

- `cln`
  - switches `@1.3:` to the upstream CMake build system
  - keeps autotools available only for older `@:1.2` releases
  - avoids the broken autoreconf path that left
    `AX_CXX_COMPILE_STDCXX([11], ...)` unexpanded in the shared env install
- `ginac`
  - keeps the CMake-based package model from `numpex`
  - uses the Codeberg Git repository after the upstream move in February 2026
  - uses the current official `1.8.10` release tarball from `ginac.de`
  - keeps the reachable MacPorts distfiles mirror as the fallback base for
    older tarball URLs because `https://www.ginac.de/ginac-1.8.8.tar.bz2` now
    returns `404`
- `petsc`
  - keeps the builtin PETSc recipe but overrides the install target for
    `+hpddm`
  - uses the full PETSc `install` target instead of `install-lib` when HPDDM
    is enabled
  - ensures `libhpddm_petsc` is installed so `PCHPDDM` works at runtime in the
    shared Feel++ environment

Layout:

- `spack-repo-index.yaml`
  - advertises the repository roots when this directory is used as a Git-based
    mono-repo later
- `spack_repo/feelpp/`
  - local repository root used by the shared env today

Important rule:

- this overlay is an incubation mechanism
- it is not the intended final public home of the Feel++ Spack package

The long-term target remains an upstream recipe in `spack/spack-packages`.
