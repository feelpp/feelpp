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

- `prrte`
  - restores three binding-warning messages missing from PRRTE 4.1.0
  - regenerates its prebuilt C message table without purging warning topics
  - prevents a launcher crash when Docker denies NUMA memory binding
  - preserves CPU binding and reports the actual memory-binding warning
- `precice`
  - adds the checksum-pinned 3.4.1 release to the builtin recipe
  - builds the version probe and disables the upstream test-suite build
  - used by the CPU OpenMPI 5 and MPICH environments with MPI/PETSc enabled
    and preCICE Python actions disabled
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
- `opencascade`
  - keeps the builtin OpenCASCADE recipe
  - passes Tcl's library directory explicitly on Darwin so OCCT can configure
    against Spack-built Tcl

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

The CPU OpenMPI 5 and MPICH manifests pin the builtin repository to
`093bed090925a0c327ba4fd927f1f1d5e36f2b3e`. Spack v1.0.0's default
`releases/v2025.07` catalog lacks their PETSc 3.25 Python recipes. The pin uses
a separate checkout below `$user_cache_path/package_repos/`, preserving other
user repositories. The preCICE overlay supplies 3.4.1 on top of that catalog.
The two manifests explicitly select PETSc/petsc4py 3.25.6 and SLEPc/slepc4py
3.25.2, keeping the solver stack on the 3.25 release series. The PETSc overlay
continues to install the HPDDM runtime library.
