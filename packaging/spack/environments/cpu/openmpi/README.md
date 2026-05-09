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
- the shared manifest currently also pins `gmsh@4.13.1`; the bundled Spack
  `v1.0.0` builtin recipe otherwise selects `4.15.x`, which is not yet
  validated for this environment
- the shared manifest includes `py-gmsh@4.13.1` so Python tests can import the
  official `gmsh` module, matching the role of `python3-gmsh` in the apt-based
  environments
- the shared manifest also forces `mesa~llvm`; this avoids the current Mesa
  `llvm-config` tool mismatch in the GLX path pulled by `opencascade`
- the shared manifest includes `mesa-glu` explicitly so OpenGL-enabled
  Feel++ binaries can resolve `libGLU.so.1` during full-image CTest/runtime
- the shared manifest includes `pugixml` explicitly so the FMI/XML-related
  CMake checks resolve the same way they do in the apt-based environments
- the shared manifest includes `rsync` explicitly because the shared Feel++
  CMake modules require it during configuration and testcase setup
- Feel++ CTest runs should propagate `PETSC_DIR` and `SLEPC_DIR` from the
  configured CMake paths so PETSc can locate `libhpddm_petsc` at runtime
- this environment intentionally uses only the repository-owned `feelpp`
  overlay plus Spack builtin packages; external third-party package repos are
  excluded here so shared env installs do not accidentally shadow core packages
  such as `gmsh`
- the community `spack/spack-packages` repository should not be listed here as
  a raw repo URL because Spack already provides it as the builtin package repo
- the long-term goal remains the same: keep Feel++-specific packaging in the
  repository-owned overlay until it is either upstreamed or no longer needed
- site-local compiler or mirror preferences should be expressed through local
  Spack configuration or example files under `packaging/spack/includes/site/`
