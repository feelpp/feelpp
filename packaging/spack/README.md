# Spack Packaging Metadata

`packaging/spack/` is the repository-owned home for Feel++ Spack support.

Current scope:

- keep shared Spack environment manifests under version control
- expose supported shared CPU MPI environments for Feel++ development
- preserve the imported `openmpi4/` manifest as a visible legacy bootstrap
  reference
- reserve separate locations for:
  - shared environments
  - site-local include templates
- a small in-repo Spack overlay repository for package fixes

Directory roles:

- `environments/`
  - shared, repository-owned environment manifests
  - these are the manifests that Feel++ can document and validate in CI
- `includes/site/`
  - examples and templates for site-specific overrides
  - these should not become the default shared repository configuration
- `repo/`
  - optional temporary Spack overlay repository
  - this is for incubation only, not the long-term target for the Feel++ package

Design rules:

- do not keep supported shared metadata under a hidden `.spack/` directory
- do not commit user-local Spack state here
- do not commit generated lockfiles here until the repository decides to manage
  them intentionally
- do not treat an in-repo overlay repository as the final public home of the
  Feel++ package
- keep generated environment views outside the repository checkout

The long-term operational interface for this metadata is:

- `fpp-pkg spack ...`
- `fpp-spack ...`

Current CI-oriented image generation support:

- `fpp-pkg image targets`
  - list the OCI image targets declared in `.github/plan-ci.json`
- `fpp-pkg image bake --target spack:openmpi`
  - generate a self-contained Docker context plus `docker-bake.json`
- `fpp-spack image targets`
  - list the Spack-backed OCI image targets declared in `.github/plan-ci.json`
- `fpp-spack image bake --target spack:openmpi`
  - generate a self-contained Docker context plus `docker-bake.json`
  - the generated image installs the repository-owned Spack environment inside
    an Ubuntu 24.04 base image by default
  - `spack:openmpi` remains the full-stack baseline
- `fpp-pkg image bake --target spack:openmpi5 --component feelpp --from-image ghcr.io/feelpp/feelpp-env:spack-openmpi5`
  - generate a Feel++ core builder and runtime image from the OpenMPI 5 Spack
    environment, without building Toolboxes or MOR
  - `toolboxes` and `mor` can also be selected explicitly; each uses the
    preceding component runtime image as its base
- `fpp-pkg image bake --target spack:mpich --component env`
  - generate the MPICH 3.4.3/CH3 environment image; build jobs default to 15
    and concurrent packages to three
- `fpp-pkg image bake --target spack:mpich --component feelpp --from-image ghcr.io/feelpp/feelpp-env:spack-mpich`
  - generate only the Feel++ core builder and runtime image from that environment

The Feel++ CI workflow keeps the full Spack build by default. To request the
OpenMPI 5 core component, dispatch `.github/workflows/ci.yml` with
`targets=spack:openmpi5`, `only=feelpp`, and `spack_components=true`. Set
`mode=components` explicitly when combining it with other targets. The
`spack:openmpi` target remains full-build only.
For the MPICH core, dispatch with `targets=spack:mpich`, `only=feelpp`, and
`spack_components=true`. The MPICH target has no full, Testsuite, Quickstart,
Toolboxes, or MOR image jobs.

The OCI environment sets `SPACK_USER_CONFIG_PATH=/opt/spack-user` and
`SPACK_USER_CACHE_PATH=/opt/spack-user-cache` inside the image. Locally
configured host paths are used by local Spack commands; Docker does not mount
those directories into this image build.

Current environment roles:

- `environments/cpu/openmpi/`
  - the first supported shared environment
  - focused on the immediate CPU/OpenMPI PETSc/HPDDM workflow
- `environments/cpu/openmpi5/`
  - portable OpenMPI 5 validation environment for Gaya and LUMI-C
  - enables UCX and OFI; LUMI supplies its HPE libfabric/CXI provider through
    site-local Spack configuration
- `environments/cpu/mpich/`
  - MPICH 3.4.3 with CH3 and a unified MPI dependency graph; see its README
- `environments/cpu/openmpi-macosx/`
  - macOS CPU/OpenMPI environment without UCX/CMA
  - uses OpenBLAS for BLAS/LAPACK to avoid the macOS Accelerate/MUMPS crash path
- `environments/openmpi4/`
  - imported from the old hidden path during Phase 0
  - retained as legacy bootstrap material
  - not the recommended baseline for new shared work

Local I/O note:

- the shared environments place their generated view under Spack's user cache
  path, not under the repository checkout
- for local NVMe-backed work, set `SPACK_USER_CACHE_PATH` to a directory on
  fast storage before activating or installing the environment
- `SPACK_USER_CONFIG_PATH` controls the user config scope, while
  `SPACK_USER_CACHE_PATH` controls caches, cloned package repos, and generated
  env views
- keep `SPACK_USER_CONFIG_PATH` stable across concretize/install/build runs;
  otherwise Spack may silently switch back to `~/.spack` package-repo caches
  and expose an older builtin package set than the one validated for Feel++
- backend-specific runtime variables such as `PETSC_DIR` and `SLEPC_DIR` may
  still need to be propagated by the consuming build or test workflow
