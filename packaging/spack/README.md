# Spack Packaging Metadata

`packaging/spack/` is the repository-owned home for Feel++ Spack support.

Current scope:

- keep shared Spack environment manifests under version control
- expose one supported shared CPU/OpenMPI environment for Feel++ development
- preserve the imported `openmpi4/` manifest as a visible legacy bootstrap
  reference
- reserve separate locations for:
  - shared environments
  - site-local include templates
  - a temporary in-repo Spack overlay repository if one is needed later

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
  - the first supported target is `spack:openmpi`

Current environment roles:

- `environments/cpu/openmpi/`
  - the first supported shared environment
  - focused on the immediate CPU/OpenMPI PETSc/HPDDM workflow
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
