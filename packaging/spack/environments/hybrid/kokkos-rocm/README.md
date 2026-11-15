# Feel++ OpenMPI 5 and Kokkos HIP environment

This environment extends [cpu/openmpi5](../../cpu/openmpi5/README.md) with
Kokkos 5.1.0 and Kokkos Kernels 5.1.0. Kokkos enables Serial, OpenMP and HIP;
both Kokkos roots require the HIP-enabled dependency. OpenMPI remains on the
OpenMPI 5 UCX/OFI/Slurm baseline. PETSc 3.25.6 and SLEPc 3.25.2 remain CPU
solver dependencies; device PETSc assembly and GPU-aware MPI are separate work.

Provision the CPU dependency environment first, then concretize this environment
with `--reuse` so the installed OpenMPI 5 solver stack can be shared.

The shared manifest contains no machine-specific compiler externals, ROCm
prefix or GPU architecture. Configure those before concretizing. The
[MI210 example](../../../includes/site/kokkos-rocm-gfx90a.example.yaml) selects
`gfx90a` and a coherent Spack-built ROCm 7.2.0 toolchain, matching the existing
prototype's local ROCm alternative. Adjust this policy for other hardware.
An external ROCm installation must be complete and version compatible; a
recorded prefix alone does not establish that its compiler still exists.

From the repository root, with Spack and `fpp-spack` available:

```bash
export FEELPP_REPO_ROOT="$(git rev-parse --show-toplevel)"
eval "$(fpp-spack activate hybrid/kokkos-rocm)"
export FEELPP_ROCM_SITE="$SPACK_USER_CACHE_PATH/site/kokkos-rocm"
mkdir -p "$FEELPP_ROCM_SITE"
cp packaging/spack/includes/site/kokkos-rocm-gfx90a.example.yaml "$FEELPP_ROCM_SITE/packages.yaml"
spack -C "$FEELPP_ROCM_SITE" concretize --reuse
spack -C "$FEELPP_ROCM_SITE" install --fail-fast kokkos kokkos-kernels rocminfo
spack find -lv kokkos kokkos-kernels hip llvm-amdgpu rocminfo openmpi
```

The copy command initializes a new local policy directory. Preserve and review
an existing site policy instead of overwriting its changes. Keep
`SPACK_USER_CONFIG_PATH` consistent with the CPU baseline; use a dedicated
cache for each checkout. See [ops](../../../../../ops/README.md) for CLI setup.
For a fresh dependency store, provision all declared roots with
`spack -C "$FEELPP_ROCM_SITE" install --fail-fast`.

A Spack-built ROCm toolchain can compile on a host without a visible GPU.
Keep the environment view activated when consuming its CMake packages: it
provides transitive package lookup paths, including rocPRIM. Setting
`Kokkos_DIR` alone does not provide those transitive paths.
Load HIP's package environment so its split compiler/runtime prefixes are
configured correctly:

```bash
eval "$(spack load --sh hip)"
hipcc --version
```

On an allocated compatible GPU node, check the actual device architecture:

```bash
rocminfo
```

For the MI210 policy, the reported GPU must include `gfx90a`. Record the
allocation, compiler version, architecture, resolved specs and package
prefixes with the lockfile. This checks dependency/toolchain readiness; it does
not establish that a Feel++ FE kernel or a GPU solve has run.

The generated view is
`$SPACK_USER_CACHE_PATH/views/feelpp/hybrid-kokkos-rocm`, separate from the host
Kokkos view and ordinary CPU baseline. Retain the generated local
`spack.lock` for exact reproduction; machine-specific locks are ignored here.
Feel++ compiler integration and runtime ownership are implemented separately.
