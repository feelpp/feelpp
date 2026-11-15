# Feel++ CPU/OpenMPI 5 Kokkos environment

This environment extends [cpu/openmpi5](../openmpi5/README.md) with Kokkos
5.1.0 and Kokkos Kernels 5.1.0, using Serial and OpenMP host execution.
CUDA, HIP and SYCL are disabled in this environment.

It retains the OpenMPI 5 baseline's UCX/OFI fabrics, Slurm support, external
PMIx, repaired `feelpp.prrte` runtime, preCICE 3.4.1, PETSc/petsc4py 3.25.6,
SLEPc/slepc4py 3.25.2, and pinned builtin package catalog
`093bed090925a0c327ba4fd927f1f1d5e36f2b3e`.
The Kokkos versions are pinned for reproducibility.
Kokkos uses C++20; the consuming Feel++ build retains its own
C++ standard.

Compiler externals, site fabric providers and mirrors remain in local Spack
configuration. The ordinary `cpu/openmpi5` environment remains independent of
Kokkos.

From the repository root, with Spack and `fpp-spack` available:

```bash
export FEELPP_REPO_ROOT="$(git rev-parse --show-toplevel)"
eval "$(fpp-spack activate cpu/openmpi5-kokkos)"
spack concretize --reuse
spack install --fail-fast kokkos kokkos-kernels
spack find -lv kokkos kokkos-kernels openmpi petsc slepc
mpiexec --version
```

`fpp-spack` honors `SPACK_USER_CONFIG_PATH` and `SPACK_USER_CACHE_PATH` when
they are set; otherwise it selects defaults scoped to this checkout. Keep the
config path consistent with the provisioned CPU baseline. Select a dedicated
cache before activation when another checkout already uses these view names,
so its generated views remain independent. A dedicated NVMe cache is suitable.
Before runtime checks, remove inherited library paths for another MPI stack;
verify that the selected launcher reports OpenMPI 5.
For the operational CLI bootstrap, see [ops](../../../../../ops/README.md).

When provisioning a fresh store, install the entire declared environment:

```bash
spack install --fail-fast
```

Installing only the Kokkos roots does not provision missing MPI, PETSc, Python
or preCICE roots. Check their installed status before configuring Feel++.

The generated view is
`$SPACK_USER_CACHE_PATH/views/feelpp/cpu-openmpi5-kokkos`. Package prefixes and
CMake export locations can be inspected independently of a Feel++ build:

```bash
spack location -i kokkos
spack location -i kokkos-kernels
spack location -i openmpi
```

Retain the generated `spack.lock` with the recorded compiler, site configuration
and Spack revision when reproducing a validation run. Machine-specific locks
are ignored here; shared source pins live in the manifest. Reuse an existing
compatible lock instead of routinely forcing reconcretization.

This change provisions dependencies. Feel++ Kokkos detection, installed-package
exports and smoke targets are integrated separately. Use
[hybrid/kokkos-rocm](../../hybrid/kokkos-rocm/README.md) for HIP dependencies.
