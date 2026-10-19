# Feel++ Shared CPU/OpenMPI 5 Environment

This validation environment mirrors `cpu/openmpi` while selecting OpenMPI 5.
It is intended to provide one portable OpenMPI 5 dependency stack for Gaya and
LUMI-C rather than maintaining a separate full manifest for each machine.

The environment enables both communication paths needed by those systems:

- UCX for systems such as Gaya
- OFI/libfabric for HPE Slingshot 11 systems such as LUMI-C

Site-specific compiler, scheduler, libfabric, CXI provider, and mirror settings
must remain in local Spack configuration. In particular, LUMI-C should expose
the HPE-provided libfabric installation with its CXI provider as an external
package instead of rebuilding a generic libfabric inside this manifest.

Activate the environment with:

```console
eval "$(fpp-spack activate openmpi5)"
```

Configure, build, and test in an environment-specific tree with:

```console
cmake --preset release-clang-spack-cpu-openmpi5
cmake --build --preset release-clang-spack-cpu-openmpi5
ctest --preset release-clang-spack-cpu-openmpi5
```

Before installing on LUMI-C, verify that the site configuration exposes the
CXI provider and that Spack selects the intended external libfabric. Runtime
validation should include at least:

- `fi_info -p cxi`
- `ompi_info --param mtl ofi --level 9`
- one multi-node MPI smoke test through the site-supported Slurm launch path

No lockfile is committed initially. Generate it only after the environment has
been concretized and validated on the intended Spack version.

## Solver stack update — 2026-10-11

The manifest selects PETSc/petsc4py 3.25.6, SLEPc/slepc4py 3.25.2 and preCICE
3.4.1. Its pinned builtin catalog is
`093bed090925a0c327ba4fd927f1f1d5e36f2b3e`, checked out separately below
`$SPACK_USER_CACHE_PATH/package_repos/`. The repository's PETSc HPDDM install
overlay remains enabled. These source pins require rebuilding the dependency
image and then the Feel++ core image; they do not update an existing image.

The PRRTE overlay restores binding-warning messages missing from version
4.1.0 and regenerates the release's embedded message table without purging
warning topics. Without that repair, a normal Docker container's denied NUMA
memory-binding request crashes the launcher while formatting its warning.
The manifest selects external PMIx and an explicit `feelpp.prrte` root so
Spack uses the repaired runtime. CPU pinning stays enabled; the launcher
reports the memory-binding restriction.
