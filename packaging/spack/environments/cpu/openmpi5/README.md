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

Before installing on LUMI-C, verify that the site configuration exposes the
CXI provider and that Spack selects the intended external libfabric. Runtime
validation should include at least:

- `fi_info -p cxi`
- `ompi_info --param mtl ofi --level 9`
- one multi-node MPI smoke test through the site-supported Slurm launch path

No lockfile is committed initially. Generate it only after the environment has
been concretized and validated on the intended Spack version.
