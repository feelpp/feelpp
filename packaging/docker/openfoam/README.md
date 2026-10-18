# OpenFOAM 2606 / preCICE recipes

This directory owns the source recipes for a separate OpenFOAM participant image.
KUB continues to build from the Feel++ core image with preCICE support. The two
images communicate through preCICE sockets; each participant uses its own MPI
installation. No OpenFOAM layer is added to the Feel++ environments or KUB image.

## Versions and scope

| Dependency | Selection |
|---|---|
| Base | Ubuntu 24.04 Noble, `linux/amd64` initially |
| OpenFOAM repository | <https://dl.openfoam.com/repos/deb>, signed `noble/main` |
| Requested packages | `openfoam2606-default` **and** `openfoam2606-dev` |
| Exact package version | `2606.0~rc2-1`, observed in the repository on 2026-10-11 |
| preCICE | Official Noble `libprecice3` package, **3.4.1**, verified SHA-256 |
| OpenFOAM adapter | **1.4.0**, source archive with verified SHA-256 |
| MPI | Noble's system OpenMPI for OpenFOAM, preCICE, and the adapter |

The repository currently exposes an **RC2 package version** for 2606. Inspect the
signed APT index before changing that pin. This is not a claim that the package
is a final 2606 release. Other distributions and arm64 remain later work; this
recipe does not change the platforms of existing Feel++ images.

The adapter supports OpenFOAM 2606 and preCICE v3. It is built using `Allwmake`
against the installed OpenFOAM development files and the image's own preCICE.
[Compatibility](https://precice.org/adapter-openfoam-support.html),
[adapter 1.4.0](https://github.com/precice/openfoam-adapter/releases/tag/v1.4.0),
[build instructions](https://precice.org/adapter-openfoam-get.html).

Both requested OpenFOAM packages stay installed, including their source/tools/
tutorial dependencies. This is a development-capable image. The adapter goes
into `/usr/local/lib/openfoam-precice`; the entrypoint sources OpenFOAM's bashrc
and adds that directory to the library path. Build logs, package versions,
archive checksums, and the adapter licence are saved under
`/usr/local/share/openfoam-precice`.

These targets are a standalone Docker Bake recipe. They are **not yet registered**
as `fpp-pkg image bake` targets or in `.github/plan-ci.json`. The existing Feel++
environment/core workflow remains unchanged. Build/publish CI for this separate
image is a subsequent integration step recorded in the overall plan.

## Build and smoke test

From this directory:

```bash
docker buildx bake --print openfoam2606-precice
docker buildx bake --load openfoam2606-precice
image=ghcr.io/feelpp/openfoam:2606-precice3.4.1-ubuntu24.04
docker image inspect "$image" --format '{{.Id}} {{json .RepoDigests}}'
docker run --rm "$image" openfoam-smoke-test --mpi
```

Use `BASE_IMAGE=ubuntu:24.04@sha256:<verified-digest>` to fix the base digest for
the image campaign. Set `BUILD_JOBS` to bound adapter compilation. A package pin
update can be passed with `OPENFOAM_PACKAGE_VERSION` after reviewing the index.
The optional `openfoam2606` target installs only the OpenFOAM packages.

The build checks the adapter with `ldd -r` and a load probe. The runtime smoke
test reports OpenFOAM/preCICE/OpenMPI versions, checks solver and adapter library
resolution, calls solver help, and optionally runs a pinned two-rank MPI probe.
It does not run a CFD or KUB simulation. The `Allwmake` exit code alone is
insufficient: upstream can print linking errors and return success, so the
recipe independently rejects missing libraries and undefined symbols.

## Published artifacts — 2026-10-11

Both formats are available in `ghcr.io/feelpp/openfoam`:

| Format | Tag | Registry manifest digest |
|---|---|---|
| Docker | `2606-precice3.4.1-ubuntu24.04` | `sha256:9a11e4fa1cb4506885b44125e8f58cb1a2c55aa25585d8c0d45e6d422c583627` |
| SIF / ORAS | `2606-precice3.4.1-ubuntu24.04-sif` | `sha256:500aef4a8593ed056ba8044845c11fcc4d854d334b54b45895c3da3f4f38ad0e` |

The same artifacts also have dated tags
`2606-precice3.4.1-ubuntu24.04-20261011-155729` and
`2606-precice3.4.1-ubuntu24.04-20261011-155729-sif`.
The package is **public**. Anonymous Docker and SIF manifest requests return
HTTP 200, and an anonymous SIF download passed the file checksum check on
2026-10-11. No registry login is needed:

```bash
apptainer pull openfoam2606-precice.sif \
    oras://ghcr.io/feelpp/openfoam:2606-precice3.4.1-ubuntu24.04-sif
printf '%s  %s\n' \
    c3b626e5e0161d5b5caa4c52363fe11b4f623b14546e55d9351c9814ae369e70 \
    openfoam2606-precice.sif | sha256sum --check
```

The SIF is **638,783,488 bytes**, built with Apptainer 1.4.3 directly from the
validated local Docker image. It passed the version, `ldd -r`, adapter-load,
solver-help and pinned two-rank MPI gates. Pulling the published ORAS artifact
back produced the identical SHA-256. Use `oras://` for this already converted
SIF; `docker://` selects the Docker image and performs another conversion.
[Apptainer registry protocols](https://apptainer.org/docs/user/latest/registry.html).

For the CPU-only smoke on Gaya, hwloc's optional GL discovery produced X11
authorization messages. Clearing `DISPLAY` alone did not resolve them;
disabling that component did, with all runtime gates still passing:

```bash
apptainer run --cleanenv --env HWLOC_COMPONENTS=-gl \
    openfoam2606-precice.sif openfoam-smoke-test --mpi
```

No missing libraries, undefined symbols or memory-binding warnings occurred in
that SIF run. This is a local runtime gate; Slurm/PMIx and inter-node operation
remain to be validated independently. Logs, exact commands, image/file hashes,
registry manifests and authenticated/anonymous round-trip verification are retained under
`build/openfoam-image-20261011-mkr1TU/publish-20261011-155729/`.

## Compose contract

`compose.example.yaml` is a template for **two actual coupling participants**.
An existing `feelpp_kub_cem` invocation does not become a preCICE participant
merely because its image links Feel++ with preCICE. KUB still needs a coupling
entry point and case configuration. The example therefore requires explicit
commands rather than starting a KUB simulation by default.

Prepare one absolute host case path containing:

```text
case/
  precice-config.xml
  exchange/             # writable by both participants; fresh for each run
  kub/                  # KUB input/output and coupling configuration
  openfoam/             # OpenFOAM case, controlDict and preciceDict
```

Both services mount it at `/coupling`. Use the same preCICE XML, participant and
mesh/data names, and write these socket settings in that XML:

```xml
<m2n:sockets acceptor="KUB" connector="OpenFOAM"
             network="eth0" exchange-directory="/coupling/exchange" />
```

Loopback is not shared by the containers. The bridge interface must be selected
explicitly. TCP ports need no host publication when both containers use the same
Compose network. The shared directory holds discovery files; it is not the
transport for the coupled fields. Use one independent `mpiexec` per participant,
with each image's MPI and binding flags. Cross-image MPI ABI compatibility is not
required by this socket arrangement. [preCICE communication configuration](https://precice.org/configuration-communication.html).

Set `KUB_IMAGE` to the rebuilt image (preferably a digest), `CASE_ROOT`,
`KUB_COMMAND`, and `OPENFOAM_COMMAND`. Optionally set `OPENFOAM_IMAGE`, `RUN_UID`,
and `RUN_GID`. Ensure case permissions match these IDs. Quote shell commands in
the environment file and use `$$` for variables to be expanded inside a
container rather than by Compose.

```bash
docker compose --env-file coupling.env -f compose.example.yaml config
# After both images and a small coupling case have passed their gates:
docker compose --env-file coupling.env -f compose.example.yaml up --abort-on-container-exit
```

No `depends_on` completion dependency is used: both preCICE participants must run
together. Place ranks and give the services disjoint CPU sets for measurements;
record their resolved commands and bindings. For EuroHPC, use Apptainer and
Slurm allocations/job steps instead of Compose, retaining the participant and
socket/shared-path contract with the cluster's reachable network interface.

## Apptainer and Slurm

Apptainer can convert and run the Docker/OCI images; it does **not natively run
`compose.yaml`**. On a supercomputer, Slurm supplies the allocation, concurrent
job steps, placement and lifecycle management. Keep the two participant images.
Neither Docker nor Compose is needed on the compute nodes.

| Compose configuration | Apptainer / Slurm equivalent |
|---|---|
| Two services | Two concurrent Slurm steps, one SIF per participant |
| Service MPI command | Independent MPI launch for each participant |
| Shared bind mount | Bind the same shared case directory at `/coupling` |
| Bridge network | Default Apptainer host network; choose a reachable cluster interface |
| CPU sets | Slurm CPU binding and disjoint resources |
| Stop on participant failure | Terminate the peer step and retain both logs |

### Prepare images and the case

After the Docker smoke tests pass, convert the local images on a build machine
that has Docker and Apptainer. The Feel++ tag below is our locally built core;
replace it with the rebuilt KUB image when testing a KUB coupling participant.

```bash
mkdir -p build/coupling-sif
apptainer build build/coupling-sif/feelpp-openmpi5.sif \
    docker-daemon:ghcr.io/feelpp/feelpp:spack-openmpi5-feature-eurohpc
apptainer build build/coupling-sif/openfoam2606-precice.sif \
    docker-daemon:ghcr.io/feelpp/openfoam:2606-precice3.4.1-ubuntu24.04
sha256sum build/coupling-sif/*.sif
```

Transfer the SIFs and record their hashes. Alternatively, after publishing an
image, use `apptainer pull image.sif docker://ghcr.io/...@sha256:<digest>` on a
transfer node. Conversion from a local daemon and from a registry are supported
by [Apptainer](https://apptainer.org/docs/user/latest/docker_and_oci.html).

Use the same case layout as above, with a fresh writable `exchange/` for every
run. On multiple nodes, that directory must be visible from both participants;
node-local `/tmp` is insufficient. Meshes and bulk output may use participant
scratch storage independently. Select an interface that provides reachable IP
addresses on the allocated nodes, checking `ip -br addr` there. For example,
**only if that cluster actually provides IP over `ib0`**:

```xml
<m2n:sockets acceptor="KUB" connector="OpenFOAM"
             network="ib0" exchange-directory="/coupling/exchange" />
```

Omit Apptainer `--net`: its default host network allows access to the cluster
network. For two participants entirely on one node, loopback is also usable.
For participants on different nodes, select their reachable interface instead
of loopback or a Docker bridge. preCICE uses the shared directory for connection
discovery and TCP for field transfer; it supports different numbers of ranks in
the participants. [Apptainer networking](https://apptainer.org/docs/user/latest/docker_and_oci.html),
[preCICE communication](https://precice.org/configuration-communication.html).

### Validate each MPI launch before coupling

Our Feel++ image uses OpenMPI **5**, while this OpenFOAM image uses Noble's
OpenMPI **4**. The TCP connection between participants does not require MPI ABI
compatibility across those images. Each image's MPI must nevertheless bootstrap
correctly through the cluster launcher and use the site's network devices.
Docker's successful two-rank test does not validate Slurm/PMIx on a cluster.

Inside a valid allocation, inspect `srun --mpi=list`, then test a supported
mode separately for each SIF. An OpenFOAM example, after selecting a compatible
mode in `OPENFOAM_MPI_MODE`, is:

```bash
srun --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpu-bind=verbose,cores \
    --mpi="$OPENFOAM_MPI_MODE" apptainer exec "$OPENFOAM_SIF" \
    /opt/feelpp-openfoam-recipe/entrypoint.sh openfoam-mpi-probe
```

Require `rank=0/1`, `size=2` and `sum=3`, with different hostnames in this
inter-node test. Two independent `rank=0 size=1` outputs are a failed launch.
Repeat with a probe linked to the Feel++ image's MPI; retain `feelpp_info`,
`ompi_info`, `ldd` and OpenFOAM/plugin checks for the SIFs too. A different MPI
mode may be needed for each image. Preserve the Slurm/PMI/PMIx bootstrap
environment; `--cleanenv` is inappropriate unless those variables are restored
explicitly. If direct `srun` launch is incompatible, validate the site's
supported host-launcher/container combination or rebuild against its MPI stack.
[Apptainer MPI models](https://apptainer.org/docs/user/latest/mpi.html),
[Slurm MPI configuration](https://slurm.schedmd.com/mpi_guide.html).

### Launch two participants

[`slurm-coupling.example.sbatch`](slurm-coupling.example.sbatch) provides a
two-node template, with one participant on each node and pinned ranks. It uses
two separate `srun` commands, retains separate logs and stops the peer on an
error. Adapt the rank counts, allocation and time limit to the case. This
template has **not been validated on Karolina or MeluXina**.

Create two scripts in the shared case directory. Each script runs **one rank**,
as Slurm creates the MPI processes. `run-feelpp.sh` must `exec` an actual Feel++
or KUB preCICE participant with its case configuration. The current one-way
Laplacian quickstart can serve as a separate Feel++ correctness example; the
existing KUB CEM executable still needs its coupling implementation.
For `run-openfoam.sh`, source the image's setup through the explicit entrypoint:

```bash
#!/usr/bin/env bash
set -euo pipefail
exec /opt/feelpp-openfoam-recipe/entrypoint.sh \
    "${OPENFOAM_SOLVER:?Set the solver configured for this coupling case}" \
    -parallel -case /coupling/openfoam
```

`apptainer exec` bypasses the Docker entrypoint, which is why this setup call is
explicit. The OpenFOAM case must have the adapter enabled/configured and be
decomposed for `OPENFOAM_RANKS` beforehand. Both scripts must use matching
participant names and `/coupling/precice-config.xml`. Do not put another
`mpiexec` in these scripts, or use a single unsplit MPI world for both solvers.

Set absolute `CASE_ROOT`, `FEELPP_SIF`, `OPENFOAM_SIF`, and validated
`FEELPP_MPI_MODE` / `OPENFOAM_MPI_MODE`; also set `OPENFOAM_SOLVER`. Then submit:

```bash
sbatch --account="$CURRENT_ACCOUNT" --partition="$CURRENT_PARTITION" \
    --export=ALL packaging/docker/openfoam/slurm-coupling.example.sbatch
```

The template defaults to 32 ranks per participant; adapt the SBATCH resource
requests along with `FEELPP_RANKS` / `OPENFOAM_RANKS`. It requires Bash 5.1 or
newer for `wait -n -p`. For more than two nodes, assign disjoint node subsets to
the two steps. Smaller cases can share a node with disjoint CPU allocations;
validate both placements before comparing timings.

The checked `hpc-systems` inventory selects Karolina `eu-26-92 / qcpu` and
MeluXina `p201438 / cpu`, but **both records end on 2026-10-01**. Confirm a
current allocation before using those values; no cluster job was submitted
during this documentation/build repair.
