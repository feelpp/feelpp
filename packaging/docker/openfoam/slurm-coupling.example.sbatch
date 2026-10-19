#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: LGPL-3.0-or-later
#SBATCH --job-name=feelpp-openfoam
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00

# Template only: validate each image's MPI bootstrap before running a coupling.
# Supply the current account/partition with sbatch; adapt ranks and time above.
# Both participant scripts launch ONE rank; they must not invoke mpiexec/srun.
set -euo pipefail
if (( BASH_VERSINFO[0] < 5 || ( BASH_VERSINFO[0] == 5 && BASH_VERSINFO[1] < 1 ) )); then
    echo 'This template requires Bash 5.1 or newer.' >&2
    exit 1
fi
: "${SLURM_JOB_NODELIST:?Run this script inside a two-node Slurm allocation}"
: "${CASE_ROOT:?Set the absolute shared case directory}"
: "${FEELPP_SIF:?Set the absolute Feel++ or KUB SIF path}"
: "${OPENFOAM_SIF:?Set the absolute OpenFOAM SIF path}"
: "${FEELPP_MPI_MODE:?Set an srun MPI mode validated with the Feel++ SIF}"
: "${OPENFOAM_MPI_MODE:?Set an srun MPI mode validated with the OpenFOAM SIF}"
apptainer_command=${APPTAINER_COMMAND:-apptainer}
feelpp_ranks=${FEELPP_RANKS:-32}
openfoam_ranks=${OPENFOAM_RANKS:-32}
case_root=$(realpath "$CASE_ROOT")
feelpp_sif=$(realpath "$FEELPP_SIF")
openfoam_sif=$(realpath "$OPENFOAM_SIF")
for input in "$feelpp_sif" "$openfoam_sif" "$case_root/precice-config.xml" \
    "$case_root/run-feelpp.sh" "$case_root/run-openfoam.sh"; do
    test -r "$input"
done
test -d "$case_root/exchange" && test -w "$case_root/exchange"
mapfile -t nodes < <(scontrol show hostnames "$SLURM_JOB_NODELIST")
if (( ${#nodes[@]} != 2 )); then
    echo 'This example requires exactly two allocated nodes.' >&2
    exit 1
fi
log_dir="$case_root/logs/$SLURM_JOB_ID"
mkdir -p "$log_dir"
srun_pids=()
cleanup()
{
    if (( ${#srun_pids[@]} )); then
        kill "${srun_pids[@]}" 2>/dev/null || true
    fi
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
step_options=(--nodes=1 --exact --exclusive --cpus-per-task=1
    --cpu-bind=verbose,cores --kill-on-bad-exit=1 --export=ALL --label)

# Distinct nodes avoid resource competition and leave two MPI_COMM_WORLD groups.
# Host networking and MPI/PMIx environment variables are retained intentionally.
srun "${step_options[@]}" --nodelist="${nodes[0]}" --ntasks="$feelpp_ranks" \
    --ntasks-per-node="$feelpp_ranks" --mpi="$FEELPP_MPI_MODE" \
    "$apptainer_command" exec --bind "$case_root:/coupling" \
    --pwd /coupling/kub "$feelpp_sif" \
    bash -lc 'unset BASH_ENV; exec bash /coupling/run-feelpp.sh' \
    >"$log_dir/feelpp.log" 2>&1 &
srun_pids+=("$!")
srun "${step_options[@]}" --nodelist="${nodes[1]}" --ntasks="$openfoam_ranks" \
    --ntasks-per-node="$openfoam_ranks" --mpi="$OPENFOAM_MPI_MODE" \
    "$apptainer_command" exec --bind "$case_root:/coupling" \
    --pwd /coupling/openfoam "$openfoam_sif" \
    bash /coupling/run-openfoam.sh >"$log_dir/openfoam.log" 2>&1 &
srun_pids+=("$!")

# Stop the peer if the first completed participant fails instead of leaving it waiting.
first_status=0
wait -n -p completed_pid "${srun_pids[@]}" || first_status=$?
if [[ -n ${completed_pid:-} ]] && (( first_status != 0 )); then
    exit "$first_status"
fi
# If both steps already finished, wait -n reports no eligible job (127).
# Explicit waits below still recover each completed step's exit status.
for participant_pid in "${srun_pids[@]}"; do
    if [[ $participant_pid != ${completed_pid:-} ]]; then
        wait "$participant_pid"
    fi
done
srun_pids=()
