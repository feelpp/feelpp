import argparse
import os
import sys
import traceback

import pytest


def finalize_mpi():
    try:
        from mpi4py import MPI
    except ImportError:
        return

    if not MPI.Is_initialized() or MPI.Is_finalized():
        return

    try:
        MPI.COMM_WORLD.Barrier()
    except Exception:
        pass
    MPI.Finalize()


def main():
    parser = argparse.ArgumentParser(
        description="Run a CFPDE pytest file and exit without Python interpreter teardown."
    )
    parser.add_argument("--pytest-file", required=True)
    args = parser.parse_args()
    return pytest.main(["-s", args.pytest_file])


if __name__ == "__main__":
    status = 0
    try:
        status = main()
    except Exception:
        traceback.print_exc()
        status = 1
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        finalize_mpi()
        os._exit(status)
