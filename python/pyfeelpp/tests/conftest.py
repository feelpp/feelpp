from logging import getLogger
import sys

import py
import pytest
import feelpp
import gmsh

# Attempt to import feelpp.toolboxes.core dynamically
try:
    import feelpp.toolboxes.core as tb
    has_toolboxes_core = True
except ImportError:
    has_toolboxes_core = False

log = getLogger(__name__)
MPI_ARGS = ("mpirun", "-n")
PYTEST_ARGS = (sys.executable, "-mpytest")

@pytest.fixture
def has_mpi4py():
    try:
        import mpi4py
        return True
    except ImportError:
        return False

@pytest.fixture
def has_petsc4py():
    try:
        import petsc4py
        return True
    except ImportError:
        return False

class InitFeelpp:
    def __init__(self, config):
        try:
            sys.argv = ['test_feelpp']
            if has_toolboxes_core:
                # Use toolboxes.core if available
                self.feelpp_env = feelpp.Environment(sys.argv, config=config, opts=tb.toolboxes_options("heat"))
            else:
                # Proceed without toolboxes.core specific functionality
                self.feelpp_env = feelpp.Environment(sys.argv, config=config)
        except Exception as e:
            log.error(f"Failed to initialize Feel++ environment: {e}")
            return None

@pytest.fixture(scope="session")
def init_feelpp():
    return InitFeelpp(feelpp.globalRepository("pyfeelpp-tests")).feelpp_env

@pytest.fixture(scope="session")
def init_feelpp_config_local():
    return InitFeelpp(feelpp.localRepository("feelppdb")).feelpp_env
