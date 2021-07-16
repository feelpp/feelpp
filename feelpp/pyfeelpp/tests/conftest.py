from logging import getLogger
import sys

import py
import pytest
import feelpp

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
    def __init__(self,config):
        try:
            sys.argv=['test_feelpp']
            self.e = feelpp.Environment(sys.argv,config=config)
        except Exception:
            return 


@pytest.fixture(scope="session")
def init_feelpp():
    return InitFeelpp(feelpp.globalRepository("pyfeelpp-tests"))

@pytest.fixture(scope="session")
def init_feelpp_config_local():
    return InitFeelpp(feelpp.localRepository("feelppdb"))
