from logging import getLogger
import sys

import py
import pytest
import feelpp.core as fppc
#import fppc.toolboxes.core as tb
import gmsh

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
            self.feelpp_env = fppc.Environment(sys.argv,config=config)#, 
                                                       #opts= tb.toolboxes_options("heat"))
        except Exception:
            return 


@pytest.fixture(scope="session")
def init_feelpp():
    return InitFeelpp(fppc.globalRepository("pyfeelpp-tests")).feelpp_env

@pytest.fixture(scope="session")
def init_feelpp_config_local():
    return InitFeelpp(fppc.localRepository("feelppdb")).feelpp_env


