from logging import getLogger
import sys

import py
import pytest
import feelpp
import feelpp.mor as fppmor
import feelpp.toolboxes.core as tb
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
            print('xxx call init_feelpp;', __file__)
            sys.argv=['test_pyfeelppmor']
            self.e = fppc.Environment(
                sys.argv, opts= fppc.backend_options("Iv")
                                .add(tb.toolboxes_options("heat"))
                                .add(tb.toolboxes_options("fluid"))
                                .add(fppmor.makeToolboxMorOptions()),
                config=config
                                )
            print('is master? ', fppc.Environment.worldCommPtr().isMasterRank())
        except Exception as err:
            print('Exception caught while initializing Feel++: '.format(err))
            return 


@pytest.fixture(scope="session")
def init_feelpp():
    return InitFeelpp(fppc.globalRepository("pyfeelppmor-tests"))

@pytest.fixture(scope="session")
def init_feelpp_config_local():
    return InitFeelpp(fppc.localRepository("feelppdb"))
