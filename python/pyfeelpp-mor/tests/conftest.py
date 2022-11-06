from logging import getLogger
import sys, shutil

import py
import pytest
import feelpp
import feelpp.mor as mor
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
    def __init__(self, config):

        try:
            shutil.rmtree(feelpp.Environment.rootRepository() + '/pyfeelppmor-tests')
            print(f'Directory {feelpp.Environment.rootRepository()}/pyfeelppmor-tests removed')
        except FileNotFoundError:
            print(f"Deletion of {feelpp.Environment.rootRepository()}/pyfeelppmor-tests did not succeded : Directory doesn't exist")

        try:
            print('xxx call init_feelpp;', __file__)
            sys.argv=['test_pyfeelppmor']
            self.e = feelpp.Environment(
                sys.argv, opts= feelpp.backend_options("Iv")
                                .add(tb.toolboxes_options("heat"))
                                .add(mor.makeToolboxMorOptions()),
                config=config
                                )
            print('is master? ', feelpp.Environment.worldCommPtr().isMasterRank())
        except Exception as err:
            print('Exception caught while initializing Feel++: '.format(err))
            return


@pytest.fixture(scope="session")
def init_feelpp():
    return InitFeelpp(feelpp.globalRepository("pyfeelppmor-tests"))

@pytest.fixture(scope="session")
def init_feelpp_config_local():
    return InitFeelpp(feelpp.localRepository("feelppdb"))
