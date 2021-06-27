from logging import getLogger
import sys

import py
import pytest
import feelpp
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



class InitFeelpp:
    def __init__(self):
        try:
            print('xxx call init_feelpp;', __file__)
            sys.argv=['test_pyfeelpptoolboxes']
            self.e = feelpp.Environment(
                sys.argv, opts=tb.toolboxes_options("electric")
                                .add(tb.toolboxes_options("fluid"))
                                .add(tb.toolboxes_options("heat"))
                                .add(tb.toolboxes_options("coefficient-form-pdes", "cfpdes"))
                                )
            print('is master? ', feelpp.Environment.worldCommPtr().isMasterRank())
        except Exception as err:
            print('Exception caucht while initializing Feel++: '.format(err))
            return 


@pytest.fixture(scope="session", autouse=True)
def init_feelpp():
    return InitFeelpp()
