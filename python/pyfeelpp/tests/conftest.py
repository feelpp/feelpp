from logging import getLogger
import sys

import py
import pytest
import feelpp.core as fppc
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
        self.feelpp_env = None
        try:
            sys.argv = ['test_feelpp']
            if has_toolboxes_core:
                # Use toolboxes.core if available
                self.feelpp_env = fppc.Environment(sys.argv, config=config, opts=tb.toolboxes_options("heat"))
            else:
                # Proceed without toolboxes.core specific functionality
                self.feelpp_env = fppc.Environment(sys.argv, config=config)
            log.info("Feel++ environment initialized successfully")
        except Exception as e:
            log.error(f"Failed to initialize Feel++ environment: {e}")
            # Try a simpler initialization without config
            try:
                self.feelpp_env = fppc.Environment(sys.argv)
                log.info("Feel++ environment initialized with minimal config")
            except Exception as e2:
                log.error(f"Complete failure to initialize Feel++ environment: {e2}")
                self.feelpp_env = None

@pytest.fixture(scope="session")
def init_feelpp():
    init_obj = InitFeelpp(fppc.globalRepository("pyfeelpp-tests"))
    if init_obj.feelpp_env is None:
        log.warning("Feel++ environment not properly initialized")
        # Try to create a minimal environment for testing
        try:
            minimal_env = fppc.Environment(['test_feelpp'])
            return minimal_env
        except Exception as e:
            log.error(f"Could not create minimal Feel++ environment: {e}")
            pytest.skip("Feel++ environment could not be initialized")
    return init_obj.feelpp_env

@pytest.fixture(scope="session")
def init_feelpp_config_local():
    init_obj = InitFeelpp(fppc.localRepository("feelppdb"))
    if init_obj.feelpp_env is None:
        log.warning("Feel++ environment not properly initialized")
        # Try to create a minimal environment for testing
        try:
            minimal_env = fppc.Environment(['test_feelpp'])
            return minimal_env
        except Exception as e:
            log.error(f"Could not create minimal Feel++ environment: {e}")
            pytest.skip("Feel++ environment could not be initialized")
    return init_obj.feelpp_env
