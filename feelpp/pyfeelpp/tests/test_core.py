from mpi4py import MPI
import feelpp 
import sys
import pytest



def test_core(init_feelpp):
    feelpp.Environment.changeRepository(directory="pyfeelpp-tests/core/test_core")
    if feelpp.Environment.isMasterRank():
        print("pid:",feelpp.Environment.worldComm().localRank() )
        print("isMasterRank:", feelpp.Environment.isMasterRank())


def test_mpi_bcast(init_feelpp):
    feelpp.Environment.changeRepository(directory="pyfeelpp-tests/core/test_core_bcast")
    
    if feelpp.Environment.isMasterRank():
        data={"key":"test"}
    else:
        data = None
    data=feelpp.Environment.worldComm().localComm().bcast(data,root=0)
    assert(data["key"] == "test")

def test_mpi_numpy(init_feelpp):
    import numpy as np
    if feelpp.Environment.isMasterRank():
        data = np.arange(100, dtype='i')
    else:
        data = np.empty(100, dtype='i')
    feelpp.Environment.worldComm().to_comm().Bcast(data, root=0)
    for i in range(100):
        assert(data[i] == i)

#def test_config_local(init_feelpp_config_local):
#    feelpp.Environment.changeRepository(
#        directory="pyfeelpp-tests/core/test_config_local")

