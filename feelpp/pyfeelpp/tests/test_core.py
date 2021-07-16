import feelpp 
import sys
import pytest



def test_core(init_feelpp):
    feelpp.Environment.changeRepository(directory="pyfeelpp-tests/core/test_core")
    if feelpp.Environment.isMasterRank():
        print("pid:",feelpp.Environment.worldComm().localRank() )
        print("isMasterRank:", feelpp.Environment.isMasterRank())


def test_config_local(init_feelpp_config_local):
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/core/test_config_local")

