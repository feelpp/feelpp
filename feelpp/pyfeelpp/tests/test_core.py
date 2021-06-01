import feelpp 
import sys
import pytest



def test_core(init_feelpp):
    if feelpp.Environment.isMasterRank():
        print("pid:",feelpp.Environment.worldComm().localRank() )
        print("isMasterRank:", feelpp.Environment.isMasterRank())

