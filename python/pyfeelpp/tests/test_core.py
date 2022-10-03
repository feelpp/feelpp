from mpi4py import MPI
import feelpp 
import sys
import pytest
import os,shutil
from pathlib import Path
import subprocess

config_cases=[  ("pyfeelpp-tests/core/test_config",feelpp.Location.standard,False),
                ("pyfeelpp-tests/core/test_config",feelpp.Location.relative,True),
                #("pyfeelpp-tests/core/test_config",feelpp.Location.git,True),
                ("/tmp/toto/pyfeelpp-tests/core/test_config",feelpp.Location.absolute,True),
            ]

@pytest.mark.parametrize("dir,location,rm", config_cases)
def test_config(init_feelpp,dir,location,rm):
    e=init_feelpp
    if location == feelpp.Location.git :
        if e.isMasterRank():
            shutil.rmtree("/tmp/test_config")
            os.mkdir("/tmp/test_config")
            os.chdir("/tmp/test_config")
            subprocess.call(["git", "init"])
            os.mkdir("/tmp/test_config/tutu")
        e.worldComm().globalComm().Barrier()
        os.chdir("/tmp/test_config/tutu")
    
    e.changeRepository(directory=dir,location=location)
    thedir=Path(feelpp.Environment.rootRepository()) / Path(dir)
    assert thedir.is_dir()
    assert (Path(os.getcwd())/Path("logs")).is_dir()
    e.worldComm().globalComm().Barrier()
    if rm and e.isMasterRank():
        os.chdir(Path(feelpp.Environment.rootRepository()).parent)
        shutil.rmtree(feelpp.Environment.rootRepository())


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

def test_worldcomm_split(init_feelpp):
    e = init_feelpp
    if e.numberOfProcessors() > 1 and e.numberOfProcessors()%2 == 0:
        color,w,wglob=e.worldCommPtr().split(2)
        assert(wglob.globalSize() == e.numberOfProcessors())
        assert(wglob.localSize() == e.numberOfProcessors()/2)
        assert(w.localSize() == e.numberOfProcessors()/2)
        assert(w.globalSize() == e.numberOfProcessors()/2)
          


#def test_config_local(init_feelpp_config_local):
#    feelpp.Environment.changeRepository(
#        directory="pyfeelpp-tests/core/test_config_local")

def test_config_parser(init_feelpp):
    e=init_feelpp
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/core/test_config_parser")
    config = feelpp.readcfg(os.path.dirname(__file__)+'/test.cfg')
    print("sections: {}".format(config.sections()))
    d = config['feelpp']['directory']
    #assert(d == 'toolboxes/fluid/TurekHron/cfd1/P2P1G1')
    j = config['fluid']['filename']
    assert(j == "$cfgdir/cfd1.json")
