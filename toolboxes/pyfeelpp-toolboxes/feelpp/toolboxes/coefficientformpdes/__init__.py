import feelpp
from feelpp.toolboxes import *
from ._cfpdes import *

_cfpdes={
    'cfpdes(2)':cfpdes_2D,
    'cfpdes(3)':cfpdes_3D
}

def cfpdes( dim=2, worldComm=None ):
    """create a cfpde toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    worldComm -- the parallel communicator for the mesh (default: feelpp.Environment::worldCommPtr())
    """
    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()
    key='cfpdes('+str(dim)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _cfpdes:
        raise RuntimeError('cfpde solver '+key+' not existing')
    return _cfpdes[key]( "cfpdes", "cfpdes", worldComm )

def simulate(dim=2):
    f=cfpdes(dim=dim)
    f.init()
    # f.printAndSaveInfo()
    f.solve()
    f.exportResults()


    