import feelpp
import feelpp.toolboxes
from ._fsi import *

_fsis={
    'fsi(2,1)':Fsi_2DP1,
    'fsi(2,2)':Fsi_2DP2,
    'fsi(3,1)':Fsi_3DP1,
    'fsi(3,2)':Fsi_3DP2,
}

def fsi( dim=2, order=1, buildMesh=True, worldComm=feelpp.Environment.worldCommPtr() ):
    """create a fsi toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: feelpp.Environment::worldCommPtr())
    """
    key='fsi('+str(dim)+','+str(orderPotential)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _fsis:
        raise RuntimeError('Fsi solver '+key+' not existing')
    return _fsis[key]( "fsi", buildMesh, worldComm )
