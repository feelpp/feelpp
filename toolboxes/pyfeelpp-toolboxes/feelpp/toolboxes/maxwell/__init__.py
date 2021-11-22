from pyfeelpp import core
from pyfeelpptoolboxes.modelcore import *
from _maxwell import *

_maxwells={
    'maxwell(2,1)':Maxwell_2DP1,
    'maxwell(3,1)':Maxwell_3DP1,
}

def maxwell( dim=2, order=1, buildMesh=True, worldComm=core.Environment.worldCommPtr() ):
    """create a maxwell toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the potential vector (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    key='maxwell('+str(dim)+','+str(order)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _maxwells:
        raise RuntimeError('Maxwell solver '+key+' not existing')
    return _maxwells[key]( "maxwell", buildMesh, worldComm )
