from pyfeelpp import core
from pyfeelpp.toolboxes.modelcore import *
from _hdg import *

_hdgs={
    'hdg(2,1)':HDGPoisson_2DP1,
    'hdg(2,2)':HDGPoisson_2DP2,
    'hdg(3,1)':HDGPoisson_3DP1,
    'hdg(3,2)':HDGPoisson_3DP2,
}

def hdg( dim=2, order=1, buildMesh=True, worldComm=core.Environment.worldCommPtr() ):
    """create a hdg toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the fields : potential, flux, displacement, stress and associated traces (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    key='hdg('+str(dim)+','+str(order)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _hdgs:
        raise RuntimeError('HDG solver '+key+' not existing')
    return _hdgs[key]( "hdg", buildMesh, worldComm )
