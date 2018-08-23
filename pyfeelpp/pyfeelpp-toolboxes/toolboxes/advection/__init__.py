from pyfeelpp import core
from pyfeelpp.toolboxes.modelcore import *
from _advection import *

_advections={
    'advection(2,1)':Advection_2DP1,
    'advection(2,2)':Advection_2DP2,
    'advection(3,1)':Advection_3DP1,
    'advection(3,2)':Advection_3DP2,
}

def advection( dim=2, order=1, buildMesh=True, worldComm=core.Environment.worldCommPtr() ):
    """create a advection toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    key='advection('+str(dim)+','+str(order)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _advections:
        raise RuntimeError('Advection solver '+key+' not existing')
    return _advections[key]( "advection", buildMesh, worldComm )
