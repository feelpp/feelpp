import feelpp
import feelpp.toolboxes as *
from ._solid import *

_csms={
    'solid(2,1)':Solid_2DP1,
    'solid(2,2)':Solid_2DP2,
    'solid(3,1)':Solid_3DP1,
    'solid(3,2)':Solid_3DP2,
}

def solid( dim=2, orderDisp=1, worldComm=None):
    """create a solid toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderDisp -- the polynomial order for the displacement (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    if worldComm is None:
        worldComm=feelpp.Environment.worldCommPtr()
    key='solid('+str(dim)+','+str(orderDisp)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _csms:
        raise RuntimeError('Solid solver '+key+' not existing')
    return _csms[key]( "solid", "solid", worldComm )
