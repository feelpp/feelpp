import feelpp
from feelpp.toolboxes.core import *

has_csm = False
_csms = None
try:
    from ._solid import *

    _csms={
        'solid(2,1)': Solid_2DP1,
        'solid(2,2)': Solid_2DP2,
        'solid(2,3)': Solid_2DP3,
        'solid(3,1)': Solid_3DP1,
        'solid(3,2)': Solid_3DP2,
        'solid(3,3)': Solid_3DP3,
    }
    has_csm = True
except ImportError as e:
    print('Import feelpp.toolboxes.solid failed: Feel++ Toolbox Solid is not available')
    pass  # module doesn't exist, deal with it.

def solid( dim=2, orderDisp=1, worldComm=None):
    """create a solid toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderDisp -- the polynomial order for the displacement (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    if not has_csm:
        raise Exception('Solid toolbox is not enabled in Feel++')
    if worldComm is None:
        worldComm=feelpp.Environment.worldCommPtr()
    key='solid('+str(dim)+','+str(orderDisp)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _csms:
        raise RuntimeError('Solid solver '+key+' not existing')
    return _csms[key]( "solid", "solid", worldComm )
