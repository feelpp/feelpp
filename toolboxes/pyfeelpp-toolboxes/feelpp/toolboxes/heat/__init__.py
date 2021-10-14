import feelpp

from feelpp.toolboxes.core import *

has_heat = False
_heats = None
try:
    from ._heat import *

    _heats={
        'heat(2,1)':Heat_2DP1,
        'heat(2,2)':Heat_2DP2,
        'heat(3,1)':Heat_3DP1,
        'heat(3,2)':Heat_3DP2,
    }
    has_heat = True
except ImportError as e:
    print('Import feelpp.toolboxes.heat failed: Feel++ Toolbox Heat is not available')
    pass  # module doesn't exist, deal with it.

def heat( dim=2, order=1, worldComm=None ):
    """create a heat toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the temperature (default: 1)
    worldComm -- the parallel communicator for the mesh (default: feelpp.Environment::worldCommPtr())
    """
    if not has_heat:
        raise Exception('Heat toolbox is not enabled in Feel++')
    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()
    key='heat('+str(dim)+','+str(order)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _heats:
        raise RuntimeError('Heat solver '+key+' not existing')
    return _heats[key]( "heat", "heat", worldComm )
