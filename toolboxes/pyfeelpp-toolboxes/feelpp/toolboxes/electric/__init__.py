import feelpp
from feelpp.toolboxes.core import *
has_electric = False
_electrics = None
try:
    from ._electric import *

    _electrics={
        'electric(2,1)':Electric_2DP1,
        'electric(2,2)':Electric_2DP2,
        'electric(3,1)':Electric_3DP1,
        'electric(3,2)':Electric_3DP2,
    }
    has_electric = True
except ImportError as e:
    print('Import feelpp.toolboxes.electric failed: Feel++ Toolbox electric is not available')
    pass  # module doesn't exist, deal with it.


def electric(dim=2, orderPotential=1, worldComm=None, subprefix="", modelRep=None ):
    """create a electric toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    if not has_electric:
        raise Exception('Electric toolbox is not enabled in Feel++')
    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()
    key='electric('+str(dim)+','+str(orderPotential)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _electrics:
        raise RuntimeError('Electric solver '+key+' not existing')
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _electrics[key]( "electric", "electric", worldComm, "", modelRep  )
