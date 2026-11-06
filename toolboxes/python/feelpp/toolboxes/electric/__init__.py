import feelpp.core as fppc
from feelpp.toolboxes.core import *
has_electric = False
_electrics = None
_electrics_dynamic = None
try:
    from ._electric import *

    _electrics={
        'electric(2,1)':Electric_2DP1,
        'electric(3,1)':Electric_3DP1,
    }
    _electrics_dynamic={
        2: Electric_2DDynamic,
        3: Electric_3DDynamic,
    }
    has_electric = True
except ImportError as e:
    print('Import feelpp.toolboxes.electric failed: Feel++ Toolbox electric is not available')
    pass  # module doesn't exist, deal with it.


def electric(dim=2, orderPotential=1, orderGeometry=1, worldComm=None, keyword="electric", prefix="electric", subprefix="", modelRep=None):
    """create a electric toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    orderGeometry -- the polynomial order for the geometry (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    keyword -- the json keyword for the toolbox    (default: "electric")
    prefix -- the prefix for the toolbox for the command line and .cfg options  (default: "electric")
    subprefix -- the subprefix for the toolbox for the command line and .cfg options (default: "")
    """
    if not has_electric:
        raise Exception('Electric toolbox is not enabled in Feel++')
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    if modelRep is None:
        modelRep = ModelBaseRepository()

    if orderPotential == 1 and orderGeometry == 1:
        key = 'electric(' + str(dim) + ',1)'
        if key in _electrics:
            return _electrics[key]( prefix=prefix, keyword=keyword, worldComm=worldComm, subprefix=subprefix, modelRep=modelRep )

    if orderGeometry != 1:
        raise RuntimeError('Electric solver currently supports only geometry order 1')

    if dim not in _electrics_dynamic:
        raise RuntimeError('Electric solver electric(' + str(dim) + ',' + str(orderPotential) + ')/G' + str(orderGeometry) + ' not existing')

    return _electrics_dynamic[dim](
        prefix=prefix,
        keyword=keyword,
        worldComm=worldComm,
        subprefix=subprefix,
        modelRep=modelRep,
        orderPotential=orderPotential,
        orderGeometry=orderGeometry,
    )
