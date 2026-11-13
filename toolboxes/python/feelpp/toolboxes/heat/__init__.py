import feelpp.core as fppc

from feelpp.toolboxes.core import *

has_heat = False
_heats_dynamic = None
try:
    from ._heat import *

    _heats_dynamic = {2: Heat_2DDynamic, 3: Heat_3DDynamic}
    has_heat = True
except ImportError as e:
    print('Import feelpp.toolboxes.heat failed: Feel++ Toolbox Heat is not available')
    pass  # module doesn't exist, deal with it.


def heat(dim=2, order=1, worldComm=None, keyword="heat", prefix="heat", subprefix="", modelRep=None,verbose=0,orderGeometry=1):
    """create a heat toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the temperature (default: 1)
    orderGeometry -- the geometry order (default: 1)
    worldComm -- the parallel communicator for the mesh (default: fppc.Environment::worldCommPtr())
    keyword -- the json keyword for the toolbox (default: heat)
    prefix -- the prefix for the toolbox command line and .cfg options (default: heat)
    subprefix -- the subprefix for the toolbox command line and .cfg options (default: "")
    """
    if not has_heat:
        raise Exception('Heat toolbox is not enabled in Feel++')
    for value in (order, orderGeometry):
        if not isinstance(value, int) or not 1 <= value <= 65535:
            raise ValueError("Heat orders must be integers between 1 and 65535")
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key='heat('+str(dim)+','+str(order)+')'
    if worldComm.isMasterRank() and verbose > 0:
        print(f"Toolbox {key} created.")
    if modelRep is None:
        modelRep = ModelBaseRepository()
    if dim not in _heats_dynamic:
        raise RuntimeError('Heat solver '+key+' not existing')
    return _heats_dynamic[dim](prefix=prefix, keyword=keyword, worldComm=worldComm,
                               subprefix=subprefix, modelRep=modelRep, order=order, orderGeometry=orderGeometry)
