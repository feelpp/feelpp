import feelpp
from feelpp.toolboxes.core import *

has_te = False
_heatfluids = None
try:
    from _heatfluid import *

    _heatfluids={
        'heatfluid(2,1)':HeatFluid_2D_P1_P1P1,
        'heatfluid(2,2)':HeatFluid_2D_P1_P2P1,
        'heatfluid(3,1)':HeatFluid_3D_P1_P1P1,
        'heatfluid(3,2)':HeatFluid_3D_P1_P2P1,
    }
    has_te = True
except ImportError as e:
    print('Import feelpp.toolboxes.heatfluid failed: Feel++ Toolbox heatfluid is not available')
    pass  # module doesn't exist, deal with it.


def heatfluid(dim=2, orderPotential=1, worldComm=None, keyword="thermo-electric", subprefix="", modelRep=None):
    """create a heatfluid toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()
    key='heatfluid('+str(dim)+','+str(orderPotential)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _heatfluids:
        raise RuntimeError('heatfluid solver '+key+' not existing')
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _heatfluids[key](prefix="heatfluid", keyword=keyword, worldComm=worldComm, subprefix="", modelRep=modelRep)
