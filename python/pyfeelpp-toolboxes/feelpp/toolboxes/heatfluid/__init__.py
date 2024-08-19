import feelpp.core as fppc
from feelpp.toolboxes.core import *

has_hf = False
_heatfluids = None
try:
    from ._heatfluid import *

    _heatfluids={
        'heatfluid(2,1,1,1)':HeatFluid_2D_P1_P1P1,
        'heatfluid(2,1,2,1)':HeatFluid_2D_P1_P2P1,
        'heatfluid(3,1,1,1)':HeatFluid_3D_P1_P1P1,
        'heatfluid(3,1,2,1)':HeatFluid_3D_P1_P2P1,
    }
    print(_heatfluids)
    has_hf = True
except ImportError as e:
    print('Import feelpp.toolboxes.heatfluid failed: Feel++ Toolbox heatfluid is not available')
    pass  # module doesn't exist, deal with it.


def heatfluid(dim=2, orderTemperature=1, orderVelocity=1, orderPressure=1, worldComm=None, keyword="heatfluid", prefix="heat-fluid", subprefix="", modelRep=None):
    """create a heatfluid toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderTemperature -- the polynomial order for the temperature (default: 1)
    orderVelocity -- the polynomial order for the velocity (default: 1)
    orderPressure -- the polynomial order for the pressure (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    keyword -- the json keyword for the toolbox (default: "heatfluid")
    prefix -- the prefix for the toolbox for the command line and .cfg options (default: "heat-fluid")
    subprefix -- the subprefix for the toolbox for the command line and .cfg options (default: "")
    """
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key=f'heatfluid({dim},{orderTemperature},{orderVelocity},{orderPressure})'
    if worldComm.isMasterRank():
        print(f"heatfluid:: key:{key}, has_hf:{has_hf}")
    if key not in _heatfluids:
        raise RuntimeError('heatfluid solver '+key+' not existing')
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _heatfluids[key](prefix=prefix, keyword=keyword, worldComm=worldComm, subprefix="", modelRep=modelRep)
