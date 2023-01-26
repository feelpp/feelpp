import feelpp
from feelpp.toolboxes.core import *

has_te = False
_thermoelectrics = None
try:
    from ._thermoelectric import *

    _thermoelectrics={
        'thermoelectric(2,1)':Thermoelectric_2DP1,
        'thermoelectric(2,2)':Thermoelectric_2DP2,
        'thermoelectric(3,1)':Thermoelectric_3DP1,
        'thermoelectric(3,2)':Thermoelectric_3DP2,
    }
    has_te = True
except ImportError as e:
    print('Import feelpp.toolboxes.thermoelectric failed: Feel++ Toolbox ThermoElectric is not available')
    pass  # module doesn't exist, deal with it.


def thermoelectric(dim=2, orderPotential=1, worldComm=None, keyword="thermoelectric", subprefix="", modelRep=None):
    """create a thermoelectric toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()
    key='thermoelectric('+str(dim)+','+str(orderPotential)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _thermoelectrics:
        raise RuntimeError('Thermoelectric solver '+key+' not existing')
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _thermoelectrics[key](prefix="thermo-electric", keyword=keyword, worldComm=worldComm, subprefix="", modelRep=modelRep)
