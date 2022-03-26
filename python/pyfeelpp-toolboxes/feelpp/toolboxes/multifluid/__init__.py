from pyfeelpp import core
from pyfeelpp.toolboxes.modelcore import *
from _thermoelectric import *

_thermoelectrics={
    'thermoelectric(2,1)':Thermoelectric_2DP1,
    'thermoelectric(2,2)':Thermoelectric_2DP2,
    'thermoelectric(3,1)':Thermoelectric_3DP1,
    'thermoelectric(3,2)':Thermoelectric_3DP2,
}

def thermoelectric( dim=2, orderPotential=1, buildMesh=True, worldComm=core.Environment.worldCommPtr() ):
    """create a thermoelectric toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    key='thermoelectric('+str(dim)+','+str(orderPotential)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _thermoelectrics:
        raise RuntimeError('Thermoelectric solver '+key+' not existing')
    return _thermoelectrics[key]( "thermoelectric", buildMesh, worldComm )
