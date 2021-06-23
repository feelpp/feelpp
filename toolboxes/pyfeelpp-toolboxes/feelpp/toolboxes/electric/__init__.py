from ._electric import *

_electrics={
    'electric(2,1)':Electric_2DP1,
    'electric(2,2)':Electric_2DP2,
    'electric(3,1)':Electric_3DP1,
    'electric(3,2)':Electric_3DP2,
}

def electric( dim=2, orderPotential=1, worldComm=core.Environment.worldCommPtr() ):
    """create a electric toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    key='electric('+str(dim)+','+str(orderPotential)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _electrics:
        raise RuntimeError('Electric solver '+key+' not existing')
    return _electrics[key]( "electric", "electric", worldComm )
