import feelpp

from feelpp.toolboxes.core import *

has_hdg = False
_hdgs = None
try:
    from ._hdgpoisson import *

    _hdgs={
        'mixedpoisson(2,1)':mixedpoisson_2DP1,
        'mixedpoisson(2,2)':mixedpoisson_2DP2,
        'mixedpoisson(3,1)':mixedpoisson_3DP1,
        'mixedpoisson(3,2)':mixedpoisson_3DP2,
    }
    has_hdg = True
except ImportError as e:
    print('has_hdg:', has_hdg)
    pass

def mixedpoisson( dim=2, order=1, prefix="", prefix_toolbox="hdg.poisson", worldComm=None ):
    """create a hdg toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the fields : potential, flux, displacement, stress and associated traces (default: 1)
    prefix -- application prefix for the HDG poisson
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    if worldComm is None:
        worldComm=feelpp.Environment.worldCommPtr()
    key='mixedpoisson('+str(dim)+','+str(order)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _hdgs:
        raise RuntimeError('HDG solver '+key+' not existing')
    _prefix= prefix_toolbox+"."+prefix if prefix else prefix_toolbox
    return _hdgs[key]( _prefix, worldComm=worldComm )
