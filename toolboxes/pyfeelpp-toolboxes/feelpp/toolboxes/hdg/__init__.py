import feelpp
from feelpp.toolboxes.core import *
has_hdg = False
_hdgs = None
try:
    from ._hdg import *

    _hdgs={
        'hdg(2,1)':HDGPoisson_2DP1,
        'hdg(2,2)':HDGPoisson_2DP2,
        'hdg(3,1)':HDGPoisson_3DP1,
        'hdg(3,2)':HDGPoisson_3DP2,
    }
    has_hdg = True
except ImportError as e:
    print('Import feelpp.toolboxes.hdg failed: Feel++ Toolbox HDG is not available')
    pass  # module doesn't exist, deal with it.


def hdgpoisson(dim=2, order=1, prefix="", prefix_toolbox="hdg.poisson", physic=None, worldComm=None, subprefix="", modelRep=None):
    """create a hdg toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the fields : potential, flux, displacement, stress and associated traces (default: 1)
    prefix -- application prefix for the HDG poisson
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    if not has_hdg:
        raise Exception('HDG toolbox is not enabled in Feel++')
    if worldComm is None:
        worldComm=feelpp.Environment.worldCommPtr()
    key='hdg('+str(dim)+','+str(order)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _hdgs:
        raise RuntimeError('HDG solver '+key+' not existing')
    _prefix= prefix_toolbox+"."+prefix if prefix else prefix_toolbox
    if modelRep is None:
        modelRep = ModelBaseRepository()
    if physic is None:
        physic = MixedPoissonPhysics.none
    return _hdgs[key](prefix=_prefix, physic=physic, worldComm=worldComm, subprefix="", modelRep=modelRep)
