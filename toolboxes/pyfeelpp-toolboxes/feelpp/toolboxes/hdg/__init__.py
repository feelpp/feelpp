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
    print('Import feelpp.toolboxes.hdg failed: Feel++ Toolbox HDG is not available')
    pass  # module doesn't exist, deal with it.

def mixedpoisson(dim=2, order=1, prefix="", prefix_toolbox="hdg.poisson", physic=None, worldComm=None, subprefix="", modelRep=None):
    """create a hdg toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the fields : potential, flux, displacement, stress and associated traces (default: 1)
    prefix -- application prefix for the HDG poisson
    prefix_toolbox -- toolbox prefix
    physic -- physic to use
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    subprefix -- subprefix
    modelRep -- model repository
    """
    if not has_hdg:
        raise Exception('HDG toolbox is not enabled in Feel++')
    if worldComm is None:
        worldComm=feelpp.Environment.worldCommPtr()
    key='mixedpoisson('+str(dim)+','+str(order)+')'
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
