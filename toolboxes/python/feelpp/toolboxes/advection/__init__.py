import feelpp.core as fppc
from feelpp.toolboxes.core import *

has_advection = False
_advections = None
try:
    from ._advection import *

    _advections = {
        "advection(2,1)": Advection_2DP1,
        "advection(2,2)": Advection_2DP2,
        "advection(3,1)": Advection_3DP1,
        "advection(3,2)": Advection_3DP2,
    }
    has_advection = True
except ImportError:
    print("Import feelpp.toolboxes.advection failed: Feel++ Toolbox Advection is not available")
    pass


def advection(dim=2, order=1, buildMesh=True, worldComm=None, prefix="advection", subprefix="", modelRep=None):
    """create an advection toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: fppc.Environment::worldCommPtr())
    """
    if not has_advection:
        raise Exception("Advection toolbox is not enabled in Feel++")
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key = "advection(" + str(dim) + "," + str(order) + ")"
    if worldComm.isMasterRank():
        print(key)
    if key not in _advections:
        raise RuntimeError("Advection solver " + key + " not existing")
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _advections[key](prefix=prefix, buildmesh=buildMesh, worldComm=worldComm, subprefix=subprefix, modelRep=modelRep)
