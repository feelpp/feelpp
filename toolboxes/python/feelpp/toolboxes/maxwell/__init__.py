import feelpp.core as fppc
from feelpp.toolboxes.core import *

has_maxwell = False
_maxwells = None
try:
    from ._maxwell import *

    _maxwells = {
        "maxwell(2,1)": Maxwell_2DP1,
        "maxwell(3,1)": Maxwell_3DP1,
    }
    has_maxwell = True
except ImportError:
    print("Import feelpp.toolboxes.maxwell failed: Feel++ Toolbox Maxwell is not available")
    pass


def maxwell(dim=2, order=1, buildMesh=True, worldComm=None, prefix="maxwell", subprefix="", modelRep=None):
    """create a maxwell toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    order -- the polynomial order for the potential vector (default: 1)
    worldComm -- the parallel communicator for the mesh (default: fppc.Environment::worldCommPtr())
    """
    if not has_maxwell:
        raise Exception("Maxwell toolbox is not enabled in Feel++")
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key = "maxwell(" + str(dim) + "," + str(order) + ")"
    if worldComm.isMasterRank():
        print(key)
    if key not in _maxwells:
        raise RuntimeError("Maxwell solver " + key + " not existing")
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _maxwells[key](prefix=prefix, buildmesh=buildMesh, worldComm=worldComm, subprefix=subprefix, modelRep=modelRep)
