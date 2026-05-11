import feelpp.core as fppc
from feelpp.toolboxes.core import *

has_fsi = False
_fsis = None
try:
    from ._fsi import *

    _fsis = {
        "fsi(2,2,1,1)": Fsi_2DP2P1G1,
        "fsi(2,3,2,1)": Fsi_2DP3P2G1,
        "fsi(3,2,1,1)": Fsi_3DP2P1G1,
        "fsi(3,3,2,1)": Fsi_3DP3P2G1,
    }
    has_fsi = True
except ImportError:
    print("Import feelpp.toolboxes.fsi failed: Feel++ Toolbox FSI is not available")
    pass


def fsi(dim=2, orderU=2, orderP=1, orderGeo=1, orderDisp=None, worldComm=None, keyword="fsi", prefix="fsi", subprefix="", modelRep=None):
    """create a fsi toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderU -- the polynomial order for the fluid velocity space (default: 2)
    orderP -- the polynomial order for the fluid pressure space (default: 1)
    orderGeo -- the geometric order used when orderDisp is not set (default: 1)
    worldComm -- the parallel communicator for the mesh (default: fppc.Environment::worldCommPtr())
    """
    if not has_fsi:
        raise Exception("FSI toolbox is not enabled in Feel++")
    if orderDisp is None:
        orderDisp = orderGeo
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key = "fsi(" + str(dim) + "," + str(orderU) + "," + str(orderP) + "," + str(orderDisp) + ")"
    if worldComm.isMasterRank():
        print(key)
    if key not in _fsis:
        raise RuntimeError("Fsi solver " + key + " not existing")
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _fsis[key](prefix=prefix, keyword=keyword, worldComm=worldComm, subprefix=subprefix, modelRep=modelRep)
