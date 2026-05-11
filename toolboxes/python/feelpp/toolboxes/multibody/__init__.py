import feelpp.core as fppc
from feelpp.toolboxes.core import *

has_multibody = False
_multibodies = None
try:
    from ._multibody import *

    _multibodies = {
        "multibody(2,1)": Multibody_2DG1,
        "multibody(2,2)": Multibody_2DG2,
        "multibody(3,1)": Multibody_3DG1,
        "multibody(3,2)": Multibody_3DG2,
    }
    has_multibody = True
except ImportError:
    print("Import feelpp.toolboxes.multibody failed: Feel++ Toolbox Multibody is not available")
    pass


def multibody(dim=2, orderGeo=1, worldComm=None, keyword="multibody", prefix="multibody", modelRep=None):
    """create a multibody toolbox solver

    Keyword arguments:
    dim -- the dimension (default: 2)
    orderGeo -- the geometric order (default: 1)
    worldComm -- the parallel communicator (default: fppc.Environment.worldCommPtr())
    keyword -- the json keyword for the toolbox (default: multibody)
    prefix -- the prefix for toolbox options (default: multibody)
    """
    if not has_multibody:
        raise Exception("Multibody toolbox is not enabled in Feel++")
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key = "multibody(" + str(dim) + "," + str(orderGeo) + ")"
    if key not in _multibodies:
        raise RuntimeError("Multibody solver " + key + " not existing")
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _multibodies[key](prefix=prefix, keyword=keyword, worldComm=worldComm, modelRep=modelRep)
