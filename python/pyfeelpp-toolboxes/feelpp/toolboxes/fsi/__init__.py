import feelpp.core as fppc
from feelpp.toolboxes.core import *

has_fsi = False
_fsis=None
try :
    from ._fsi import *

    _fsis={
        'fsi(2,2,1,1)':Fsi_2DP2P1G1,
        'fsi(2,3,2,1)':Fsi_2DP3P2G1,
        'fsi(3,2,1,1)':Fsi_3DP2P1G1,
        'fsi(3,3,2,1)':Fsi_3DP3P2G1,
    }
    has_fsi =True
except ImportError as e:
    print('Import feelpp.toolboxes.fsi failed: Feel++ Toolbox Fsi is not avaible')
    pass #module doesn't exist, deal with it

def fsi( dim=2, orderU=2, orderP=1, orderGeo=1, orderDisp=None, buildMesh=True, worldComm=None, modelRep = None ):
    """create a fsi toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: fppc.Environment::worldCommPtr())
    """
    if not has_fsi:
        raise Exception('Fsi toolbox is not enabled in Feel++')
    if orderDisp is None:
        orderDisp=orderGeo
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key='fsi('+str(dim)+','+str(orderU)+','+str(orderP)+','+str(orderGeo)+')'
    if worldComm.isMasterRank():
        print(f"Instantiate fsi toolbox {key}")
    if key not in _fsis:
        raise RuntimeError('Fsi solver '+key+' not existing')
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _fsis[key]( "fsi", "fsi", worldComm, modelRep)


