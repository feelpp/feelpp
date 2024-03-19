import feelpp.core as fppc
import feelpp.toolboxes as fppt
from ._fsi import *

_fsis={
    'fsi(2,2,1,1)':Fsi_2DP1,
    'fsi(2,3,2,1)':Fsi_2DP2,
    'fsi(3,2,1,1)':Fsi_3DP1,
    'fsi(3,3,2,1)':Fsi_3DP2,
}

def fsi( dim=2, orderU=2, orderP=1, orderGeo=1, orderDisp=None, buildMesh=True, worldComm=None ):
    """create a fsi toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    worldComm -- the parallel communicator for the mesh (default: fppc.Environment::worldCommPtr())
    """
    if orderDisp is None:
        orderDisp=orderGeo
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key='fsi('+str(dim)+','+str(orderPotential)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _fsis:
        raise RuntimeError('Fsi solver '+key+' not existing')
    return _fsis[key]( "fsi", buildMesh, worldComm )
