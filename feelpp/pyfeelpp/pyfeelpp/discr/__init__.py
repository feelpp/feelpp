import pyfeelpp.core as core
import pyfeelpp.mesh as mesh
from _discr import *


_spaces={
    # Pch
    'Pch(1,1)':Pch_1D_P1,
    'Pch(1,2)':Pch_1D_P2,
    'Pch(1,3)':Pch_1D_P3,

    'Pch(2,1)':Pch_2D_P1,
    'Pch(2,2)':Pch_2D_P2,
    'Pch(2,3)':Pch_2D_P3,
    
    'Pch(3,1)':Pch_3D_P1,
    'Pch(3,2)':Pch_3D_P2,
    'Pch(3,3)':Pch_3D_P3,

    # Pdh
    #'Pdh(1,0)':Pdh_1D_P0,
    #'Pdh(1,1)':Pdh_1D_P1,
    #'Pdh(1,2)':Pdh_1D_P2,
    #'Pdh(1,3)':Pdh_1D_P3,

    'Pdh(2,0)':Pdh_2D_P0,
    'Pdh(2,1)':Pdh_2D_P1,
    'Pdh(2,2)':Pdh_2D_P2,
    'Pdh(2,3)':Pdh_2D_P3,

    'Pdh(3,0)':Pdh_3D_P0,
    'Pdh(3,1)':Pdh_3D_P1,
    'Pdh(3,2)':Pdh_3D_P2,
    'Pdh(3,3)':Pdh_3D_P3,
}

def functionSpace( mesh, space="Pch", order=1, worldscomm=None):
    """create a function space
    """
    if worldscomm is None:
        worldscomm=core.Environment.worldsComm(1)
    key=space+'('+str(mesh.dimension())+','+str(order)+')'
    if key not in _spaces:
        raise RuntimeError('FunctionSpace '+key+' not existing in dictionary')
    return _spaces[key]( mesh=mesh, worldsComm=worldscomm )


