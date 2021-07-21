import sys,time
from ._feelpp import  *
from ._core import *
from ._alg import *
from ._mesh import *
from ._discr import *
from ._exporter import  *
from ._ts import  *
from ._vf import  *
from ._measure import  *
from ._models import *


def download(data,worldComm):
    """Download remote data file"""
    rd = RemoteData(data,worldComm)
    if rd.canDownload():
        d=Environment.downloadsRepository()
        return rd.download( d )
    else:
        raise RuntimeError("Remote Data " + data + " cannot be downloaded")

_meshes={
    'mesh(3,1,3)':Mesh_S3DG1R3,
    #'mesh(3,2,3)':Mesh_3DG2R3,
    #'mesh(2,1,3)':Mesh_2DG1R3,
    #'mesh(2,2,3)':Mesh_2DG2R3,
    #'mesh(1,1,3)':Mesh_1DG1R3,
    #'mesh(1,2,3)':Mesh_1DG2R3,
    'mesh(2,1,2)':Mesh_S2DG1R2,
    'meshstructured(2,1,2)':Mesh_H2DG1R2,
    'mesh(2,2,2)':Mesh_S2DG2R2,
    #'mesh(1,1,2)':Mesh_1DG1R2,
    #'mesh(1,2,2)':Mesh_1DG2R2,
    'mesh(1,1,1)':Mesh_S1DG1R1,
    #'mesh(1,2,1)':Mesh_1DG2R1
}

def mesh( dim=2, geo=1, realdim=2, worldComm=None ):
    """create a mesh 
    The mesh is of topological dimension 'dim', in real dimension 'realdim' with geometric order 'geo'.
    The mesh is configured with a WorldComm which provides the parallel process layout.
    Keyword arguments:
    dim -- the topological dimension (default: 2)
    geo -- the geometrical order (default: 1)
    realdim -- the real dimension (default: 2)
    worldComm -- the parallel communicator for the mesh (default: Environment::worldCommPtr())
    """
    if worldComm is None:
        worldComm=Environment.worldCommPtr()
    key='mesh('+str(dim)+','+str(geo)+','+str(realdim)+')'
    if key not in _meshes:
        raise RuntimeError('Mesh '+key+' no existing')
    return _meshes[key]( worldComm )

def meshstructured( dim=2, geo=1, realdim=2, worldComm=None ):
    """create a structured mesh 
    The mesh is of topological dimension 'dim', in real dimension 'realdim' with geometric order 'geo'.
    The mesh is configured with a WorldComm which provides the parallel process layout.
    Keyword arguments:
    dim -- the topological dimension (default: 2)
    geo -- the geometrical order (default: 1)
    realdim -- the real dimension (default: 2)
    worldComm -- the parallel communicator for the mesh (default: Environment::worldCommPtr())
    """
    if worldComm is None:
        worldComm=Environment.worldCommPtr()
    key='meshstructured('+str(dim)+','+str(geo)+','+str(realdim)+')'
    if key not in _meshes:
        raise RuntimeError('Mesh '+key+' no existing')
    return _meshes[key]( worldComm )



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
        worldscomm=Environment.worldsComm(1)
    key=space+'('+str(mesh.dimension())+','+str(order)+')'
    if key not in _spaces:
        raise RuntimeError('FunctionSpace '+key+' not existing in dictionary')
    return _spaces[key]( mesh=mesh, worldsComm=worldscomm )


def expr(e,filename=None,worldComm=None,directory=None):
    if filename is None:
        filename=""
    if directory is None:
        directory=""
    if worldComm is None:
        worldComm=Environment.worldComm()
    return expr_(e,filename=filename,dir=directory,worldComm=worldComm)
