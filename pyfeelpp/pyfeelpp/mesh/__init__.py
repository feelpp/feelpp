import pyfeelpp.core as core
from _mesh import *

_meshes={
    'mesh(3,1,3)':Mesh_3DG1R3,
    #'mesh(3,2,3)':Mesh_3DG2R3,
    #'mesh(2,1,3)':Mesh_2DG1R3,
    #'mesh(2,2,3)':Mesh_2DG2R3,
    #'mesh(1,1,3)':Mesh_1DG1R3,
    #'mesh(1,2,3)':Mesh_1DG2R3,
    'mesh(2,1,2)':Mesh_2DG1R2,
    'mesh(2,2,2)':Mesh_2DG2R2,
    #'mesh(1,1,2)':Mesh_1DG1R2,
    #'mesh(1,2,2)':Mesh_1DG2R2,
    'mesh(1,1,1)':Mesh_1DG1R1,
    #'mesh(1,2,1)':Mesh_1DG2R1
}

def mesh( dim=2, geo=1, realdim=2, worldComm=core.Environment.worldCommPtr() ):
    """create a mesh 
    The mesh is of topological dimension 'dim', in real dimension 'realdim' with geometric order 'geo'. 
    The mesh is configured with a WorldComm which provides the parallel process layout.
    Keyword arguments:
    dim -- the topological dimension (default: 2)
    geo -- the geometrical order (default: 1)
    realdim -- the real dimension (default: 2)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    """
    key='mesh('+str(dim)+','+str(geo)+','+str(realdim)+')'
    if key not in _meshes:
        raise RuntimeError('Mesh '+key+' no existing')
    return _meshes[key]( worldComm )
