import feelpp
from feelpp.toolboxes.core import *

has_cfpde = False
_cfpdes = None
try:
    from ._cfpdes import *

    _cfpdes={
        'cfpdes(2)':cfpdes_2D,
        'cfpdes(3)':cfpdes_3D
    }
    has_cfpde = True
except ImportError as e:
    print('Import feelpp.toolboxes.cfpde failed: Feel++ Toolbox cfpde is not available')
    pass  # module doesn't exist, deal with it.


def cfpdes(dim=2, worldComm=None, keyword="cfpdes", subprefix="", modelRep=None ):
    """create a cfpde toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    worldComm -- the parallel communicator for the mesh (default: feelpp.Environment::worldCommPtr())
    """
    if not has_cfpde:
        raise Exception('CFPDE toolbox is not enabled in Feel++')
    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()
    key='cfpdes('+str(dim)+')'
    if key not in _cfpdes:
        raise RuntimeError('cfpde solver '+key+' not existing')
    if modelRep is None:
        modelRep = ModelBaseRepository()
    return _cfpdes[key](prefix="cfpdes", keyword=keyword, worldComm=worldComm, subprefix="", modelRep=modelRep)



    
