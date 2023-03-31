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
    print('Import feelpp.toolboxes.cfpdes failed: Feel++ Toolbox cfpdes is not available')
    pass  # module doesn't exist, deal with it.


def cfpdes(dim=2, worldComm=None, keyword="cfpdes", prefix="cfpdes", subprefix="", modelRep=None ):
    """create a cfpde toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    worldComm -- the parallel communicator for the mesh (default: feelpp.Environment::worldCommPtr())
    keyword -- the json keyword for the toolbox    (default: "cfpdes")
    prefix -- the prefix for the toolbox for the command line and .cfg options  (default: "cfpdes")
    subprefix -- the subprefix for the toolbox for the command line and .cfg options (default: "")
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
    pdes = _cfpdes[key](prefix=prefix, keyword=keyword, worldComm=worldComm, subprefix=subprefix, modelRep=modelRep)
    pdes.pde = lambda nameeq: pde(pdes,nameeq)
    return pdes


    
def dispatch_call(funcs, *args, **kwargs):
    """
    Dispatches a call to a list of functions and returns the result of the first
    function that does not return None.

    :param funcs: A list of functions to call.
    :param args: Arguments to pass to the functions.
    :param kwargs: Keyword arguments to pass to the functions.
    :return: The result of the first function that does not return None, or None if all
             functions return None.
    """
    for func in funcs:
        result = func(*args, **kwargs)
        if result is not None:
            return result
    return None

def pde(toolbox,nameeq):
    """get the pde from toolbox with nameeq
    Keyword arguments:
    toolbox -- the toolbox
    nameeq -- the name of the equation
    """
#    if not has_cfpde:
#        raise Exception('CFPDE toolbox is not enabled in Feel++')
    result = dispatch_call([toolbox.pdePch1, toolbox.pdePch2, toolbox.pdePchv1, toolbox.pdePchv2], nameeq)
    return result