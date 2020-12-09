import feelpp
from ._mor import *
from ._toolboxmor import *

_toolboxmor={
    'toolboxmor(2)':ToolboxMor_2D,
    'toolboxmor(3)':ToolboxMor_3D,
}

def toolboxmor( dim=2, prefix="" ):
    """create a model toolbox mor
    Keyword arguments:
    dim -- the dimension (default: 2)
    prefix -- the prefix of the model
    """
    key='toolboxmor('+str(dim)+')'
    if key not in _toolboxmor:
        raise RuntimeError('ToolboxMor model '+key+' not existing')
    return _toolboxmor[key]( "" )
