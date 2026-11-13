import feelpp.core as fppc
from feelpp.toolboxes.core import *

has_te = False
_thermoelectrics_dynamic = None
try:
    from ._thermoelectric import *

    _thermoelectrics_dynamic = {2: Thermoelectric_2DDynamic, 3: Thermoelectric_3DDynamic}
    has_te = True
except ImportError as e:
    print('Import feelpp.toolboxes.thermoelectric failed: Feel++ Toolbox ThermoElectric is not available')
    pass  # module doesn't exist, deal with it.


def thermoelectric(dim=2, orderPotential=1, worldComm=None, keyword="thermoelectric", prefix="thermo-electric", subprefix="", modelRep=None, orderTemperature=None, orderGeometry=1):
    """create a thermoelectric toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderPotential -- the polynomial order for the potential (default: 1)
    orderTemperature -- the temperature order (default: orderPotential)
    orderGeometry -- the shared geometry order (default: 1)
    worldComm -- the parallel communicator for the mesh (default: core.Environment::worldCommPtr())
    keyword -- the json keyword for the toolbox    (default: "thermoelectric")
    prefix -- the prefix for the toolbox for the command line and .cfg options  (default: "thermo-electric")
    subprefix -- the subprefix for the toolbox for the command line and .cfg options (default: "")
    """
    if not has_te:
        raise Exception('Thermoelectric toolbox is not enabled in Feel++')
    if orderTemperature is None:
        orderTemperature = orderPotential
    for value in (orderTemperature, orderPotential, orderGeometry):
        if not isinstance(value, int) or not 1 <= value <= 65535:
            raise ValueError("Thermoelectric orders must be integers between 1 and 65535")
    if worldComm is None:
        worldComm = fppc.Environment.worldCommPtr()
    key='thermoelectric('+str(dim)+','+str(orderPotential)+')'
    if modelRep is None:
        modelRep = ModelBaseRepository()
    if dim not in _thermoelectrics_dynamic:
        raise RuntimeError('Thermoelectric solver '+key+' not existing')
    return _thermoelectrics_dynamic[dim](prefix=prefix, keyword=keyword, worldComm=worldComm,
                                        subprefix=subprefix, modelRep=modelRep, orderTemperature=orderTemperature,
                                        orderPotential=orderPotential, orderGeometry=orderGeometry)
