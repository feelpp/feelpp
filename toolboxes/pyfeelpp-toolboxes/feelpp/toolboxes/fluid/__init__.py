import feelpp

has_fluid = False
_cfds=None
try:
    from ._fluid import *

    _cfds={
        'fluid(2,2,1,1)':Fluid_2DP2P1G1,
        'fluid(2,3,2,1)':Fluid_2DP3P2G1,
        'fluid(3,2,1,1)':Fluid_3DP2P1G1,
        'fluid(3,3,2,1)':Fluid_3DP3P2G1,
    }
    has_fluid=True
except ImportError as e:
    print('has_fluid', has_fluid)
    pass  # module doesn't exist, deal with it.

def fluid( dim=2, orderVelocity=2, orderPressure=1, orderGeometry=1, worldComm=None ):
    """create a fluid toolbox solver
    Keyword arguments:
    dim -- the dimension (default: 2)
    orderVelocity -- the polynomial order for the velocity (default: 2)
    orderPressure -- the polynomial order for the pressure (default: 1)
    orderGeometry -- the polynomial order for the geometry (default: 1)
    worldComm -- the parallel communicator for the mesh (default: feelpp.Environment::worldCommPtr())
    """
    if not has_fluid:
        raise Exception('Fluid toolbox is not enabled in Feel++')
    if worldComm is None:
        worldComm=feelpp.Environment.worldCommPtr()
    key='fluid('+str(dim)+','+str(orderVelocity)+','+str(orderPressure)+','+str(orderGeometry)+')'
    if worldComm.isMasterRank():
        print(key)
    if key not in _cfds:
        raise RuntimeError('Fluid solver '+key+' not existing')
    return _cfds[key]( "fluid", "fluid", worldComm )
