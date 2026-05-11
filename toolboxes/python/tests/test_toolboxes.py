def test_init_cfpdes():
    from feelpp.toolboxes.cfpdes import cfpdes

    f = cfpdes(dim=2)

def test_init_solid():
    from feelpp.toolboxes.solid import solid

    f = solid(dim=2)

def test_init_hdg():
    from feelpp.toolboxes.hdg import mixedpoisson

    f = mixedpoisson(dim=2)

def test_init_fluid():
    from feelpp.toolboxes.fluid import fluid

    f = fluid(dim=2)

def test_init_fsi():
    from feelpp.toolboxes.fsi import fsi

    f = fsi(dim=2)

def test_init_electric():
    from feelpp.toolboxes.electric import electric

    f = electric(dim=2)

def test_init_thermoelectric():
    from feelpp.toolboxes.thermoelectric import thermoelectric

    f = thermoelectric(dim=2)

def test_init_heatfluid():
    from feelpp.toolboxes.heatfluid import heatfluid

    f = heatfluid(dim=2)
