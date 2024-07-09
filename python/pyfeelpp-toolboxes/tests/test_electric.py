import sys
import pytest
import feelpp.core as fppc
from feelpp.toolboxes.core import *
from feelpp.toolboxes.electric import *

electric_cases = [('electric/quarter-turn/2d.cfg', 2,1), ('electric/quarter-turn/2d.cfg', 2,2),
                  ('electric/quarter-turn/3d.cfg', 3,1),('electric/quarter-turn/3d.cfg', 3,2),
                  #('electric/busbar/2d.cfg', 2, 1), ('electric/busbar/2d.cfg', 2, 2),
                  ('electric/busbar/3d.cfg', 3,1), ('electric/busbar/3d.cfg', 3,2)]


@pytest.mark.parametrize("casefile,dim,order", electric_cases)
def test_electric(casefile,dim,order):
    fppc.Environment.setConfigFile(casefile)
    f = electric(dim=dim, orderPotential=order)
    simulate(f)

    return not f.checkResults()
