import sys
import pytest
import feelpp.core as fppc
from feelpp.toolboxes.core import *
from feelpp.toolboxes.thermoelectric import *

thermoelectric_cases = [('thermoelectric/ElectroMagnets/HL-31_H1/HL-31_H1.cfg', 3,1), 
                        ('thermoelectric/ElectroMagnets/HL-31_H1/HL-31_H1.cfg', 3,2)]


@pytest.mark.parametrize("casefile,dim,order", thermoelectric_cases)
def test_thermoelectric(casefile,dim,order):
    fppc.Environment.setConfigFile(casefile)
    f = thermoelectric(dim=dim, orderPotential=order)
    simulate(f)

    return not f.checkResults()
