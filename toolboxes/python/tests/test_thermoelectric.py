import sys
import pytest
import feelpp.core as fppc
from feelpp.toolboxes.core import *
from feelpp.toolboxes.thermoelectric import *
from _case_paths import toolbox_case

thermoelectric_cases = [(toolbox_case('thermoelectric/ElectroMagnets/HL-31_H1/HL-31_H1.cfg'), 3,1),
                        (toolbox_case('thermoelectric/ElectroMagnets/HL-31_H1/HL-31_H1.cfg'), 3,2)]


@pytest.mark.parametrize("casefile,dim,order", thermoelectric_cases)
def test_thermoelectric(casefile,dim,order):
    fppc.Environment.setConfigFile(casefile)
    f = thermoelectric(dim=dim, orderPotential=order)
    simulate(f)
    assert f.checkResults()
