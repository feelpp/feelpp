import sys
import pytest
import feelpp.core as fppc
from feelpp.toolboxes.core import *
from feelpp.toolboxes.heatfluid import *
from _case_paths import toolbox_case

heatfluid_cases = [(toolbox_case('heatfluid/NaturalConvection/cavity/2d_laminar.cfg'), 2,1,1,1), (toolbox_case('heatfluid/NaturalConvection/cavity/2d_laminar.cfg'), 2,1,2,1)]


@pytest.mark.parametrize("casefile,dim,orderT,orderV,orderP", heatfluid_cases)
def test_heatfluid(casefile,dim,orderT,orderV,orderP):
    fppc.Environment.setConfigFile(casefile)
    f = heatfluid(dim=dim, orderTemperature=orderT,orderVelocity=orderV,orderPressure=orderP)
    simulate(f)
    assert f.checkResults()
