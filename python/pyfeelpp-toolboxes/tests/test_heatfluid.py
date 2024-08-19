import sys
import pytest
import feelpp.core as fppc
from feelpp.toolboxes.core import *
from feelpp.toolboxes.heatfluid import *

heatfluid_cases = [('heatfluid/NaturalConvection/cavity/2d_laminar.cfg', 2,1,1,1), ('heatfluid/NaturalConvection/cavity/2d_laminar.cfg', 2,1,2,1)]


@pytest.mark.parametrize("casefile,dim,orderT,orderV,orderP", heatfluid_cases)
def test_heatfluid(casefile,dim,orderT,orderV,orderP):
    fppc.Environment.setConfigFile(casefile)
    f = heatfluid(dim=dim, orderTemperature=orderT,orderVelocity=orderV,orderPressure=orderP)
    simulate(f)

    return not f.checkResults()
