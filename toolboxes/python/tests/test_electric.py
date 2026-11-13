import sys
from pathlib import Path
import pytest
import feelpp.core as fppc
from feelpp.toolboxes.core import *
from feelpp.toolboxes.electric import *
from _case_paths import toolbox_case

electric_cases = [
    (toolbox_case('electric/quarter-turn/2d.cfg'), 2, order, geo)
    for order, geo in [(1, 1), (1, 3), (2, 1), (2, 3), (3, 2), (4, 1)]
]


@pytest.mark.parametrize("casefile,dim,order,geo_order", electric_cases)
def test_electric(casefile, dim, order, geo_order, case_repository):
    fppc.Environment.setConfigFile(casefile)
    fppc.Environment.changeRepository(directory=case_repository)
    f = electric(dim=dim, orderPotential=order, orderGeometry=geo_order)
    f.setModelProperties(str(Path(casefile).with_suffix(".json")))
    simulate(f,export=True)
    assert f.spaceElectricPotential().order() == order
    assert f.mesh().order() == geo_order
    assert f.checkResults()
