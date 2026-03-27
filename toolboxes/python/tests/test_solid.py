import feelpp.core as fppc
import pytest
from feelpp.toolboxes.core import *
from feelpp.toolboxes.solid import *
from _case_paths import toolbox_case

solid_cases = [
    (toolbox_case("solid/NAFEMS-LE1/le1.cfg"), 2, 2),
]


@pytest.mark.parametrize("casefile,dim,order_disp", solid_cases)
def test_solid(casefile, dim, order_disp):
    fppc.Environment.setConfigFile(casefile)
    f = solid(dim=dim, orderDisp=order_disp)
    simulate(f, export=False)
    assert f.checkResults()
