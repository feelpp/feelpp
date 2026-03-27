import feelpp.core as fppc
import pytest
from feelpp.toolboxes.core import *
from feelpp.toolboxes.hdg import *
from _case_paths import toolbox_case

hdg_cases = [
    (toolbox_case("hdg/poisson/convection-diffusion/convection-diffusion-2d-square.cfg"), 2, 1),
]


@pytest.mark.parametrize("casefile,dim,order", hdg_cases)
def test_hdg_poisson(casefile, dim, order):
    fppc.Environment.setConfigFile(casefile)
    f = mixedpoisson(dim=dim, order=order)
    if not f.isStationary():
        f.setTimeFinal(10 * f.timeStep())
    simulate(f)
    measures = f.postProcessMeasures().values()

    try:
        import pandas as pd

        pd.DataFrame([measures])
    except ImportError:
        print("cannot import pandas, no problem it was just a test")
    assert f.checkResults()
