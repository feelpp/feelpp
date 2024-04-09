import sys,os
import feelpp.core as fppc
import pytest
from pathlib import Path
from feelpp.toolboxes.core import *
from feelpp.toolboxes.hdg import *

heat_cases = [
    ('hdg/poisson/convection-diffusion/convection-diffusion-2d-square.cfg', 2, 1)]


@pytest.mark.parametrize("casefile,dim,order", heat_cases)
def test_heat(casefile,dim,order):
    fppc.Environment.setConfigFile(casefile)
    f = mixedpoisson(dim=dim, order=order)
    if not f.isStationary():
        f.setTimeFinal(10*f.timeStep())
    simulate(f)
    meas = f.postProcessMeasures().values()

    try:
        import pandas as pd

        df=pd.DataFrame([meas])
    except ImportError:
        print("cannot import pandas, no problem it was just a test")

    return not f.checkResults()
