import sys,os
import pytest
import feelpp.core as fppc 
from feelpp.toolboxes.core import *
from feelpp.toolboxes.electric import *
from feelpp.toolboxes.fluid import *

def test_multipletoolbox(init_feelpp):
    # create the application
    fppc.Environment.setConfigFile("python/pyfeelpp-toolboxes/tests/electric/quarter-turn/2d.cfg")
    s = electric(dim=2)
    # # get displacement and von-mises measures from the model
    ok,meas=simulate(s)
    assert( not meas )
    
    fppc.Environment.setConfigFile("python/pyfeelpp-toolboxes/tests/fluid/TurekHron/cfd1.cfg")
    s = fluid(dim=2)
    # # get displacement and von-mises measures from the model
    ok,meas=simulate(s)

    import pandas as pd
    df=pd.DataFrame(meas)
    print(df.head())

    assert(meas )