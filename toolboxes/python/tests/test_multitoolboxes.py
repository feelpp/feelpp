import sys,os
import pytest
import feelpp.core as fppc 
from feelpp.toolboxes.core import *
from feelpp.toolboxes.electric import *
from feelpp.toolboxes.fluid import *
from _case_paths import toolbox_case

def test_multipletoolbox(init_feelpp):
    # create the application
    electric_cfg = toolbox_case("electric/quarter-turn/2d.cfg")
    fppc.Environment.setConfigFile(electric_cfg)
    s = electric(dim=2)
    # # get displacement and von-mises measures from the model
    ok,meas=simulate(s)
    assert( not meas )
    
    fluid_cfg = toolbox_case("fluid/TurekHron/cfd1.cfg")
    fppc.Environment.setConfigFile(fluid_cfg)
    s = fluid(dim=2)
    # # get displacement and von-mises measures from the model
    ok,meas=simulate(s)

    import pandas as pd
    df=pd.DataFrame(meas)
    print(df.head())

    assert(meas )
