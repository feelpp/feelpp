import sys
import pytest
import feelpp
import feelpp.toolboxes as tb
from feelpp.toolboxes.electric import *

def test_electric():
    feelpp.Environment.setConfigFile('electric/quarter-turn/2d.cfg')
    f = electric(dim=2, orderPotential=1)
    f.init()
    #f.printAndSaveInfo()
    f.solve()
    f.exportResults()
