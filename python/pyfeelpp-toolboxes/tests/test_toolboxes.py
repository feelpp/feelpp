import sys
import pytest
import feelpp.core as fppc
import feelpp.toolboxes as fppt 
from feelpp.toolboxes.cfpdes import *
from feelpp.toolboxes.electric import *
from feelpp.toolboxes.fluid import *
from feelpp.toolboxes.hdg import *
from feelpp.toolboxes.solid import *
from feelpp.toolboxes.thermoelectric import *
from feelpp.toolboxes.heatfluid import *

def test_init_cfpdes():
    f = cfpdes(dim=2)

def test_init_solid():
    f = solid(dim=2)

def test_init_fluid():
    f = fluid(dim=2)

#def test_init_hdg():
#    f = hdgpoisson(dim=2)

def test_init_electric():
    f = electric(dim=2)

def test_init_thermoelectric():
    f = thermoelectric(dim=2)

def test_init_heatfluid():
    f = heatfluid(dim=2)
