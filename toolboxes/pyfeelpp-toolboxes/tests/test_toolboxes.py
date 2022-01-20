import sys
import pytest
import feelpp
import feelpp.toolboxes as tb
from feelpp.toolboxes.cfpdes import *
from feelpp.toolboxes.electric import *
from feelpp.toolboxes.fluid import *
from feelpp.toolboxes.hdg import *
from feelpp.toolboxes.solid import *


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

