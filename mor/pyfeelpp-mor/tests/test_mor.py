import sys
import pytest
import feelpp
import feelpp.toolboxes as tb
from feelpp.toolboxes.heat import *
from feelpp.mor import *


def test_init_heat():
    f = heat(dim=2)

