#!/usr/bin/python

##
# \file myfunctionspace.py
# \author Thomas Lantz <lantz.thomas0@gmail.com>
# \date 2014-08-26

from mpi4py import MPI
import libPyFunctSpace
import sys


libPyFunctSpace.wrap(sys.argv)
