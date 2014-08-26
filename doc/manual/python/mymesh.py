#!/usr/bin/python

##
# \file mymesh.py
# \author Thomas Lantz <lantz.thomas0@gmail.com>
# \date 2014-08-26



from mpi4py import MPI
from libPyFeelpp import *
import sys

z=Environment(sys.argv)
m=MeshS2()
l=loadMesh(m)
w=Environment.worldComm()

x=export(l);




