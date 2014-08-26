#!/usr/bin/python

##
# \file feelpp-test.py
# \author Thomas Lantz <lantz.thomas0@gmail.com>
# \date 2014-08-26



from mpi4py import MPI
import libPyFeelpp
import sys

z=libPyFeelpp.Environment(sys.argv)

s=libPyFeelpp.Simplex1()
m=libPyFeelpp.MeshS2()
l=libPyFeelpp.loadMesh(m)
Vh=libPyFeelpp.newPch2(m)
w=libPyFeelpp.Environment.worldComm()

x=libPyFeelpp.export(l);


