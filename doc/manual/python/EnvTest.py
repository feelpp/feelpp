#!/usr/bin/python

from mpi4py import MPI
import libPyFeelpp
import sys

z=libPyFeelpp.Environment(sys.argv)
s=libPyFeelpp.Simplex()
m=libPyFeelpp.Mesh.new()
l=libPyFeelpp.loadMesh(m)
w=libPyFeelpp.Environment.worldComm()

x=libPyFeelpp.export(l);



print "Thomas"
