#!/usr/bin/python

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


