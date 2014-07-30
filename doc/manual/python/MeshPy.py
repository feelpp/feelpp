#!/usr/bin/python

from mpi4py import MPI
#import libEnvironment
#import libMesh
import libPyMesh
import sys

z=libPyMesh.Environment(sys.argv)

s=libPyMesh.Simplex1()
m=libPyMesh.MeshS2()
l=libPyMesh.loadMesh(m)
Vh=libPyMesh.newPch2(m)
w=libPyMesh.Environment.worldComm()

x=libPyMesh.export(l);


