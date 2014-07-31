## to move into the folders where Python libraries are created

#!/usr/bin/python

from mpi4py import MPI
import libPyMesh
import sys

z=libPyMesh.Environment(sys.argv)
s=libPyMesh.Simplex()
m=libPyMesh.Mesh.new()
l=libPyMesh.loadMesh(m)
w=libPyMesh.Environment.worldComm()

x=libPyMesh.export(l);




