## to move into the folders where Python libraries are created

#!/usr/bin/python

from mpi4py import MPI
from libPyFeelpp import *
import sys

z=Environment(sys.argv)
m=MeshS2
l=loadMesh(m)
w=Environment.worldComm()

x=export(l);




