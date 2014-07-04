# to place where the module are created

#!/usr/bin/python

from mpi4py import MPI
import libEnvironment
import libMesh
import sys

#from libEnvironment import *

#class Test(Argv):
#    def __init__(self):
#        Argv.__init__(self)

#n=create(Test)

#x=libEnvironment.Argv.create()
#x=libEnvironment.Argv();
#x=libEnvironment.Argv(sys.argv);
#print x.test()

#libEnvironment.Environment(len(sys.argv),sys.argv)

x=libEnvironment.Environment(sys.argv)
k=libMesh.Simplex()
print k.dim()

m=libMesh.Mesh();

#k=libMesh.ConstruSimplex()
print "Thomas"
#libEnvironment.Environment()
