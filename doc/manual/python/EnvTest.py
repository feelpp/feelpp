#!/usr/bin/python

from mpi4py import MPI
#import libEnvironment
#import libMesh
import libPyFeelpp
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

#x=libEnvironment.Environment(sys.argv)
#k=libMesh.Simplex()
#print k.dim()

#m=libMesh.Mesh()
#m=libMesh.Mesh.new()

#s=libMesh.loadMesh(m);
#B=libMesh.loadMesh(libMesh.Mesh())

#m.clear()

z=libPyFeelpp.Environment(sys.argv)
#g=libPyFeelpp.Environment.worldComm()
s=libPyFeelpp.Simplex()
m=libPyFeelpp.Mesh.new()
l=libPyFeelpp.loadMesh(m)
#g=libPyFeelpp.exporter(s)
#g=libPyFeelpp.new()
w=libPyFeelpp.Environment.worldComm()
#e=libPyFeelpp.new();
#e=libPyFeelpp.Exporter(w)
#e.setMesh(s)

#libPyFeelpp.setMesh1(e,s)

#e.addRegions()
#e.save()


x=libPyFeelpp.export(l);


#d=libPyFeelpp.export2(m,w)


print "Thomas"
