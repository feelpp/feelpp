#!/usr/bin/python

##
# \file myintegrals.py
# \author Thomas Lantz <lantz.thomas0@gmail.com>
# \date 2014-08-26



from mpi4py import MPI
import libPyInteg
from libPyFeelpp import *
import sys


env=Environment(sys.argv)
w=Environment.worldComm()

m=MeshS3()
l=loadMesh(m)

p=newPch2(l)

g=libPyInteg.expr(libPyInteg.soption("functions.g"))

gradg=libPyInteg.grad(g)

elem=libPyInteg.elements(l)
belem=libPyInteg.boundaryelements(l)
bfaces=libPyInteg.boundaryfaces(l)

integra=libPyInteg.integrate(elem,g)
integra1=libPyInteg.integrate(belem,g)
integra2=libPyInteg.integrate(elem,gradg)
integra3=libPyInteg.integrate(bfaces,g)

sol=libPyInteg.evaluate(integra)
sol1=libPyInteg.evaluate(integra1)
sol2=libPyInteg.evaluate(integra2)
sol3=libPyInteg.evaluate(integra3)

libPyInteg.printExpr(g)
print ("elements:" ) 
libPyInteg.printSol(sol)

print ("boundaryfaces:" ) 
libPyInteg.printSol(sol3)

print ("boundaryelements:" ) 
libPyInteg.printSol(sol1)

libPyInteg.printExpr(gradg)
print ("elements:" ) 
libPyInteg.printSol(sol2)

x=export(l);





