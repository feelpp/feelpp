#!/usr/bin/python

##
# \file mylaplacian.py
# \author Thomas Lantz <lantz.thomas0@gmail.com>
# \date 2014-08-26



from mpi4py import MPI
import libPyLapla
from libPyFeelpp import *
import sys


env=Environment(sys.argv)
w=Environment.worldComm()

mesh=libPyLapla.unitSquare()
Vh=newPch1(mesh)

u=libPyLapla.element(Vh)
v=libPyLapla.element(Vh)

a0=libPyLapla.form2(Vh,Vh)
l0=libPyLapla.form1(Vh)

a0=libPyLapla.integrate(a0,mesh,u,v)
l0=libPyLapla.integrate(l0,mesh,v)

a0=libPyLapla.on(a0,l0,mesh,u)

us=libPyLapla.solve(a0,l0,u)

libPyLapla.exporter(mesh,us)
