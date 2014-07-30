#!/usr/bin/python

from mpi4py import MPI
import libPyLapla
import sys


env=libPyLapla.Environment(sys.argv)
#w=libPyLapla.Environment.worldComm()

mesh=libPyLapla.unitSquare()
Vh=libPyLapla.newPch(mesh,False)

u=libPyLapla.element(Vh)
v=libPyLapla.element(Vh)

a0=libPyLapla.form2(Vh,Vh)
l0=libPyLapla.form1(Vh)

a0=libPyLapla.integrate(a0,mesh,u,v)
l0=libPyLapla.integrate(l0,mesh,v)

a0=libPyLapla.on(a0,l0,mesh,u)

us=libPyLapla.solve(a0,l0,u)

libPyLapla.exporter(mesh,us)
