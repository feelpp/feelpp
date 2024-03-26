import math
import sys
import feelpp
from feelpp.operators import *
import pytest


def fem(m, mesh_name, h=0.1):
    # load the mesh
    mesh = feelpp.load(m, mesh_name, h)
    # function space
    Xh = feelpp.functionSpace(mesh=mesh)

    # a function
    v = Xh.element()
    v.on(range=feelpp.elements(mesh), expr=feelpp.expr("1"))

    k = []
    Bi=0.1
    Aq = []
    for q in [0,1,2,3,4,5]:
        if q <= 4:
            k.append(1.0+q)
            Aq.append( stiffness(test=Xh, trial=Xh, range=feelpp.markedelements(mesh,""), coeff=str(k[-1])) )
        else:
            Aq.append( mass(test=Xh, trial=Xh, range=feelpp.markedfaces(
                            mesh, "Gamma_ext"), coeff=str(Bi)) )
    E = Aq[-1].energy(v, v)/Bi
    print(f"perimeter fin\Root:{E})")
    print("Done!")

feelppenv = feelpp.Environment(['fin'],config=feelpp.localRepository('.'))
fem(feelpp.mesh(dim=2), "/nvme0/prudhomm/Devel/feelpp/python/pyfeelpp/tests/fin.geo")




