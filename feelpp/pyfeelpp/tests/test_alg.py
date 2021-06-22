import feelpp
import sys
import pytest
from petsc4py import PETSc

def test_alg():
    wc=feelpp.Environment.worldCommPtr()
    v = feelpp.VectorPetscDouble(10,wc)
    if feelpp.Environment.isMasterRank():
        print("vector size:",v.size())
    assert(v.size()==10)

    
    
    dm = feelpp.DataMap(wc)
    if feelpp.Environment.isMasterRank():
        print("dm size:", dm.nDof())
    assert(dm.nDof()==0)

    w=v.vec()

    from math import sqrt
    w.set(1)
    n1 = w.norm(PETSc.NormType.NORM_1)
    n2 = w.norm(PETSc.NormType.NORM_2)
    ni = w.norm(PETSc.NormType.NORM_INFINITY)

    assert(n1 == w.getSize())
    assert(n2 == sqrt(w.getSize()))

    # scalar product with itself
    d=w.dot(w)
    assert(abs(d) == w.getSize())


    return
    dmrow=feelpp.DataMap(10, 10, wc)
    if feelpp.Environment.isMasterRank():
        print("dm row size:", dmrow.nDof())
        print("dm row local size:", dmrow.nLocalDof())
    assert(dmrow.nDof()==10)
    assert(dmrow.nLocalDof() == 10)
    M = feelpp.MatrixPetscDouble(dmrow, dmrow, wc)
#    return
    if feelpp.Environment.isMasterRank():
        print("matrix rows:{}, cols:{}".format(M.size1(),M.size2()))
    
#    assert(M.size1() == 10)
