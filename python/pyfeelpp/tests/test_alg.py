import feelpp
import sys
import pytest
from petsc4py import PETSc

def test_alg(init_feelpp):
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/alg/test_alg")
    wc=feelpp.Environment.worldCommPtr()
    if feelpp.Environment.isSequential():
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
 

def test_backend_vector(init_feelpp):
    e = init_feelpp
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/alg/test_vector")
    wc = feelpp.Environment.worldCommPtr()
    if feelpp.Environment.isSequential():
        b = feelpp.backend(worldcomm=wc)
        dm = feelpp.DataMap(10,10,wc)
        v=b.newVector(dm)
        assert(v.size()==10)


def test_backend_matrix(init_feelpp):
    e = init_feelpp
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/alg/test_backend_matrix")
    wc = feelpp.Environment.worldCommPtr()
    if feelpp.Environment.isSequential():
        b = feelpp.backend(worldcomm=wc)
        dm = feelpp.DataMap(10, 10, wc)
        v = b.newMatrix(dmrow=dm,dmcol=dm)
        assert(v.size1() == 10 and v.size2() == 10)
        M=v.mat()
        M.zeroEntries()

