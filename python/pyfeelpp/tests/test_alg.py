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

def test_alg_vector(init_feelpp):
    feelpp.Environment.changeRepository(
        directory = "pyfeelpp-tests/alg/test_alg_vector")
    wc = feelpp.Environment.worldCommPtr()
    u = feelpp.VectorPetscDouble(10, wc)
    u.setConstant(1)
    us = []
    for i in range(3):
        us.append( feelpp.VectorPetscDouble(10, wc) )
        us[i].setConstant(i+1)
    if feelpp.Environment.isMasterRank():
        print("Test mDot 3 vectors")
    r = u.mDot( us )
    assert( r.size() == 3 )
    for i in range(3):
        assert( r[i] == us[i].sum() )

    
   

def test_backend_vector():
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/alg/test_vector")
    wc = feelpp.Environment.worldCommPtr()
    if feelpp.Environment.isSequential():
        b = feelpp.backend(worldcomm=wc)
        dm = feelpp.DataMap(10,10,wc)
        v=b.newVector(dm)
        assert(v.size()==10)


def test_backend_matrix():
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


#def test_createFromPETSc(init_feelpp):
#    e = init_feelpp
#    feelpp.Environment.changeRepository(
#        directory="pyfeelpp-tests/alg/test_createFromPETSc")
#    from mpi4py import MPI
#    wc = MPI.COMM_WORLD
#
#    N = 150
#    v = PETSc.Vec().create(comm=wc)
#    v.setSizes(N)
#    v.setUp()
#    v.setFromOptions()
#    v.set(67)
#    v[45] = 450
#
#    v_petsc = feelpp.VectorPetscDouble(v)
#    v_petsc.vec().assemble()
#    v_petsc.vec().view()
#    assert( v_petsc.size()     == N   )
#    assert( v_petsc.vec()[45]  == 450 )
#    assert( v_petsc.vec()[148] == 67  )
#    
#     v_ublas = feelpp._alg.VectorUBlas.createFromPETSc(v_petsc)
#     assert( v_ublas.size() == N   )
#     assert( v_ublas[45]    == 450 )
#     assert( v_ublas[148]   == 67  )
