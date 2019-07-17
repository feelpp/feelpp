/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-18

  Copyright (C) 2005,2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_matrix.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-18
 */

#define BOOST_TEST_MODULE matrix testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{
namespace detail
{
template <typename T>
double myLocalProcessSum( Feel::Vector<T> const& vec )
{
    double res = 0;
    for ( size_type k=0;k<vec.map().nLocalDofWithGhost();++k )
        res += vec( k );
    return res;
}
}
}

FEELPP_ENVIRONMENT_NO_OPTIONS
using namespace Feel;

BOOST_AUTO_TEST_SUITE( matrix )

BOOST_AUTO_TEST_CASE( test_matrix_petsc_base )
{
    //#if defined(FEELPP_HAS_PETSC_H)

    int nprocs = Environment::worldComm().size();
    int m = 8;//10*Application::nProcess();
    int n = 8;//10*Application::nProcess();

    VectorPetsc<double> vec;
    MatrixPetsc<double> mat;
    std::cout << "is initialized ? " << mat.isInitialized() << "\n";
    mat.init( m*n, m*n, m*n/nprocs, m*n );
    std::cout << "is initialized ? " << mat.isInitialized() << "\n";

    /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
     */
    int Istart;
    int Iend;
    MatGetOwnershipRange( mat.mat(),&Istart,&Iend );
    VLOG(1) << "Istart = "<< Istart << "\n";
    VLOG(1) << "Iend = "<< Iend << "\n";

    /*
     Set matrix elements for the 2-D, five-point stencil in parallel.
     - Each processor needs to insert only elements that it owns
     locally (but any non-local elements will be sent to the
     appropriate processor during matrix assembly).
     - Always specify global rows and columns of matrix entries.

     Note: this uses the less common natural ordering that orders first
     all the unknowns for x = h then for x = 2h etc; Hence you see J = I +- n
     instead of J = I +- m as you might expect. The more standard ordering
     would first do all variables for y = h, then y = 2h etc.

     */

    for ( int I=Istart; I<Iend; I++ )
    {
        int J;
        double v = -1.0;
        int i = I/n;
        int j = I - i*n;
        VLOG(1) << "I= " << I << "\n";

        if ( i>0 )
        {
            J = I - n+1;
            VLOG(1) << "1 J= " << J << "\n";
            mat.set( I,J,v );
        }

        if ( i<m-1 )
        {
            J = I + n-1;
            VLOG(1) << "2 J= " << J << "\n";
            mat.set( I,J,v );
        }

        if ( j>0 )
        {
            J = I - 1;
            VLOG(1) << "3 J= " << J << "\n";
            mat.set( I,J,v );
        }

        if ( j<n-1 )
        {
            J = I + 1;
            VLOG(1) << "4 J= " << J << "\n";
            mat.set( I,J,v );
        }

        v = 4.0;
        mat.set( I,I,v );
    }

    VLOG(1) << "closing petsc matrix\n";
    mat.close();
    VLOG(1) << "closing petsc matrix done\n";

    VLOG(1) << "saving petsc matrix in matlab\n";
    //mat.printMatlab("m");
    //mat.printMatlab(std::string("/tmp/mat.m") );
    mat.printMatlab( "mat.m" );
    VLOG(1) << "saving petsc matrix in matlab done\n";
    //#endif
}
BOOST_AUTO_TEST_CASE( test_matrix_petsc_operations )
{
    double tolCheck = 1e-9;
    rank_type myrank = Environment::rank();
    rank_type worldsize = Environment::numberOfProcessors();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    double meshMeasure = mesh->measure();
    auto Vh1 = Pch<2>( mesh );
    size_type nDofVh1 = Vh1->nDof();
    size_type nLocalDofWithGhostVh1 = Vh1->nLocalDofWithGhost();
    auto u1 = Vh1->element();
    auto u1b = Vh1->element();
    auto backendPetsc = backend(_kind="petsc");
    auto vec1a = backendPetsc->newVector( Vh1 );
    auto vec1b = backendPetsc->newVector( Vh1 );
    auto mat1a = backendPetsc->newMatrix( _test=Vh1,_trial=Vh1 );
    auto mat1b = backendPetsc->newMatrix( _test=Vh1,_trial=Vh1 );
    auto mat1c = backendPetsc->newMatrix( _test=Vh1,_trial=Vh1 );
    size_type nDofActivesMat1Row = mat1a->mapRow().nLocalDofWithoutGhost();
    size_type nDofGhostsMat1aRow = mat1a->mapRow().nLocalGhosts();

    // l1Norm, linftyNorm
    for ( size_type k=0;k<nDofActivesMat1Row;++k )
        mat1a->set( k,k, 3.);
    mat1a->close();
    BOOST_CHECK_SMALL( mat1a->l1Norm() - 3., tolCheck );
    BOOST_CHECK_SMALL( mat1a->linftyNorm() - 3., tolCheck );
    // add scalar (on ghost only)
    for ( size_type k=0;k<nDofActivesMat1Row;++k )
        mat1a->add( k,k, 5.);
    mat1a->close();
    BOOST_CHECK_SMALL( mat1a->l1Norm() - 8., tolCheck );
    // set values
    for ( size_type k=0;k<nDofActivesMat1Row;++k )
        mat1a->set( k,k, myrank);
    mat1a->close();
    BOOST_CHECK_SMALL( mat1a->l1Norm() - (worldsize-1), tolCheck );
    // scale
    mat1a->scale(1./3.);
    BOOST_CHECK_SMALL( mat1a->l1Norm() - (worldsize-1)/3., tolCheck );
    mat1a->scale(3.);
    // add values (on ghost only)
    for ( size_type k=nDofActivesMat1Row;k<(nDofActivesMat1Row+nDofGhostsMat1aRow);++k )
        mat1a->add( k,k, 4.);
    mat1a->close();
    double valExact = worldsize-1;
    if ( Environment::worldComm().size() > 1 )
        valExact += 4.;
    BOOST_CHECK( mat1a->l1Norm() >= (valExact-tolCheck) );

    // energy
    vec1a->setConstant( 2. );
    vec1b->setConstant( 3. );
    form2(_test=Vh1,_trial=Vh1,_matrix=mat1a ) =
        integrate(_range=elements(mesh),_expr=idt(u1)*id(u1) );
    mat1a->close();
    BOOST_CHECK_CLOSE( mat1a->energy(vec1a,vec1b), 2*3*meshMeasure , tolCheck );
    // energy with ublas vector
    u1.setConstant( 3. );
    u1b.setConstant( 4. );
    BOOST_CHECK_CLOSE( mat1a->energy(u1,u1b), 3*4*meshMeasure , tolCheck );

    // add matrix
    form2(_test=Vh1,_trial=Vh1,_matrix=mat1b ) =
        integrate(_range=elements(mesh),_expr=idt(u1)*id(u1) );
    mat1b->close();
    mat1a->addMatrix(3,mat1b );
    BOOST_CHECK_SMALL( mat1a->energy(vec1a,vec1b) - (1+3)*2*3*meshMeasure , tolCheck );

    // diagonal
    mat1a->zero();
    for ( size_type k=0;k<nDofActivesMat1Row;++k )
        mat1a->set( k,k, 3.);
    mat1a->close();
    mat1a->diagonal(vec1a);
    BOOST_CHECK_SMALL( vec1a->sum() - 3*nDofVh1 , tolCheck );
    size_type diagsum=0;
    for ( size_type k=0;k<nLocalDofWithGhostVh1;++k)
        diagsum += vec1a->operator()(k);
    BOOST_CHECK( diagsum == 3*nLocalDofWithGhostVh1 );
    for ( size_type k=0;k<nLocalDofWithGhostVh1;++k)
        vec1a->set( k,4. );
    vec1a->close();
    BOOST_CHECK_SMALL( vec1a->sum() - 4*nDofVh1 , tolCheck );
    // diagonal as new vector
    auto vec1ad = mat1a->diagonal();
    BOOST_CHECK_SMALL( vec1ad->sum() - 3*nDofVh1 , tolCheck );
    diagsum=0;
    for ( size_type k=0;k<nLocalDofWithGhostVh1;++k)
        diagsum += vec1ad->operator()(k);
    BOOST_CHECK( diagsum == 3*nLocalDofWithGhostVh1 );

    // transpose
    mat1a->zero();
    for ( size_type k=0;k<nDofActivesMat1Row;++k )
        mat1a->set( k,k, 7.);
    mat1a->close();
    mat1b->zero();
    mat1a->transpose( mat1b );
    BOOST_CHECK_SMALL( mat1b->l1Norm() - 7., tolCheck );
    // transpose with matrix clear
    mat1b->clear();
    mat1a->transpose( mat1b );
    BOOST_CHECK_SMALL( mat1b->l1Norm() - 7., tolCheck );
    // transpose with new matrix
    auto mat1at = mat1a->transpose();
    BOOST_CHECK_SMALL( mat1at->l1Norm() - 7., tolCheck );

    // symmetricPart
    form2(_test=Vh1,_trial=Vh1,_matrix=mat1a ) =
        integrate(_range=elements(mesh),_expr=gradt(u1)*vec(cst(3.),cst(1.))*id(u1) );
    mat1a->close();
    mat1a->symmetricPart( mat1b );
    mat1a->transpose( mat1c );
    mat1c->addMatrix(1.,mat1a);
    mat1c->scale( 0.5 );
    BOOST_CHECK_SMALL( mat1b->l1Norm() - mat1c->l1Norm(), tolCheck );

    // matMatMult
    mat1a->zero();
    mat1b->zero();
    for ( size_type k=0;k<nDofActivesMat1Row;++k )
    {
        mat1a->set( k,k, (7.+k));
        mat1b->set( k,k, 1./(7.+k));
    }
    mat1a->close();
    mat1b->close();
    mat1c->clear(); // stencil will change
    mat1a->matMatMult( *mat1b, *mat1c );
    BOOST_CHECK_SMALL( mat1c->l1Norm() - 1., tolCheck );

    // re-check local assembly
    mat1a->zero();
    mat1b->zero();
    mat1c->zero();
    for ( size_type k=0;k<nDofActivesMat1Row;++k )
    {
        mat1a->set( k,k, 3. );
        mat1b->set( k,k, 1./(7.+k));
        mat1c->set( k,k, 1./(7.+k));
    }
    mat1a->close();
    mat1b->close();
    mat1c->close();

    // multVector
    vec1a->setConstant( 2. );
    mat1a->multVector( vec1a, vec1b );
    BOOST_CHECK_CLOSE( vec1b->sum(), 2*3*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vec1b), 2*3*nLocalDofWithGhostVh1, tolCheck );
    // multVector with vector clear
    vec1b->clear();
    mat1a->multVector( vec1a, vec1b );
    BOOST_CHECK_CLOSE( vec1b->sum(), 2*3*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vec1b), 2*3*nLocalDofWithGhostVh1, tolCheck );
    // multVector with ublas vector
    u1.setConstant( 5. );
    mat1a->multVector( u1, *vec1b );
    BOOST_CHECK_CLOSE( vec1b->sum(), 5*3*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vec1b), 5*3*nLocalDofWithGhostVh1, tolCheck );
    vec1a->setConstant( 7. );
    mat1a->multVector( *vec1a, u1 );
    BOOST_CHECK_CLOSE( u1.sum(), 7*3*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(u1), 7*3*nLocalDofWithGhostVh1, tolCheck );
    u1b.setConstant( 7. );
    mat1a->multVector( u1b, u1 );
    BOOST_CHECK_CLOSE( u1.sum(), 7*3*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(u1), 7*3*nLocalDofWithGhostVh1, tolCheck );

    // setDiagonal
    mat1a->zero();
    vec1a->setConstant( 8. );
    mat1a->setDiagonal( vec1a );
    BOOST_CHECK_CLOSE( mat1a->l1Norm(), 8., tolCheck );
    // setDiagonal with ublas vector
    u1.setConstant( 9. );
    mat1a->setDiagonal( u1 );
    BOOST_CHECK_CLOSE( mat1a->l1Norm(), 9., tolCheck );
    // addDiagonal
    vec1a->setConstant( 8. );
    mat1a->addDiagonal( vec1a );
    BOOST_CHECK_CLOSE( mat1a->l1Norm(), 9+8, tolCheck );
    // addDiagonal with ublas vector
    u1.setConstant( 7. );
    mat1a->addDiagonal( u1 );
    BOOST_CHECK_CLOSE( mat1a->l1Norm(), 9+8+7, tolCheck );

    // PtAP
    vec1b->setConstant( 3. );
    mat1b->setDiagonal( vec1b );
    mat1c->clear(); // stencil will change
    mat1a->PtAP( *mat1b, *mat1c );
    BOOST_CHECK_CLOSE( mat1c->l1Norm(), 3*(9+8+7)*3, tolCheck );
    mat1a->PtAP( *mat1b, *mat1c );// reuse stencil
    BOOST_CHECK_CLOSE( mat1c->l1Norm(), 3*(9+8+7)*3, tolCheck );
    // PAPt
    vec1b->setConstant( 5. );
    mat1b->setDiagonal( vec1b );
    mat1c->clear(); // stencil will change
    mat1a->PAPt( *mat1b, *mat1c );
    BOOST_CHECK_CLOSE( mat1c->l1Norm(), 5*(9+8+7)*5, tolCheck );
    mat1a->PAPt( *mat1b, *mat1c ); // reuse stencil
    BOOST_CHECK_CLOSE( mat1c->l1Norm(), 5*(9+8+7)*5, tolCheck );

}
BOOST_AUTO_TEST_SUITE_END()
