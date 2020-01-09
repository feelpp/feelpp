/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2015-04-18

  Copyright (C) 2015 UJF

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
   \file test_matrix_block.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2015-04-18
 */
#if 1
#define BOOST_TEST_MODULE test_matrix_block
#include <feel/feelcore/testsuite.hpp>
#else
#define USE_BOOST_TEST 0
#endif

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>


namespace test_matrix_block
{
using namespace Feel;

void run()
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>> );
    auto Vhu = Pchv<2>( mesh );
    auto Vhp = Pch<1>( mesh );
    auto Vhl = Pch<0>( mesh );
    auto Vht = Pch<1>( mesh );
    auto u = Vhu->element();
    auto p = Vhp->element();
    auto l = Vhl->element();
    auto t = Vht->element();

    // first version for build the block matrix
    auto a_uu = form2( _test=Vhu, _trial=Vhu);
    auto a_up = form2( _test=Vhu, _trial=Vhp);
    auto a_pu = form2( _test=Vhp, _trial=Vhu);
    auto a_pl = form2( _test=Vhp, _trial=Vhl);
    auto a_lp = form2( _test=Vhl, _trial=Vhp);
    auto a_tt = form2( _test=Vht, _trial=Vht);
    a_uu = integrate(_range=elements(mesh),_expr=inner( gradt(u), grad(u) ) );
    a_up = integrate(_range=elements(mesh),_expr=-div(u)*idt(p) );
    a_pu = integrate(_range=elements(mesh),_expr=divt(u)*id(p) );
    a_pl = integrate(_range=elements(mesh),_expr=id(p)*idt(l) );
    a_lp = integrate(_range=elements(mesh),_expr=idt(p)*id(l) );
    a_tt = integrate(_range=elements(mesh),_expr=inner( gradt(t), grad(t) ) );
    BlocksBaseSparseMatrix<double> myblockMat(4,4);
    myblockMat(0,0) = a_uu.matrixPtr();
    myblockMat(0,1) = a_up.matrixPtr();
    myblockMat(1,0) = a_pu.matrixPtr();
    myblockMat(1,2) = a_pl.matrixPtr();
    myblockMat(2,1) = a_lp.matrixPtr();
    myblockMat(3,3) = a_tt.matrixPtr();
    auto A = backend()->newBlockMatrix(_block=myblockMat, _copy_values=true);

    // second version for build the block matrix : more efficient and less memory
    BlocksBaseGraphCSR myblockGraph(4,4);
    myblockGraph(0,0) = stencil(_test=Vhu,_trial=Vhu, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(0,1) = stencil(_test=Vhu,_trial=Vhp, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,0) = stencil(_test=Vhp,_trial=Vhu, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,2) = stencil(_test=Vhp,_trial=Vhl, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(2,1) = stencil(_test=Vhl,_trial=Vhp, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(3,3) = stencil(_test=Vht,_trial=Vht, _diag_is_nonzero=false, _close=false)->graph();
    auto Abis = backend()->newBlockMatrix(_block=myblockGraph);
    size_type indexStart_u = 0;
    size_type indexStart_p = 1;
    size_type indexStart_l = 2;
    size_type indexStart_t = 3;
    form2( _test=Vhu, _trial=Vhu, _matrix=Abis,_rowstart=indexStart_u,_colstart=indexStart_u )
        += integrate(_range=elements(mesh),_expr=inner( gradt(u), grad(u) ) );
    form2( _test=Vhu, _trial=Vhp, _matrix=Abis,_rowstart=indexStart_u,_colstart=indexStart_p )
        += integrate(_range=elements(mesh),_expr=-div(u)*idt(p) );
    form2( _test=Vhp, _trial=Vhu, _matrix=Abis,_rowstart=indexStart_p,_colstart=indexStart_u )
        += integrate(_range=elements(mesh),_expr=divt(u)*id(p) );
    form2( _test=Vhp, _trial=Vhl, _matrix=Abis,_rowstart=indexStart_p,_colstart=indexStart_l )
        += integrate(_range=elements(mesh),_expr=id(p)*idt(l) );
    form2( _test=Vhl, _trial=Vhp, _matrix=Abis,_rowstart=indexStart_l,_colstart=indexStart_p )
        += integrate(_range=elements(mesh),_expr=idt(p)*id(l) );
    form2( _test=Vht, _trial=Vht, _matrix=Abis,_rowstart=indexStart_t,_colstart=indexStart_t )
        += integrate(_range=elements(mesh),_expr=inner( gradt(t), grad(t) ) );

    // check that matrix are identicaly
    Abis->addMatrix(-1.0, A );
    double evall1Norm = Abis->l1Norm();
    double evallinftyNorm = Abis->linftyNorm();
#if USE_BOOST_TEST
    BOOST_CHECK_SMALL( evall1Norm, 1e-9 );
    BOOST_CHECK_SMALL( evallinftyNorm, 1e-9 );
#else
    CHECK( std::abs( evall1Norm ) < 1e-9 ) << "must be zero : " << evall1Norm;
    CHECK( std::abs( evallinftyNorm ) < 1e-9 ) << "must be zero : " << evallinftyNorm;
#endif

    // extract submatrix : each block
    auto const& indiceExtract_u = A->mapRow().dofIdToContainerId( indexStart_u );
    auto const& indiceExtract_p = A->mapRow().dofIdToContainerId( indexStart_p );
    auto const& indiceExtract_l = A->mapRow().dofIdToContainerId( indexStart_l );
    auto const& indiceExtract_t = A->mapRow().dofIdToContainerId( indexStart_t );
    auto A_uu = A->createSubMatrix(indiceExtract_u,indiceExtract_u);
    auto A_up = A->createSubMatrix(indiceExtract_u,indiceExtract_p);
    auto A_pu = A->createSubMatrix(indiceExtract_p,indiceExtract_u);
    auto A_pl = A->createSubMatrix(indiceExtract_p,indiceExtract_l);
    auto A_lp = A->createSubMatrix(indiceExtract_l,indiceExtract_p);
    auto A_tt = A->createSubMatrix(indiceExtract_t,indiceExtract_t);
    // check
#if USE_BOOST_TEST
    BOOST_CHECK_CLOSE( A_uu->l1Norm(), a_uu.matrixPtr()->l1Norm(), 1e-9 );
    BOOST_CHECK_CLOSE( A_up->l1Norm(), a_up.matrixPtr()->l1Norm(), 1e-9 );
    BOOST_CHECK_CLOSE( A_pu->l1Norm(), a_pu.matrixPtr()->l1Norm(), 1e-9 );
    BOOST_CHECK_CLOSE( A_pl->l1Norm(), a_pl.matrixPtr()->l1Norm(), 1e-9 );
    BOOST_CHECK_CLOSE( A_lp->l1Norm(), a_lp.matrixPtr()->l1Norm(), 1e-9 );
    BOOST_CHECK_CLOSE( A_tt->l1Norm(), a_tt.matrixPtr()->l1Norm(), 1e-9 );
#else
    CHECK( std::abs( A_uu->l1Norm() - a_uu.matrixPtr()->l1Norm() ) < 1e-9 ) << "must be identicaly";
    CHECK( std::abs( A_up->l1Norm() - a_up.matrixPtr()->l1Norm() ) < 1e-9 ) << "must be identicaly";
    CHECK( std::abs( A_pu->l1Norm() - a_pu.matrixPtr()->l1Norm() ) < 1e-9 ) << "must be identicaly";
    CHECK( std::abs( A_pl->l1Norm() - a_pl.matrixPtr()->l1Norm() ) < 1e-9 ) << "must be identicaly";
    CHECK( std::abs( A_lp->l1Norm() - a_lp.matrixPtr()->l1Norm() ) < 1e-9 ) << "must be identicaly";
    CHECK( std::abs( A_tt->l1Norm() - a_tt.matrixPtr()->l1Norm() ) < 1e-9 ) << "must be identicaly";
#endif


    BlocksBaseSparseMatrix<double> myblockMat1(4,4);
    myblockMat1(0,0) = A_uu;
    myblockMat1(0,1) = A_up;
    myblockMat1(1,0) = A_pu;
    myblockMat1(1,2) = A_pl;
    myblockMat1(2,1) = A_lp;
    myblockMat1(3,3) = A_tt;
    auto Afromsubmat = backend()->newBlockMatrix(_block=myblockMat1, _copy_values=true);
    // check
#if USE_BOOST_TEST
    BOOST_CHECK_CLOSE( Afromsubmat->l1Norm(), A->l1Norm(), 1e-9 );
#else
    CHECK( std::abs( Afromsubmat->l1Norm() - A->l1Norm() ) < 1e-9 ) << "must be identicaly";
#endif

    A_uu->addMatrix(-1.0, a_uu.matrixPtr() );
    A_up->addMatrix(-1.0, a_up.matrixPtr() );
    A_pu->addMatrix(-1.0, a_pu.matrixPtr() );
    A_pl->addMatrix(-1.0, a_pl.matrixPtr() );
    A_lp->addMatrix(-1.0, a_lp.matrixPtr() );
    A_tt->addMatrix(-1.0, a_tt.matrixPtr() );
#if USE_BOOST_TEST
    BOOST_CHECK_SMALL( A_uu->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( A_up->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( A_pu->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( A_pl->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( A_lp->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( A_tt->l1Norm(), 1e-9 );
#else
    CHECK( std::abs( A_uu->l1Norm() ) < 1e-9 ) << "must be zero";
    CHECK( std::abs( A_up->l1Norm() ) < 1e-9 ) << "must be zero";
    CHECK( std::abs( A_pu->l1Norm() ) < 1e-9 ) << "must be zero";
    CHECK( std::abs( A_pl->l1Norm() ) < 1e-9 ) << "must be zero";
    CHECK( std::abs( A_lp->l1Norm() ) < 1e-9 ) << "must be zero";
    CHECK( std::abs( A_tt->l1Norm() ) < 1e-9 ) << "must be zero";
#endif

    // extract another submatrix : a set of block
    std::vector<uint32_type> indiceExtract_upl;
    for ( size_type k=0;k<Vhu->nLocalDofWithGhost();++k)
        indiceExtract_upl.push_back( A->mapRow().dofIdToContainerId(indexStart_u,k) );
    for ( size_type k=0;k<Vhp->nLocalDofWithGhost();++k)
        indiceExtract_upl.push_back( A->mapRow().dofIdToContainerId(indexStart_p,k) );
    for ( size_type k=0;k<Vhl->nLocalDofWithGhost();++k)
        indiceExtract_upl.push_back( A->mapRow().dofIdToContainerId(indexStart_l,k) );
    auto A_uplupl = A->createSubMatrix(indiceExtract_upl,indiceExtract_upl);
    // build another block matrix for compare with extracted matrix
    BlocksBaseSparseMatrix<double> myblockMat2(3,3);
    myblockMat2(0,0) = a_uu.matrixPtr();
    myblockMat2(0,1) = a_up.matrixPtr();
    myblockMat2(1,0) = a_pu.matrixPtr();
    myblockMat2(1,2) = a_pl.matrixPtr();
    myblockMat2(2,1) = a_lp.matrixPtr();
    auto A_upluplbis = backend()->newBlockMatrix(_block=myblockMat2, _copy_values=true);
    // check
#if USE_BOOST_TEST
    BOOST_CHECK_CLOSE( A_uplupl->l1Norm(), A_upluplbis->l1Norm(), 1e-9 );
#else
    CHECK( std::abs( A_uplupl->l1Norm() - A_upluplbis->l1Norm() ) < 1e-9 ) << "must be identicaly";
#endif
    A_upluplbis->addMatrix(-1.0, A_uplupl );
#if USE_BOOST_TEST
    BOOST_CHECK_SMALL( A_upluplbis->l1Norm(), 1e-9 );
#else
    CHECK( std::abs( A_upluplbis->l1Norm() ) < 1e-9 ) << "must be zero";
#endif

}

} // namespace test_matrix_block

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( matrix_block )

BOOST_AUTO_TEST_CASE( matrix_block )
{
    test_matrix_block::run();
}

BOOST_AUTO_TEST_SUITE_END()

#else

int
main( int argc, char* argv[] )
{
    Feel::Environment env( argc,argv );
    test_matrix_block::run();
}

#endif
