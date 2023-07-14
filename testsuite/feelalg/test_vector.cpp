/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-08-07

  Copyright (C) 2010 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_vector.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2010-08-07
 */
//#include <cmath>

#define BOOST_TEST_MODULE vector testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/traits.hpp>
#include <feel/feelcore/enumerate.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/thch.hpp>
#include <feel/feeldiscr/stencil.hpp>

#include <range/v3/view/take.hpp>

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

BOOST_AUTO_TEST_SUITE( vector )

BOOST_AUTO_TEST_CASE( test_vector_ublas_base )
{
    BOOST_TEST_MESSAGE( "test_vector_ublas_base" );
    using namespace Feel;
    VectorUblas<double> v1( 100 ), v2( 100 ), v3( 100 );
    v1.setConstant( 1 );
    v2.setConstant( 2 );
    BOOST_CHECK_CLOSE( v1.sqrt().sum(), v1.size(), 1e-10 );
    BOOST_CHECK_CLOSE( v2.sqrt().sum(), sqrt( 2 )*v1.size(), 1e-10 );
    BOOST_CHECK_CLOSE( v2.pow( 2 ).sqrt().sum(), 2*v1.size(), 1e-10 );
    v3.setZero();
    v3 = element_product( v1, v2 );
    BOOST_CHECK_CLOSE( v3.sqrt().sum(), sqrt( 2*1 )*v1.size(), 1e-10 );
    BOOST_CHECK_CLOSE( v3.sum(), 2*1*v1.size(), 1e-10 );
}
BOOST_AUTO_TEST_CASE( test_vector_ublas_operations )
{
    BOOST_TEST_MESSAGE( "test_vector_ublas_operations" );
    using namespace Feel;
    double tolCheck = 1e-9;
    auto mesh = unitSquare();
    auto Vh1 = Pchv<2>( mesh );
    auto v1 = Vh1->element();
    v1.setConstant( 2 );
    auto v1bis = Vh1->element();
    v1bis.setConstant( 3 );

    auto Vh2 = THch<2>( mesh );
    auto V2 = Vh2->element();
    auto v2a = V2.element<0>();
    auto v2b = V2.element<1>();
    v2a.setConstant( 5 );
    v2b.setConstant( 7 );

    auto Vh3 = Vh2->functionSpace<0>();
    auto v3 = Vh3->element();
    v3.setConstant( 9 );

    auto backendPetsc = backend(_kind="petsc");
    auto v_petsc1 = backendPetsc->newVector( Vh1 );
    auto v_petsc2 = backendPetsc->newVector( Vh2 );
    auto v_petsc3 = backendPetsc->newVector( Vh3 );


    size_type nDofVh1 = Vh1->nDof();
    size_type nLocalDofWithGhostVh1 = Vh1->nLocalDofWithGhost();
    size_type nDofVh2 = Vh2->nDof();
    size_type nLocalDofWithGhostVh2 = Vh2->nLocalDofWithGhost();
    size_type nDofVh2a = v2a.functionSpace()->nDof();
    size_type nLocalDofWithGhostVh2a = v2a.functionSpace()->nLocalDofWithGhost();
    size_type nDofVh2b = v2b.functionSpace()->nDof();
    size_type nLocalDofWithGhostVh2b = v2b.functionSpace()->nLocalDofWithGhost();
    size_type nDofVh3 = Vh3->nDof();
    BOOST_CHECK( nDofVh2a == nDofVh3 && v2a.nLocalDof() == v3.nLocalDof() );

    // sum
    BOOST_CHECK( v1.sum() == 2*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v1), 2*nLocalDofWithGhostVh1, tolCheck );
    BOOST_CHECK( v2a.sum() == 5*nDofVh2a );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v2a), 5*nLocalDofWithGhostVh2a, tolCheck );
    BOOST_CHECK( v2b.sum() == 7*nDofVh2b );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v2b), 7*nLocalDofWithGhostVh2b, tolCheck );
    // l1Norm
    BOOST_CHECK( v1.l1Norm() == 2*nDofVh1 );
    BOOST_CHECK( v2a.l1Norm() == 5*nDofVh2a );
    BOOST_CHECK( v2b.l1Norm() == 7*nDofVh2b );
    //l2Norm
    BOOST_CHECK_SMALL( v1.l2Norm()-std::sqrt(4*nDofVh1), 1e-9 );
    BOOST_CHECK_SMALL( v2a.l2Norm()-std::sqrt(25*nDofVh2a), 1e-9 );
    BOOST_CHECK_SMALL( v2b.l2Norm()-std::sqrt(49*nDofVh2b), 1e-9 );
    // element_product
    auto prod1 = element_product(v1,v1);
    BOOST_CHECK( prod1.sum() == 4*nDofVh1 );
    auto prod2a = element_product(v2a,v2a);
    BOOST_CHECK( prod2a.sum() == 25*nDofVh2a );
    auto prod2b = element_product(v2b,v2b);
    BOOST_CHECK( prod2b.sum() == 49*nDofVh2b );
    auto prod2a3 = element_product(v2a,v3);
    BOOST_CHECK( prod2a3.sum() == 45*nDofVh2a );
    // add scalar
    v1.add( 3 );
    BOOST_CHECK( v1.sum() == 5*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v1), 5*nLocalDofWithGhostVh1, tolCheck );
    v2a.add( 3 );
    BOOST_CHECK( v2a.sum() == 8*nDofVh2a );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v2a), 8*nLocalDofWithGhostVh2a, tolCheck );
    v2b.add( 3 );
    BOOST_CHECK( v2b.sum() == 10*nDofVh2b );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v2b), 10*nLocalDofWithGhostVh2b, tolCheck );
    // operator=
    v1 = v1bis;
    BOOST_CHECK( v1.sum() == 3*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v1), 3*nLocalDofWithGhostVh1, tolCheck );
    v3 = v2a;
    BOOST_CHECK( v3.sum() == 8*nDofVh2a );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v3), 8*nLocalDofWithGhostVh2a, tolCheck );
    v3.setConstant( 9 );
    v2a = v3;
    BOOST_CHECK( v2a.sum() == 9*nDofVh2a );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v2a), 9*nLocalDofWithGhostVh2a, tolCheck );
    // operator= (with petsc vector)
    v_petsc1->setConstant( 2 );
    v_petsc3->setConstant( 4 );
    v1 = *v_petsc1;
    BOOST_CHECK( v1.sum() == 2*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v1), 2*nLocalDofWithGhostVh1, tolCheck );
    v3 = *v_petsc3;
    BOOST_CHECK( v3.sum() == 4*nDofVh2a );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v3), 4*nLocalDofWithGhostVh2a, tolCheck );
    v2a = *v_petsc3;
    BOOST_CHECK( v2a.sum() == 4*nDofVh2a );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v2a), 4*nLocalDofWithGhostVh2a, tolCheck );
    // operator= (from petsc vector to element type not init)
    decltype(Vh1)::element_type::element_type u1NotInit;
    u1NotInit = *v_petsc1;
    BOOST_CHECK( u1NotInit.sum() == 2*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v1), 2*nLocalDofWithGhostVh1, tolCheck );
    decltype(Vh2)::element_type::element_type u2NotInit;
    v_petsc2->setConstant( 4 );
    u2NotInit = *v_petsc2;
    BOOST_CHECK( u2NotInit.sum() == 4*nDofVh2 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(u2NotInit), 4*nLocalDofWithGhostVh2, tolCheck );
    decltype(Vh3)::element_type::element_type u3NotInit;
    v2a.setConstant( 5 );
    u3NotInit = v2a;
    BOOST_CHECK( u3NotInit.sum() == 5*nDofVh3 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(u3NotInit), 5*nLocalDofWithGhostVh2a, tolCheck );

    // add vector
    v1.setConstant( 2 );
    v2a.setConstant( 5 );
    v2b.setConstant( 7 );
    v3.setConstant( 8 );
    v1.add(-1.,v1 );
    BOOST_CHECK_CLOSE( v1.sum(), 0., tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v1), 0., tolCheck );
    v2a.add( 2., v3 );
    BOOST_CHECK_CLOSE( v2a.sum(), (5+2*8)*nDofVh2a, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v2a), (5+2*8)*nLocalDofWithGhostVh2a, tolCheck );
    v3.add( -3., v2a );
    BOOST_CHECK_CLOSE( v3.sum(), (8.-3.*(5+2*8))*nDofVh3, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v3), (8.-3.*(5+2*8))*nLocalDofWithGhostVh2a, tolCheck );
    // add vector (with petsc vector)
    v1.setConstant( 2 );
    v2a.setConstant( 5 );
    v3.setConstant( 7 );
    v_petsc1->setConstant( 4 );
    v_petsc3->setConstant( 6 );
    v1.add(3.,*v_petsc1 );
    BOOST_CHECK( v1.sum() == (2+3*4)*nDofVh1 );
    v2a.add(3.,*v_petsc3 );
    BOOST_CHECK( v2a.sum() == (5+3*6)*nDofVh2a );
    v3.add(3.,*v_petsc3 );
    BOOST_CHECK( v3.sum() == (7+3*6)*nDofVh2a );

    // linftyNorm
    rank_type myrank = Environment::rank();
    rank_type worldsize = Environment::numberOfProcessors();
    v1.setConstant( 2 + myrank );
    sync(v1);
    BOOST_CHECK( v1.linftyNorm() == (2+worldsize-1) );
    v2a.setConstant( 5 + myrank );
    sync(v2a);
    BOOST_CHECK( v2a.linftyNorm() == (5+worldsize-1) );
    v2b.setConstant( 7 + myrank );
    sync(v2b);
    BOOST_CHECK( v2b.linftyNorm() == (7+worldsize-1) );
    // max
    BOOST_CHECK( v1.max() == (2+worldsize-1) );
    BOOST_CHECK( v2a.max() == (5+worldsize-1) );
    BOOST_CHECK( v2b.max() == (7+worldsize-1) );
    // min
    BOOST_CHECK( v1.min() == 2 );
    BOOST_CHECK( v2a.min() == 5 );
    BOOST_CHECK( v2b.min() == 7 );
    // dot
    v1.setConstant( 2 );
    v1bis.setConstant( 3 );
    v2a.setConstant( 5 );
    v3.setConstant( 7 );
    BOOST_CHECK( v1.dot(v1bis) == (3*2)*nDofVh1 );
    BOOST_CHECK( v2a.dot(v3) == (5*7)*nDofVh2a );
    BOOST_CHECK( v3.dot(v2a) == (5*7)*nDofVh2a );
    // dot (with petsc vector)
    v_petsc1->setConstant( 4 );
    v_petsc3->setConstant( 6 );
    BOOST_CHECK( v1.dot( *v_petsc1 ) == (2*4)*nDofVh1 );
    BOOST_CHECK( v2a.dot( *v_petsc3 ) == (5*6)*nDofVh2a );
    BOOST_CHECK( v3.dot( *v_petsc3 ) == (7*6)*nDofVh2a );
    //inner_product
    BOOST_CHECK( inner_product( v1, v1bis ) == (3*2)*nDofVh1 );
    BOOST_CHECK( inner_product( v1,*v_petsc1 ) == (2*4)*nDofVh1 );
}
BOOST_AUTO_TEST_CASE( test_vector_ublas_extarray )
{
    BOOST_TEST_MESSAGE( "test_vector_ublas_extarray" );
    using namespace Feel;
    double tolCheck = 1e-9;
    auto mesh = unitSquare();
    auto Vh = Pchv<2>( mesh );
    auto Wh = THch<2>( mesh );
    auto Xh = Wh->functionSpace<1>();
    size_type nDofVh = Vh->nDof(), nDofWh = Wh->nDof(), nDofXh = Xh->nDof();
    size_type nLocalDofWithGhostVh = Vh->nLocalDofWithGhost();
    size_type nLocalDofWithGhostWh = Wh->nLocalDofWithGhost();
    size_type nLocalDofWithGhostXh = Xh->nLocalDofWithGhost();

    auto backendPetsc = backend(_kind="petsc");
    auto v_petsc1 = backendPetsc->newVector( Vh );
    auto v_petsc2 = backendPetsc->newVector( Vh );
    auto w_petsc1 = backendPetsc->newVector( Wh );
    auto x_petsc1 = backendPetsc->newVector( Xh );
    auto x_petsc2 = backendPetsc->newVector( Xh );

    auto v_ublas1 = Vh->element( v_petsc1 );
    auto v_ublas2 = Vh->element( v_petsc2 );
    auto w_ublas1 = Wh->element( w_petsc1 );
    auto w1_ublas1 = w_ublas1.element<1>();
    auto x_ublas1 = Xh->element( x_petsc1 );

    v_petsc1->setConstant( 2 );
    w_petsc1->setConstant( 3 );
    x_petsc1->setConstant( 4 );
    // sum
    BOOST_CHECK( v_ublas1.sum() == 2*nDofVh );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(v_ublas1), 2*nLocalDofWithGhostVh, tolCheck );
    BOOST_CHECK( w_ublas1.sum() == 3*nDofWh );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(w_ublas1), 3*nLocalDofWithGhostWh, tolCheck );
    BOOST_CHECK( x_ublas1.sum() == 4*nDofXh );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(x_ublas1), 4*nLocalDofWithGhostXh, tolCheck );
    // add scalar
    v_ublas1.add( 3. );
    w_ublas1.add( 4. );
    x_ublas1.add( 5. );

    BOOST_CHECK( v_ublas1.sum() == (2+3)*nDofVh );
    BOOST_CHECK( v_petsc1->sum() == (2+3)*nDofVh );
    BOOST_CHECK( w_ublas1.sum() == (3+4)*nDofWh );
    BOOST_CHECK( w_petsc1->sum() == (3+4)*nDofWh );
    BOOST_CHECK( x_ublas1.sum() == (4+5)*nDofXh );
    BOOST_CHECK( x_petsc1->sum() == (4+5)*nDofXh );
    // add vector
    v_ublas2.setConstant( 2 );
    v_ublas1.add( 3., v_ublas2 );
    BOOST_CHECK( v_ublas1.sum() == (2+3+3*2)*nDofVh );
    w1_ublas1.setConstant( 6 );
    x_ublas1.setConstant( 4 );
    w1_ublas1.add( 2, x_ublas1 );
    BOOST_CHECK( w1_ublas1.sum() == (6+2*4)*nDofXh );
    x_ublas1.add( 5, w1_ublas1 );
    BOOST_CHECK( x_ublas1.sum() == (4+5*(6+2*4))*nDofXh );

    // operator=
    v_ublas2.setConstant( 2 );
    v_ublas1 = v_ublas2;
    BOOST_CHECK( v_ublas1.sum() == 2*nDofVh );
    x_ublas1.setConstant( 4 );
    w1_ublas1 = x_ublas1;
    BOOST_CHECK( w1_ublas1.sum() == 4*nDofXh );
    w1_ublas1.setConstant( 6 );
    x_ublas1 = w1_ublas1;
    BOOST_CHECK( x_ublas1.sum() == 6*nDofXh );

    // dot
    v_ublas1.setConstant( 3 );
    v_ublas2.setConstant( 5 );
    BOOST_CHECK( v_ublas1.dot(v_ublas2) == (3*5)*nDofVh );
    w1_ublas1.setConstant( 6 );
    x_ublas1.setConstant( 4 );
    BOOST_CHECK( w1_ublas1.dot(x_ublas1) == (6*4)*nDofXh );
    BOOST_CHECK( x_ublas1.dot(w1_ublas1) == (6*4)*nDofXh );

}
BOOST_AUTO_TEST_CASE( test_vector_petsc )
{
    BOOST_TEST_MESSAGE( "test_vector_petsc" );
    using namespace Feel;
    double tolCheck = 1e-9;
    auto mesh = unitSquare();
    auto Vh1 = Pchv<2>( mesh );
    auto v_ublas1 = Vh1->element();
    size_type nDofVh1 = Vh1->nDof();
    size_type nLocalDofWithGhostVh1 = Vh1->nLocalDofWithGhost();
    auto backendPetsc = backend(_kind="petsc");
    auto v_petsc1 = backendPetsc->newVector( Vh1 );
    auto v_petsc2 = backendPetsc->newVector( Vh1 );

    auto VhB = THch<2>( mesh );
    auto VhB1 = VhB->functionSpace<1>();
    size_type nDofVhB = VhB->nDof();
    size_type nDofVhB1 = VhB1->nDof();
    size_type nLocalDofWithGhostVhB1 = VhB1->nLocalDofWithGhost();
    auto vB_ublas1 = VhB->element();
    auto vB1_ublas1 = vB_ublas1.element<1>();
    auto vB_ublas2 = VhB->element();
    auto vB1_ublas2 = vB_ublas2.element<1>();
    auto vB1_ublas3 = VhB1->element();
    auto vB1_petsc1 = backendPetsc->newVector( VhB1 );
    auto vB1_petsc2 = backendPetsc->newVector( VhB1 );

    v_petsc1->setConstant( 2 );
    v_petsc2->setConstant( 5 );
    // sum
    BOOST_CHECK( v_petsc1->sum() == 2*nDofVh1 );
    BOOST_CHECK( v_petsc2->sum() == 5*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), 2*nLocalDofWithGhostVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc2), 5*nLocalDofWithGhostVh1, tolCheck );
    // add scalar
    v_petsc1->add( 1. );
    BOOST_CHECK( v_petsc1->sum() == (2+1)*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), (2+1)*nLocalDofWithGhostVh1, tolCheck );
    // add vector
    v_petsc1->add( 3., v_petsc2 );
    BOOST_CHECK( v_petsc1->sum() == (2+1+3*5)*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), (2+1+3*5)*nLocalDofWithGhostVh1, tolCheck );
    // add vector (with ublas vector)
    v_petsc1->setConstant( 5 );
    v_ublas1.setConstant( 4 );
    v_petsc1->add( 3., v_ublas1 );
    BOOST_CHECK( v_petsc1->sum() == (5+3*4)*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), (5+3*4)*nLocalDofWithGhostVh1, tolCheck );
    vB1_petsc1->setConstant( 3 );
    vB1_ublas1.setConstant( 9 );
    vB1_petsc1->add( 4., vB1_ublas1 );
    BOOST_CHECK( vB1_petsc1->sum() == (3+4*9)*nDofVhB1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (3+4*9)*nLocalDofWithGhostVhB1, tolCheck );

    // dot
    BOOST_CHECK( v_petsc1->dot(v_petsc2) == (5+3*4)*5*nDofVh1 );
    // dot (with ublas vector)
    v_petsc1->setConstant( 8 );
    v_ublas1.setConstant( 3 );
    BOOST_CHECK( v_petsc1->dot(v_ublas1) == (8*3)*nDofVh1 );
    vB1_petsc1->setConstant( 8 );
    vB1_ublas1.setConstant( 3 );
    BOOST_CHECK( vB1_petsc1->dot(vB1_ublas1) == (8*3)*nDofVhB1 );
    //inner_product
    BOOST_CHECK( inner_product( v_petsc1, v_petsc2 ) == (8*5)*nDofVh1 );
    BOOST_CHECK( inner_product( *v_petsc1, v_ublas1 ) == (8*3)*nDofVh1 );
    // reciprocal
    v_petsc1->reciprocal();
    BOOST_CHECK_SMALL( v_petsc1->sum() - (1./8.)*nDofVh1, 1e-9 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), (1./8.)*nLocalDofWithGhostVh1, tolCheck );
    // operator=
    *v_petsc1 = *v_petsc2;
    BOOST_CHECK( v_petsc1->sum() == 5*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), 5*nLocalDofWithGhostVh1, tolCheck );
    // operator= (with ublas vector)
    v_ublas1.setConstant( 4 );
    *v_petsc1 = v_ublas1;
    BOOST_CHECK( v_petsc1->sum() == 4*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), 4*nLocalDofWithGhostVh1, tolCheck );
    vB1_ublas1.setConstant( 9 );
    *vB1_petsc1 = vB1_ublas1;
    BOOST_CHECK( vB1_petsc1->sum() == 9*nDofVhB1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), 9*nLocalDofWithGhostVhB1, tolCheck );

    // l1Norm, l2Norm
    v_petsc1->setConstant( -3 );
    BOOST_CHECK( v_petsc1->l1Norm() == 3*nDofVh1 );
    BOOST_CHECK_SMALL( v_petsc1->l2Norm() - std::sqrt(9*nDofVh1), 1e-9 );

    // linftyNorm, min, max
    rank_type myrank = Environment::rank();
    rank_type worldsize = Environment::numberOfProcessors();
    v_petsc1->setConstant( 2 + myrank );
    v_petsc1->close();//v_petsc1->localize();
    BOOST_CHECK( v_petsc1->linftyNorm() == (2+worldsize-1) );
    BOOST_CHECK( v_petsc1->max() == (2+worldsize-1) );
    BOOST_CHECK( v_petsc1->min() == 2 );
    // pointwiseMult, pointwiseDivide
    v_petsc1->setConstant( 4 );
    v_petsc2->setConstant( 5 );
    auto v_petsc3 = backendPetsc->newVector( Vh1 );
    v_petsc3->pointwiseMult( *v_petsc1,*v_petsc2 );
    BOOST_CHECK_CLOSE( v_petsc3->sum(), (4*5)*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc3), (4*5)*nLocalDofWithGhostVh1, tolCheck );
    v_petsc3->pointwiseDivide( *v_petsc1,*v_petsc2 );
    BOOST_CHECK_CLOSE( v_petsc3->sum(), (4./5.)*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc3), (4./5.)*nLocalDofWithGhostVh1, tolCheck );
    // pointwiseMult, pointwiseDivide (with ublas vector)
    v_petsc1->setConstant( 6 );
    v_ublas1.setConstant( 7 );
    v_petsc3->pointwiseMult( *v_petsc1,v_ublas1 );
    BOOST_CHECK_CLOSE( v_petsc3->sum(), (6*7)*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc3), (6*7)*nLocalDofWithGhostVh1, tolCheck );
    v_petsc3->pointwiseDivide( *v_petsc1,v_ublas1 );
    BOOST_CHECK_CLOSE( v_petsc3->sum(), (6./7.)*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc3), (6./7.)*nLocalDofWithGhostVh1, tolCheck );
    // pointwiseMult, pointwiseDivide (with ublas range vector)
    vB1_petsc2->setConstant( 6 );
    vB1_ublas1.setConstant( 4 );
    vB1_petsc1->pointwiseMult( *vB1_petsc2,vB1_ublas1 );
    BOOST_CHECK_CLOSE( vB1_petsc1->sum(), (6*4)*nDofVhB1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (6*4)*nLocalDofWithGhostVhB1, tolCheck );
    vB1_petsc1->pointwiseDivide( *vB1_petsc2,vB1_ublas1 );
    BOOST_CHECK_CLOSE( vB1_petsc1->sum(), (6./4.)*nDofVhB1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (6./4.)*nLocalDofWithGhostVhB1, tolCheck );

    auto mat = backendPetsc->newMatrix( _test=Vh1,_trial=Vh1 );
    for ( size_type k=0;k<Vh1->nLocalDof();++k )
        mat->set(k,k,3.);
    mat->close();
    auto matB1 = backendPetsc->newMatrix( _test=VhB1,_trial=VhB1 );
    for ( size_type k=0;k<VhB1->nLocalDof();++k )
        matB1->set(k,k,3.);
    matB1->close();

    // addVector with matrix
    v_petsc1->addVector( v_petsc2, mat );
    BOOST_CHECK( v_petsc1->sum() == (6+3*5)*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), (6+3*5)*nLocalDofWithGhostVh1, tolCheck );
    // addVector with matrix (with ublas vector)
    v_petsc1->setConstant( 5. );
    v_ublas1.setConstant( 8. );
    v_petsc1->addVector( v_ublas1, *mat );
    BOOST_CHECK( v_petsc1->sum() == (5+3*8)*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), (5+3*8)*nLocalDofWithGhostVh1, tolCheck );
    // addVector with matrix (with ublas range vector)
    vB1_petsc1->setConstant( 5. );
    vB1_ublas1.setConstant( 8. );
    vB1_petsc1->addVector( vB1_ublas1, *matB1 );
    BOOST_CHECK( v_petsc1->sum() == (5+3*8)*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), (5+3*8)*nLocalDofWithGhostVh1, tolCheck );


    // clone
    auto veccloned = v_petsc1->clone();
    veccloned->setConstant( 3. );
    BOOST_CHECK_SMALL( veccloned->sum() - 3*nDofVh1, 1e-9 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*veccloned), 3*nLocalDofWithGhostVh1, tolCheck );
}
BOOST_AUTO_TEST_CASE( test_vector_petsc_extarray )
{
    BOOST_TEST_MESSAGE( "test_vector_petsc_extarray" );
    using namespace Feel;
    double tolCheck = 1e-9;
    auto mesh = unitSquare();
    auto Vh1 = Pchv<2>( mesh );
    size_type nDofVh1 = Vh1->nDof();
    size_type nLocalDofWithGhostVh1 = Vh1->nLocalDofWithGhost();
    auto backendPetsc = backend(_kind="petsc");

    auto v_ublas1 = Vh1->element();
    auto v_ublas2 = Vh1->element();
    auto v_petsc1 = toPETScPtr( v_ublas1 );
    auto v_petsc2 = toPETScPtr( v_ublas2 );
    auto v_petsc3 = backendPetsc->newVector( Vh1 );

    auto VhB = THch<2>( mesh );
    auto VhB1 = VhB->functionSpace<1>();
    size_type nDofVhB = VhB->nDof();
    size_type nDofVhB0 = VhB->functionSpace<0>()->nDof();
    size_type nDofVhB1 = VhB->functionSpace<1>()->nDof();
    size_type nLocalDofWithGhostVhB0 = VhB->functionSpace<0>()->nLocalDofWithGhost();
    size_type nLocalDofWithGhostVhB1 = VhB->functionSpace<1>()->nLocalDofWithGhost();

    auto vB_ublas1 = VhB->element();
    auto vB0_ublas1 = vB_ublas1.element<0>();
    auto vB1_ublas1 = vB_ublas1.element<1>();
    auto vB_ublas2 = VhB->element();
    auto vB0_ublas2 = vB_ublas2.element<0>();
    auto vB1_ublas2 = vB_ublas2.element<1>();
    auto vB1_ublas3 = VhB1->element();

    auto vB_petsc1 = toPETScPtr( vB_ublas1 );
    auto vB0_petsc1 = toPETScPtr( vB0_ublas1 );
    auto vB1_petsc1 = toPETScPtr( vB1_ublas1 );
    auto vB_petsc2 = backendPetsc->newVector( VhB );
    auto vB1_petsc2 = toPETScPtr( vB1_ublas2 );
    auto vB1_petsc3 = toPETScPtr( vB1_ublas3 );

    // sum
    v_ublas1.setConstant( 2 );
    BOOST_CHECK( v_petsc1->sum() == 2*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), 2*nLocalDofWithGhostVh1, tolCheck );
    v_petsc1->setConstant( 3 );
    BOOST_CHECK( v_petsc1->sum() == 3*nDofVh1 );
    BOOST_CHECK( v_ublas1.sum() == 3*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), 3*nLocalDofWithGhostVh1, tolCheck );
    v_ublas1.setConstant( 2 );
    BOOST_CHECK( v_petsc1->sum() == 2*nDofVh1 );
    BOOST_CHECK( v_ublas1.sum() == 2*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), 2*nLocalDofWithGhostVh1, tolCheck );
    // sum (case two arrays : active and ghost)
    vB_petsc1->setConstant( 3 );
    BOOST_CHECK( vB_petsc1->sum() == 3*nDofVhB );
    BOOST_CHECK( vB_ublas1.sum() == 3*nDofVhB );
    vB0_ublas1.setConstant( 6 );
    BOOST_CHECK( vB0_petsc1->sum() == 6*nDofVhB0 );
    BOOST_CHECK( vB0_ublas1.sum() == 6*nDofVhB0 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB0_petsc1), 6*nLocalDofWithGhostVhB0, tolCheck );

    // add scalar
    v_petsc1->add( 1. );
    BOOST_CHECK( v_petsc1->sum() == 3*nDofVh1 );
    BOOST_CHECK( v_ublas1.sum() == 3*nDofVh1 );
    // add scalar (case two arrays : active and ghost)
    vB_petsc1->add( 1. );
    BOOST_CHECK( vB0_petsc1->sum() == ((6+1)*nDofVhB0) );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB0_petsc1), (6+1)*nLocalDofWithGhostVhB0, tolCheck );
    BOOST_CHECK( vB1_petsc1->sum() == ((3+1)*nDofVhB1) );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (3+1)*nLocalDofWithGhostVhB1, tolCheck );
    BOOST_CHECK( vB_petsc1->sum() == ( ((6+1)*nDofVhB0) + ((3+1)*nDofVhB1) ) );
    BOOST_CHECK( vB_ublas1.sum() == ( ((6+1)*nDofVhB0) + ((3+1)*nDofVhB1) ) );

    // add vector
    v_ublas2.setConstant( 7 );
    v_petsc1->add( 1., *v_petsc2 );
    BOOST_CHECK( v_petsc1->sum() == (3+7)*nDofVh1 );
    BOOST_CHECK( v_ublas1.sum() == (3+7)*nDofVh1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*v_petsc1), (3+7)*nLocalDofWithGhostVh1, tolCheck );
    v_petsc3->setConstant( 5 );
    v_petsc1->add( 1., *v_petsc3 );
    BOOST_CHECK( v_petsc1->sum() == (3+7+5)*nDofVh1 );
    BOOST_CHECK( v_ublas1.sum() == (3+7+5)*nDofVh1 );
    v_petsc3->add( 2., *v_petsc1 );
    BOOST_CHECK( v_petsc3->sum() == (5+2*(3+7+5))*nDofVh1 );
    // add vector (case two arrays : active and ghost)
    vB1_petsc1->setConstant( 3 );
    vB1_petsc2->setConstant( 7 );
    vB1_petsc1->add( 1., *vB1_petsc2 );
    BOOST_CHECK( vB1_petsc1->sum() == (3+7)*nDofVhB1 );
    BOOST_CHECK( vB1_ublas1.sum() == (3+7)*nDofVhB1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (3+7)*nLocalDofWithGhostVhB1, tolCheck );
    vB1_ublas3.setConstant( 5 );
    vB1_petsc1->add( 1., *vB1_petsc3 );
    BOOST_CHECK( vB1_petsc1->sum() == (3+7+5)*nDofVhB1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (3+7+5)*nLocalDofWithGhostVhB1, tolCheck );
    vB1_petsc3->add( 2., *vB1_petsc1 );
    BOOST_CHECK( vB1_ublas3.sum() == (5+2*(3+7+5))*nDofVhB1 );
    // add vector (case two arrays : active and ghost + ublas vector)
    vB1_ublas1.setConstant( 8 );
    vB1_ublas2.setConstant( 9 );
    vB1_ublas3.setConstant( 7 );
    vB1_petsc1->add( -4., vB1_ublas2 );
    BOOST_CHECK_CLOSE( vB1_petsc1->sum(), (8.-4.*9)*nDofVhB1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (8.-4.*9.)*nLocalDofWithGhostVhB1, tolCheck );
    vB1_petsc1->add( -2., vB1_ublas3 );
    BOOST_CHECK_CLOSE( vB1_petsc1->sum(), ((8.-4.*9)-2*7)*nDofVhB1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), ((8.-4.*9.)-2*7)*nLocalDofWithGhostVhB1, tolCheck );

    // dot
    v_ublas1.setConstant( 3 );
    v_ublas2.setConstant( 7 );
    v_petsc3->setConstant( 5 );
    BOOST_CHECK( v_petsc1->dot( *v_petsc2 ) == (3*7)*nDofVh1 );
    BOOST_CHECK( v_petsc1->dot( *v_petsc3 ) == (3*5)*nDofVh1 );
    BOOST_CHECK( v_petsc3->dot( *v_petsc1 ) == (3*5)*nDofVh1 );
    // dot (case two arrays : active and ghost)
    vB1_petsc1->setConstant( 3 );
    vB1_petsc2->setConstant( 7 );
    vB1_ublas3.setConstant( 5 );
    BOOST_CHECK( vB1_petsc1->dot( *vB1_petsc2 ) == (3*7)*nDofVhB1 );
    BOOST_CHECK( vB1_petsc1->dot( *vB1_petsc3 ) == (3*5)*nDofVhB1 );
    BOOST_CHECK( vB1_petsc3->dot( *vB1_petsc1 ) == (3*5)*nDofVhB1 );
    // dot (case two arrays : active and ghost + ublas vector)
    vB1_ublas2.setConstant( 5 );
    vB1_ublas3.setConstant( 6 );
    BOOST_CHECK( vB1_petsc1->dot( vB1_ublas2 ) == (3*5)*nDofVhB1 );
    BOOST_CHECK( vB1_petsc1->dot( vB1_ublas3 ) == (3*6)*nDofVhB1 );

    // reciprocal
    v_petsc1->reciprocal();
    BOOST_CHECK_CLOSE( v_petsc1->sum(), (1./3.)*nDofVh1, tolCheck );
    // reciprocal (case two arrays : active and ghost)
    vB1_petsc1->reciprocal();
    BOOST_CHECK_CLOSE( vB1_petsc1->sum(), (1./3.)*nDofVhB1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (1./3.)*nLocalDofWithGhostVhB1, tolCheck );

    // operator=
    (*v_petsc1) = (*v_petsc2);
    BOOST_CHECK( v_petsc1->sum() == 7*nDofVh1 );
    (*v_petsc1) = (*v_petsc3);
    BOOST_CHECK( v_petsc1->sum() == 5*nDofVh1 );
    v_petsc1->setConstant( 3 );
    (*v_petsc3) = (*v_petsc1);
    BOOST_CHECK( v_petsc3->sum() == 3*nDofVh1 );
    // operator= (case two arrays : active and ghost)
    *vB1_petsc1 = *vB1_petsc2;
    BOOST_CHECK( vB1_petsc1->sum() == 5*nDofVhB1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), 5*nLocalDofWithGhostVhB1, tolCheck );
    *vB1_petsc1 = *vB1_petsc3;
    BOOST_CHECK( vB1_petsc1->sum() == 6*nDofVhB1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), 6*nLocalDofWithGhostVhB1, tolCheck );
    vB1_petsc1->setConstant( 3 );
    *vB1_petsc3 = *vB1_petsc1;
    BOOST_CHECK( vB1_petsc3->sum() == 3*nDofVhB1 );
    // operator= (case two arrays : active and ghost + ublas vector)
    vB1_ublas1.setConstant( 5 );
    vB1_ublas3.setConstant( 6 );
    *vB1_petsc1 = vB1_ublas1;
    BOOST_CHECK( vB1_petsc1->sum() == 5*nDofVhB1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), 5*nLocalDofWithGhostVhB1, tolCheck );
    *vB1_petsc1 = vB1_ublas3;
    BOOST_CHECK( vB1_petsc1->sum() == 6*nDofVhB1 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), 6*nLocalDofWithGhostVhB1, tolCheck );

    // l1Norm, l2Norm
    v_petsc1->setConstant( -4 );
    BOOST_CHECK( v_petsc1->l1Norm() == 4*nDofVh1 );
    BOOST_CHECK( v_ublas1.l1Norm() == 4*nDofVh1 );
    BOOST_CHECK_SMALL( v_petsc1->l2Norm() - std::sqrt(4*4*nDofVh1), 1e-9 );
    BOOST_CHECK_SMALL( v_ublas1.l2Norm() - std::sqrt(4*4*nDofVh1), 1e-9 );
    // l1Norm, l2Norm (case two arrays : active and ghost)
    vB1_petsc1->setConstant( -4 );
    BOOST_CHECK( vB1_petsc1->l1Norm() == 4*nDofVhB1 );
    BOOST_CHECK( vB1_ublas1.l1Norm() == 4*nDofVhB1 );
    BOOST_CHECK_SMALL( vB1_petsc1->l2Norm() - std::sqrt(4*4*nDofVhB1), 1e-9 );
    BOOST_CHECK_SMALL( vB1_ublas1.l2Norm() - std::sqrt(4*4*nDofVhB1), 1e-9 );

    // linftyNorm, min, max
    rank_type myrank = Environment::rank();
    rank_type worldsize = Environment::numberOfProcessors();
    v_petsc1->setConstant( 2 + myrank );
    sync(v_ublas1);
    //v_petsc1->close();
    BOOST_CHECK( v_petsc1->linftyNorm() == (2+worldsize-1) );
    BOOST_CHECK( v_petsc1->max() == (2+worldsize-1) );
    BOOST_CHECK( v_petsc1->min() == 2 );
    BOOST_CHECK( v_ublas1.linftyNorm() == (2+worldsize-1) );
    BOOST_CHECK( v_ublas1.max() == (2+worldsize-1) );
    BOOST_CHECK( v_ublas1.min() == 2 );
    // linftyNorm, min, max (case two arrays : active and ghost)
    vB1_petsc1->setConstant( 2 + myrank );
    //sync(vB1_ublas1);//bug
    //sync(vB_ublas1);//bug
    vB1_petsc1->close();//work
    BOOST_CHECK( vB1_petsc1->linftyNorm() == (2+worldsize-1) );
    BOOST_CHECK( vB1_petsc1->max() == (2+worldsize-1) );
    BOOST_CHECK( vB1_petsc1->min() == 2 );
    BOOST_CHECK( vB1_ublas1.linftyNorm() == (2+worldsize-1) );
    BOOST_CHECK( vB1_ublas1.max() == (2+worldsize-1) );
    BOOST_CHECK( vB1_ublas1.min() == 2 );

    // pointwiseMult, pointwiseDivide
    v_petsc1->setConstant( 4 );
    v_petsc2->setConstant( 5 );
    v_petsc3->pointwiseMult( *v_petsc1,*v_petsc2 );
    BOOST_CHECK_SMALL( v_petsc3->sum() - 20*nDofVh1, 1e-9 );
    v_petsc3->pointwiseDivide( *v_petsc1,*v_petsc2 );
    BOOST_CHECK_SMALL( v_petsc3->sum() - (4./5.)*nDofVh1, 1e-9 );
    // pointwiseMult, pointwiseDivide (case two arrays : active and ghost)
    vB1_petsc1->setConstant( 4 );
    vB1_petsc2->setConstant( 5 );
    vB1_petsc3->pointwiseMult( *vB1_petsc1,*vB1_petsc2 );
    BOOST_CHECK_SMALL( vB1_petsc3->sum() - (4*5)*nDofVhB1, 1e-9 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc3), (4*5)*nLocalDofWithGhostVhB1, tolCheck );
    vB1_petsc3->pointwiseDivide( *vB1_petsc1,*vB1_petsc2 );
    BOOST_CHECK_SMALL( vB1_petsc3->sum() - (4./5.)*nDofVhB1, 1e-9 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc3), (4./5.)*nLocalDofWithGhostVhB1, tolCheck );
    vB1_petsc3->setConstant( 6 );
    vB1_petsc1->pointwiseMult( *vB1_petsc2,*vB1_petsc3 );
    BOOST_CHECK_SMALL( vB1_petsc1->sum() - (5*6)*nDofVhB1, 1e-9 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (5*6)*nLocalDofWithGhostVhB1, tolCheck );
    vB1_petsc1->pointwiseDivide( *vB1_petsc2,*vB1_petsc3 );
    BOOST_CHECK_SMALL( vB1_petsc1->sum() - (5./6.)*nDofVhB1, 1e-9 );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vB1_petsc1), (5./6.)*nLocalDofWithGhostVhB1, tolCheck );

    // clone
    auto vecCloned = vB1_petsc1->clone();
    vecCloned->setConstant( 3. );
    BOOST_CHECK_SMALL( vecCloned->sum() - 3*nDofVhB1, 1e-9 );
}

BOOST_AUTO_TEST_CASE( test_mdot )
{
    using namespace Feel;
    auto mesh = unitSquare();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto us = Vh->newElements(3);

    u.setConstant( 1 );
    for( auto w : us )
        u.setConstant( 1 );

    auto backendPetsc = backend( _kind = "petsc" );
    auto v = backendPetsc->newVector( Vh );
    auto vs = backendPetsc->newVectors( Vh, 3 );
    BOOST_CHECK_EQUAL( vs.size(), 3 );

    v->setConstant( 1 );
    for( auto w : vs )
        w->setConstant( 1 );
    
    tic();
    auto r = v->mDot(vs);
    toc("mdot 3 vectors");
    BOOST_TEST_MESSAGE( fmt::format("r.size()={}, r={}", r.size(),r) );
    for( auto [i,w] : enumerate(vs) )
        BOOST_CHECK_EQUAL( r[i], w->sum() );

    r = v->mDot( vs, 2 );
    BOOST_TEST_MESSAGE( fmt::format("r.size()={}, r={}", r.size(),r) );
    for( auto [i,w] : enumerate(vs | ranges::views::take(2)) )
        BOOST_CHECK_EQUAL( r[i], w->sum() );

    r = u.mDot(vs);
    toc("mdot 3 vectors petsc -> ublas");
    BOOST_TEST_MESSAGE( fmt::format("r.size()={}, r={}", r.size(),r) );
    for( auto [i,w] : enumerate(vs) )
        BOOST_CHECK_EQUAL( r[i], w->sum() );

    r = u.mDot(us);
    toc("mdot 3 vectors ublas -> ublas");
    BOOST_TEST_MESSAGE( fmt::format("r.size()={}, r={}", r.size(),r) );
    for( auto [i,w] : enumerate(us) )
        BOOST_CHECK_EQUAL( r[i], w->sum() );

    r = v->mDot(us);
    toc("mdot 3 vectors ublas -> petsc");
    BOOST_TEST_MESSAGE( fmt::format("r.size()={}, r={}", r.size(),r) );
    for( auto [i,w] : enumerate(us) )
        BOOST_CHECK_EQUAL( r[i], w->sum() );



}

BOOST_AUTO_TEST_CASE( test_maxpy )
{
    using namespace Feel;
    auto mesh = unitSquare();
    auto Vh = Pch<1>( mesh );

    auto backendPetsc = backend( _kind = "petsc" );
    auto v = backendPetsc->newVector( Vh );
    auto u = Vh->element();
    auto vs = backendPetsc->newVectors( Vh, 3 );
    BOOST_CHECK_EQUAL( vs.size(), 3 );
    auto us = Vh->newElements(3);
    BOOST_CHECK_EQUAL( us.size(), 3 );

    v->setConstant( 0 );
    for( auto w : vs )
        w->setConstant( 1 );
    for( auto w : us )
        w->setConstant( 1 );
    
    eigen_vector_type<> alpha(3);
    alpha << 1, 2, 3;
    tic();
    v->add(alpha,vs);
    toc("add 3 vectors");
    BOOST_CHECK_EQUAL( v->sum(), 6*vs[0]->size() );

    v->setConstant( 0 );
    tic();
    v->add(alpha,us);
    toc("add 3 vectors ublas -> petsc");
    BOOST_CHECK_EQUAL( v->sum(), 6*us[0]->size() );

    u.setConstant( 0 );
    tic();
    u.add(alpha,vs);
    toc("add 3 vectors petsc -> ublas");
    BOOST_CHECK_EQUAL( u.sum(), 6*vs[0]->size() );

    u.setConstant( 0 );
    tic();
    u.add(alpha,us);
    toc("add 3 vectors ublas -> ublas ");
    BOOST_CHECK_EQUAL( u.sum(), 6*us[0]->size() );
    
    auto sv = vs| ranges::views::take(2);
    v->setConstant( 0 );   
    tic();
    v->add(alpha,{sv.begin(), sv.end()});
    toc("add 2 vectors");
    BOOST_CHECK_EQUAL( v->sum(), 3 * vs[0]->size() );

    v->setConstant( 0 );    
    tic();
    v->add(alpha,vs,2);
    toc("add 2 vectors with stride");
    BOOST_CHECK_EQUAL( v->sum(), 3 * vs[0]->size() );
}
BOOST_AUTO_TEST_SUITE_END()