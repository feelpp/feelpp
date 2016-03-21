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
   \date 2010-08-07
 */
#include <cmath>

#include <boost/timer.hpp>
#define BOOST_TEST_MODULE vector testsuite
#include <testsuite/testsuite.hpp>

#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/thch.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( vector )

BOOST_AUTO_TEST_CASE( test_vector_ublas_base )
{
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
    using namespace Feel;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
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

    size_type nDofVh1 = Vh1->nDof();
    size_type nDofVh2 = Vh2->nDof();
    size_type nDofVh2a = v2a.functionSpace()->nDof();
    size_type nDofVh2b = v2b.functionSpace()->nDof();
    size_type nDofVh3 = Vh3->nDof();
    BOOST_CHECK( nDofVh2a == nDofVh3 && v2a.nLocalDof() == v3.nLocalDof() );

    // sum
    BOOST_CHECK( v1.sum() == 2*nDofVh1 );
    BOOST_CHECK( v2a.sum() == 5*nDofVh2a );
    BOOST_CHECK( v2b.sum() == 7*nDofVh2b );
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
#if 0
    auto prod2a = element_product(v2a,v2a);
    BOOST_CHECK( prod2a.sum() == 25*nDofVh2a );
    auto prod2b = element_product(v2b,v2b);
    BOOST_CHECK( prod2b.sum() == 49*nDofVh2b );
    auto prod2a3 = element_product(v2a,v3);
    BOOST_CHECK( prod2a3.sum() == 45*nDofVh2a );
#endif
    // add scalar
    v1.add( 3 );
    BOOST_CHECK( v1.sum() == 5*nDofVh1 );
    v2a.add( 3 );
    BOOST_CHECK( v2a.sum() == 8*nDofVh2a );
    v2b.add( 3 );
    BOOST_CHECK( v2b.sum() == 10*nDofVh2b );
    // operator=
    v1 = v1bis;
    BOOST_CHECK( v1.sum() == 3*nDofVh1 );
    v3 = v2a;
    BOOST_CHECK( v3.sum() == 8*nDofVh2a );
    v3.setConstant( 9 );
    v2a = v3;
    BOOST_CHECK( v3.sum() == 9*nDofVh2a );
    // add vector
    v1.setConstant( 2 );
    v2a.setConstant( 5 );
    v2b.setConstant( 7 );
    v1.add(-1.,v1 );
    BOOST_CHECK( v1.sum() == 0 );
    v2a.add( 1., v3 );
    BOOST_CHECK( v2a.sum() == 14*nDofVh2a );
    v3.add( 1., v2a );
    BOOST_CHECK( v3.sum() == 23*nDofVh3 );
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

}
BOOST_AUTO_TEST_SUITE_END()
