/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-05-08

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file test_rt.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-05-08
 */
#define BOOST_TEST_MODULE Raviar-Thomas polynomials test
// Boost.Test
#define USE_TEST 1
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <life/lifepoly/raviartthomas.hpp>

BOOST_AUTO_TEST_SUITE( rt_testsuite )

#if 0
BOOST_AUTO_TEST_CASE( rt0 )
{
    using namespace Life;
    typedef RaviartThomas<0>::apply<2>::type rt0_type;
    rt0_type rt0;

    rt0_type::points_type pts(2,1);
    pts( 0, 0 ) = -1./3.; pts( 1, 0 ) = -1./3.;

    //pts = rt0.referenceConvex().barycenterFaces();

    std::cout << "pts= " << pts << "\n";
    auto eval_at_pts = rt0.evaluate( pts );
    std::cout << "eval at pts= " << eval_at_pts << "\n";
}
#endif
BOOST_AUTO_TEST_CASE( rt0_2 )
{
    BOOST_TEST_MESSAGE( "creating RT0 on reference element :" );
    using namespace Life;
    typedef RaviartThomas<0>::apply<2>::type rt0_type;
    typedef boost::shared_ptr<rt0_type> rt0_ptrtype;
    rt0_ptrtype rt0( new rt0_type ) ;

    rt0_type::points_type pts(2,3);
    pts = rt0->referenceConvex().vertices();

    std::cout << "pts= " << pts << "\n";
    auto eval_at_pts = rt0->evaluate( pts );
    std::cout << "eval at pts= " << eval_at_pts << "\n";

    BOOST_TEST_MESSAGE( "creating another element :" );
    typedef GeoND<2,Simplex<2, 1, 2> >::point_type point_type;
    // interval
    typedef GeoND<2,Simplex<2, 1, 2> > tria_type;
    tria_type tria;
    point_type V1; V1( 0 )=1;V1( 1 )=0;
    point_type V2; V2( 0 )=3;V2( 1 )=0;
    point_type V3; V3( 0 )=2;V3( 1 )=1;
    tria.setPoint( 0, V1 );
    tria.setPoint( 1, V2 );
    tria.setPoint( 2, V3 );
    ublas::vector<double> G1( 2 );
    G1(0)=2;G1(1)=1./3.;
    tria.update();
    BOOST_TEST_MESSAGE( "evaluating the RT0 polynomialset on real element :" );

    auto bpts = tria.gm()->referenceConvex().barycenterFaces();
    auto gmpc = tria.gm()->preCompute( tria.gm(), bpts );
    auto rtpc = rt0->preCompute( rt0, bpts );

    auto geoctx = tria.gm()->context<vm::GRAD>( tria, gmpc );
#if 1
    //auto rtctx = rt0->ctx<vm::GRAD,rt0_type, decltype(*tria.gm()),decltype(*rtpc), tria_type>( rt0, tria.gm(), rtpc );
    //auto rtctx = rt0->ctx( rt0, tria.gm() );
    //auto gmctx = tria.gm()->context<vm::GRAD>( tria, gmpc );
    //auto rtctx = rt0->ctx<vm::GRAD,rt0_type,decltype(*geoctx),decltype(*rtpc),tria_type>( rt0, geoctx, rtpc, tria );
    //auto rtctx = rt0->ctx();
    //auto rtctx = rt0->context<vm::GRAD>( rt0, tria.gm(), rtpc );

    typedef rt0_type::Context<vm::GRAD,rt0_type,tria_type::gm_type,tria_type> rtctx_type;
    boost::shared_ptr<rtctx_type> rtctx( new rtctx_type( rt0, geoctx, rtpc ) );

#endif

}

BOOST_AUTO_TEST_SUITE_END()


