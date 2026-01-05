/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-17

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
     \file test_poly.cpp
     \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
     \date 2005-08-17
 */

#define BOOST_TEST_MODULE test_poly
#include <boost/test/data/test_case.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/polynomial.hpp>
#include <feel/feelpoly/operations.hpp>
#include <cmath>
#include <vector>

namespace bdata = boost::unit_test::data;
namespace ublas = boost::numeric::ublas;
using namespace Feel;
struct x2
{
    typedef double value_type;
    ublas::vector<double>
    operator()( ublas::matrix<double,ublas::column_major> const& m ) const
    {
        ublas::vector<double> v( m.size2() );
        v.clear();

        for ( int i = 0; i < m.size2(); ++i )
        {
            for ( int j = 0; j < m.size1(); ++j )
                v( i ) += ( j+2 )*m( j,i )*m( j,i );

            if ( m.size1() == 2 )
                v( i ) += 100*m( 0,i )*m( 1,i )+200*m( 1,i )*m( 0,i );

            //v( i ) = m(0,i)*m(0,i);
        }

        return v;
    }
};

template<int Dim>
ublas::matrix<double,ublas::column_major> make_points();

template<>
ublas::matrix<double,ublas::column_major> make_points<1>()
{
    ublas::matrix<double,ublas::column_major> pts( 1, 3 );
    pts( 0,0 ) = 0;
    pts( 0,1 ) = -1;
    pts( 0,2 ) = 1;
    return pts;
}

template<>
ublas::matrix<double,ublas::column_major> make_points<2>()
{
    ublas::matrix<double,ublas::column_major> pts( 2, 4 );
    pts( 0,0 ) = -1; pts( 1,0 ) = -1;
    pts( 0,1 ) =  1; pts( 1,1 ) = -1;
    pts( 0,2 ) = -1; pts( 1,2 ) =  1;
    pts( 0,3 ) =  0; pts( 1,3 ) =  0;
    return pts;
}

template<>
ublas::matrix<double,ublas::column_major> make_points<3>()
{
    ublas::matrix<double,ublas::column_major> pts( 3, 4 );
    pts( 0,0 ) = -1; pts( 1,0 ) = -1; pts( 2,0 ) = -1;
    pts( 0,1 ) =  1; pts( 1,1 ) = -1; pts( 2,1 ) = -1;
    pts( 0,2 ) = -1; pts( 1,2 ) =  1; pts( 2,2 ) = -1;
    pts( 0,3 ) =  0; pts( 1,3 ) =  0; pts( 2,3 ) =  0;
    return pts;
}

template<int Dim, int QuadOrder, int PolyOrder>
void run_polynomial_case()
{
    using polyset_type = Feel::detail::OrthonormalPolynomialSet<Dim,Dim,PolyOrder,Scalar>;

    auto pts = make_points<Dim>();
    polyset_type ps;
    Polynomial<polyset_type, Scalar> p1 = project( ps, x2(), IM<Dim,QuadOrder>() );
    PolynomialSet<polyset_type, Scalar> p( ps );
    p.insert( p1.toSet( true ) );
    p.insert( p1.toSet( true ) );

    const auto grad_vals = p.gradient().evaluate( pts );
    const auto hess_vals = p.gradient().gradient().evaluate( pts );

    BOOST_TEST_CONTEXT( "dim=" << Dim << " quad=" << QuadOrder )
    {
        BOOST_CHECK_GT( p.coeff().size1(), 0 );
        BOOST_CHECK_EQUAL( grad_vals.size2(), pts.size2() );
        BOOST_CHECK_EQUAL( hess_vals.size2(), pts.size2() );
        BOOST_CHECK( std::isfinite( grad_vals( 0, 0 ) ) );
        BOOST_CHECK( std::isfinite( hess_vals( 0, 0 ) ) );
    }
}

using Runner = void(*)();
static const std::vector<Runner> kPolyRunners = {
    &run_polynomial_case<1,3,2>,
    &run_polynomial_case<2,3,2>,
    &run_polynomial_case<3,5,3>
};

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( polynomial_suite )

BOOST_DATA_TEST_CASE( polynomial_variants, bdata::make( kPolyRunners ), run_case )
{
    run_case();
}

BOOST_AUTO_TEST_SUITE_END()
