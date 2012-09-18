/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-12-11

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file test_moment.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-12-11
 */
// Boost.Test
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
using boost::unit_test::test_suite;

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <fstream>
#include <feel/feelpoly/moment.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>

typedef boost::mpl::list<boost::mpl::int_<2> > test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_QK, T, test_types )
{
    using namespace Feel;

    Moment<2,T::value,Hypercube<2> > m;
#if 0
    std::cout << "1  :" << m.template pick<Scalar>( 0 ).evaluate( m.points() ) << "\n";
    std::cout << "x  :" << m.template pick<Scalar>( 1 ).evaluate( m.points() ) << "\n";
    std::cout << "x^2:" << m.template pick<Scalar>( 2 ).evaluate( m.points() ) << "\n";
    std::cout << "y  :" << m.template pick<Scalar>( 3 ).evaluate( m.points() ) << "\n";
    std::cout << "y^2:" << m.template pick<Scalar>( 6 ).evaluate( m.points() ) << "\n";

    auto p  = m.template pick<Scalar>( 2 ) - m.template pick<Scalar>( 6 );
    std::cout << "x^2-y^2:" << p.evaluate( m.points() ) << "\n";

    std::cout << "d 1/dx  :" << m.template pick<Scalar>( 0 ).derivate( 0,m.points() ) << "\n";
    std::cout << "d x/dx  :" << m.template pick<Scalar>( 1 ).derivate( 0,m.points() ) << "\n";
    std::cout << "d x^2/dx:" << m.template pick<Scalar>( 2 ).derivate( 0,m.points() ) << "\n";
    std::cout << "d y/dx   :" << m.template pick<Scalar>( 3 ).derivate( 0,m.points() ) << "\n";
    std::cout << "d y^2 /dx:" << m.template pick<Scalar>( 6 ).derivate( 0,m.points() ) << "\n";
    std::cout << "dx^2-y^2/x :" << p.derivate( 0, m.points() ) << "\n";
    std::cout << "\n";
    std::cout << "d 1/dy  :" << m.template pick<Scalar>( 0 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d x/dy  :" << m.template pick<Scalar>( 1 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d x^2/dy:" << m.template pick<Scalar>( 2 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d y/dy   :" << m.template pick<Scalar>( 3 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d y^2 /dy:" << m.template pick<Scalar>( 6 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d x^2-y^2 /dy:" << p.derivate( 1,m.points() ) << "\n";

    PolynomialSet<Moment<2,T::value,Hypercube<2> >, Scalar> pset( m );
    pset.insert( m.template pick<Scalar>( 0 ).toSet(), true );
    pset.insert( m.template pick<Scalar>( 1 ).toSet() );
    pset.insert( m.template pick<Scalar>( 3 ).toSet() );
    pset.insert( p.toSet() );

    std::cout << "pset :" << pset.evaluate( m.points() ) << "\n";
    std::cout << "d pset/dx :" << pset.derivate( 0, m.points() ) << "\n";
    std::cout << "d pset/dy :" << pset.derivate( 1, m.points() ) << "\n";

    fem::detail::RannacherTurekPolynomialSet<2,Scalar> RQ;
    std::cout << "rq :" << RQ.evaluate( m.points() ) << "\n";
    std::cout << "d rq/dx :" << RQ.derivate( 0, m.points() ) << "\n";
    std::cout << "d rq/dy :" << RQ.derivate( 1, m.points() ) << "\n";

    fem::CrouzeixRaviart<2,2,Scalar,double,Hypercube> cr;
    std::cout << "cr :" << cr.evaluate( m.points() ) << "\n";

    fem::detail::RannacherTurekPolynomialSet<2,Vectorial> RQv;
    std::cout << "rqv :" << RQv.evaluate( m.points() ) << "\n";
    std::cout << "d rqv/dx :" << RQv.derivate( 0, m.points() ) << "\n";
    std::cout << "d rqv/dy :" << RQv.derivate( 1, m.points() ) << "\n";
#endif
    fem::CrouzeixRaviart<2,2,Scalar,double,Hypercube> cr;
    std::cout << "cr :" << cr.evaluate( m.points() ) << "\n";
    std::cout << "d cr/dx :" << cr.derivate( 0, m.points() ) << "\n";
    std::cout << "d cr/dy :" << cr.derivate( 1, m.points() ) << "\n";

    fem::CrouzeixRaviart<2,2,Vectorial,double,Hypercube> crv;
    std::cout << "crv :" << crv.evaluate( m.points() ) << "\n";
    std::cout << "d crv/dx :" << crv.derivate( 0, m.points() ) << "\n";
    std::cout << "d crv/dy :" << crv.derivate( 1, m.points() ) << "\n";

}
