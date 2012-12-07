/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-12

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
   \file test_lag.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-12
 */
#define BOOST_TEST_MODULE lagrange polynomials test
// Boost.Test
#define USE_TEST 1
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <testlagrange.hpp>

#if 1
BOOST_AUTO_TEST_SUITE( lagrange_testsuite )

BOOST_AUTO_TEST_SUITE( lagrange_1d_simplex_testsuite )

BOOST_AUTO_TEST_CASE( lag11scr64s )
{
    TestLagrange<fem::Lagrange<1, 1, 1, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag12scr64s )
{
    TestLagrange<fem::Lagrange<1, 1, 2, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag15scr64s )
{
    TestLagrange<fem::Lagrange<1, 1, 5, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag110scr64s )
{
    TestLagrange<fem::Lagrange<1, 1, 10, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // lagrange_1d_simplex_testsuite


BOOST_AUTO_TEST_SUITE( lagrange_2d_simplex_testsuite )

BOOST_AUTO_TEST_CASE( lag21scr64s )
{
    TestLagrange<fem::Lagrange<2, 2, 1, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag21ver64s )
{
    TestLagrange<fem::Lagrange<2, 2, 1, Vectorial, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag22scr64s )
{
    TestLagrange<fem::Lagrange<2, 2, 2, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag25scr64s )
{
    TestLagrange<fem::Lagrange<2, 2, 5, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag210scr64s )
{
    TestLagrange<fem::Lagrange<2, 2, 10, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // lagrange_2d_simplex_testsuite

BOOST_AUTO_TEST_SUITE( lagrange_3d_simplex_testsuite )

BOOST_AUTO_TEST_CASE( lag31scr64s )
{
    TestLagrange<fem::Lagrange<3, 3, 1, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag32scr64s )
{
    TestLagrange<fem::Lagrange<3, 3, 2, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag35scr64s )
{
    TestLagrange<fem::Lagrange<3, 3, 5, Scalar, Continuous, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag310scr64s )
{
    TestLagrange<fem::Lagrange<3, 3, 10, Scalar, Continuous, real64_type, Simplex> > t( 1e-9 );
    t();
}

BOOST_AUTO_TEST_SUITE_END() // lagrange_3d_simplex_testsuite

BOOST_AUTO_TEST_SUITE( lagrange_1d_hypercube_testsuite )

BOOST_AUTO_TEST_CASE( lag11scr64s )
{
    TestLagrange<fem::Lagrange<1, 1, 1, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag12scr64s )
{
    TestLagrange<fem::Lagrange<1, 1, 2, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag15scr64s )
{
    TestLagrange<fem::Lagrange<1, 1, 5, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag110scr64s )
{
    TestLagrange<fem::Lagrange<1, 1, 10, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // lagrange_1d_hypercube_testsuite


BOOST_AUTO_TEST_SUITE( lagrange_2d_hypercube_testsuite )

BOOST_AUTO_TEST_CASE( lag21scr64s )
{
    TestLagrange<fem::Lagrange<2, 2,1, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag22scr64s )
{
    TestLagrange<fem::Lagrange<2, 2,2, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag25scr64s )
{
    TestLagrange<fem::Lagrange<2, 2,5, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag210scr64s )
{
    TestLagrange<fem::Lagrange<2, 2,10, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // lagrange_2d_hypercube_testsuite

BOOST_AUTO_TEST_SUITE( lagrange_3d_hypercube_testsuite )

BOOST_AUTO_TEST_CASE( lag31scr64s )
{
    TestLagrange<fem::Lagrange<3, 3, 1, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag32scr64s )
{
    TestLagrange<fem::Lagrange<3, 3, 2, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag35scr64s )
{
    TestLagrange<fem::Lagrange<3, 3, 5, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag310scr64s )
{
    TestLagrange<fem::Lagrange<3, 3, 7, Scalar, Continuous, real64_type, Hypercube> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // lagrange_3d_hypercube_testsuite


BOOST_AUTO_TEST_SUITE_END()
#else


#if defined( USE_TEST )
void test_scalar_continuous_simplex( test_suite* test )
{
#if USE_TEST
    // 1D simplex
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<1, 1, 1, Scalar, Continuous, real64_type, Simplex> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<1, 1, 2, Scalar, Continuous, real64_type, Simplex> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<1, 1, 3, Scalar, Continuous, real64_type, Simplex> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<1, 1, 4, Scalar, Continuous, real64_type, Simplex> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<1, 1, 10, Scalar, Continuous, real64_type, Simplex> >() ) ) );

    // 2D simplex
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<2, 2, 1, Scalar, Continuous, real64_type, Simplex> >() )  ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<2, 2, 2, Scalar, Continuous, real64_type, Simplex> >() )  ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<2, 2, 10, Scalar, Continuous, real64_type, Simplex> >() )  ) );

    // 3D simplex
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<3, 3, 1, Scalar, Continuous, real64_type, Simplex> >() )  ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<3, 3, 2, Scalar, Continuous, real64_type, Simplex> >() )  ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<3, 3, 3, Scalar, Continuous, real64_type, Simplex> >() )  ) );
#else

    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<3, 3, 3, Scalar, Continuous, qd_real, Simplex> >() )  ) );

#endif
}

void test_scalar_continuous_simplex_product( test_suite* test )
{
#if USE_TEST
    // 1D simplex
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<1, 1, 1, Scalar, Continuous, real64_type, Hypercube> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<1, 1, 2, Scalar, Continuous, real64_type, Hypercube> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<1, 1, 10, Scalar, Continuous, real64_type, Hypercube> >() ) ) );

    // 2D simplex
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<2, 2, 1, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<2, 2, 2, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<2, 2, 10, Scalar, Continuous, real64_type, Hypercube> >() )  ) );

    // 3D simplex
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<3, 3, 1, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<3, 3, 2, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
    test->add( BOOST_TEST_CASE( ( TestLagrange<fem::Lagrange<3, 3, 3, Scalar, Continuous, real64_type, Hypercube> >() )  ) );


#endif
}

test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
#if defined( FEELPP_HAS_QD_QD_H )
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif

    test_suite* test = BOOST_TEST_SUITE( "Finite element test suite" );


    test_scalar_continuous_simplex( test );
    //test_scalar_continuous_simplex_product( test );



    return test;
}
#else
int main( int argc, char* argv[] )
{
#if 0
    TestLagrange<fem::Lagrange<2, 2, 1, Scalar, Continuous, real64_type, Simplex> >();
    TestLagrange<fem::Lagrange<2, 2, 2, Scalar, Continuous, real64_type, Simplex> >();
    TestLagrange<fem::Lagrange<2, 2, 10, Scalar, Continuous, real64_type, Simplex> >();

    TestLagrange<fem::Lagrange<3, 3, 1, Scalar, Continuous, real64_type, Simplex> >();
    TestLagrange<fem::Lagrange<3, 3, 2, Scalar, Continuous, real64_type, Simplex> >();
    TestLagrange<fem::Lagrange<3, 3, 3, Scalar, Continuous, real64_type, Simplex> >();
#endif

#if 0
    // 1D simplex
    ( ( ( TestLagrange<fem::Lagrange<1, 1, 1, Scalar, Continuous, real64_type, Hypercube> >() ) ) );
    ( ( ( TestLagrange<fem::Lagrange<1, 1, 2, Scalar, Continuous, real64_type, Hypercube> >() ) ) );
    ( ( ( TestLagrange<fem::Lagrange<1, 1,10, Scalar, Continuous, real64_type, Hypercube> >() ) ) );

    // 2D simplex
    ( ( ( TestLagrange<fem::Lagrange<2, 2, 1, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
    ( ( ( TestLagrange<fem::Lagrange<2, 2, 2, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
    ( ( ( TestLagrange<fem::Lagrange<2, 2,10, Scalar, Continuous, real64_type, Hypercube> >() )  ) );

    // 3D simplex
    ( ( ( TestLagrange<fem::Lagrange<3, 3, 1, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
    ( ( ( TestLagrange<fem::Lagrange<3, 3, 2, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
    ( ( ( TestLagrange<fem::Lagrange<3, 3, 3, Scalar, Continuous, real64_type, Hypercube> >() )  ) );
#else
#if defined( FEELPP_HAS_QD_QD_H )

    unsigned int old_cw;
    fpu_fix_start( &old_cw );

    {
        TestLagrange<fem::Lagrange<2, 2, 1, Scalar, Continuous, qd_real, Simplex> > a;
        a();
    }
    {
        TestLagrange<fem::Lagrange<2, 2, 3, Scalar, Continuous, dd_real, Simplex> > a;
        a();
    }
    {
        TestLagrange<fem::Lagrange<3, 3, 1, Scalar, Continuous, qd_real, Simplex> > a;
        a();
    }

    fpu_fix_end( &old_cw );

#endif /* FEELPP_HAS_QD_QD_H  */

#endif
}

#endif
#endif // 0
