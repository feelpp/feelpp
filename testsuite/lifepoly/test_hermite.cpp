/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-03-04

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
   \file testhermite.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-03-04
 */
#define BOOST_TEST_MODULE hermite polynomials test
// Boost.Test
#define USE_TEST 1
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <life/lifepoly/polynomialset.hpp>
#include <life/lifepoly/orthonormalpolynomialset.hpp>
#include <test_hermite.hpp>

#if 1
BOOST_AUTO_TEST_SUITE( hermite_testsuite )

BOOST_AUTO_TEST_SUITE( hermite_simplex_testsuite )

BOOST_AUTO_TEST_CASE( hermite )
{
    TestHermite<fem::Hermite<1, 3, Scalar, real64_type, Simplex> > t;
    t();
    TestHermite<fem::Hermite<2, 3, Scalar, real64_type, Simplex> > t2;
    t2();
    //TestHermite<fem::Hermite<3, 3, Scalar, real64_type, Simplex> > t3;
    //t3();
}

#if 0
BOOST_AUTO_TEST_CASE( lag12scr64s )
{
    TestHermite<fem::Hermite<1, 2, Scalar, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag15scr64s )
{
    TestHermite<fem::Hermite<1, 5, Scalar, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag110scr64s )
{
    TestHermite<fem::Hermite<1, 10, Scalar, real64_type, Simplex> > t;
    t();
}
#endif
BOOST_AUTO_TEST_SUITE_END() // hermite_1d_simplex_testsuite

#if 0

BOOST_AUTO_TEST_SUITE( hermite_2d_simplex_testsuite )

BOOST_AUTO_TEST_CASE( lag21scr64s )
{
    TestHermite<fem::Hermite<2, 1, Scalar, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag21ver64s )
{
    TestHermite<fem::Hermite<2, 1, Vectorial, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag22scr64s )
{
    TestHermite<fem::Hermite<2, 2, Scalar, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag25scr64s )
{
    TestHermite<fem::Hermite<2, 5, Scalar, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag210scr64s )
{
    TestHermite<fem::Hermite<2, 10, Scalar, real64_type, Simplex> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // hermite_2d_simplex_testsuite

BOOST_AUTO_TEST_SUITE( hermite_3d_simplex_testsuite )

BOOST_AUTO_TEST_CASE( lag31scr64s )
{
    TestHermite<fem::Hermite<3, 1, Scalar, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag32scr64s )
{
    TestHermite<fem::Hermite<3, 2, Scalar, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag35scr64s )
{
    TestHermite<fem::Hermite<3, 5, Scalar, real64_type, Simplex> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag310scr64s )
{
    TestHermite<fem::Hermite<3, 10, Scalar, real64_type, Simplex> > t(1e-9);
    t();
}

BOOST_AUTO_TEST_SUITE_END() // hermite_3d_simplex_testsuite

BOOST_AUTO_TEST_SUITE( hermite_1d_simplexproduct_testsuite )

BOOST_AUTO_TEST_CASE( lag11scr64s )
{
    TestHermite<fem::Hermite<1, 1, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag12scr64s )
{
    TestHermite<fem::Hermite<1, 2, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag15scr64s )
{
    TestHermite<fem::Hermite<1, 5, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag110scr64s )
{
    TestHermite<fem::Hermite<1, 10, Scalar, real64_type, SimplexProduct> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // hermite_1d_simplexproduct_testsuite


BOOST_AUTO_TEST_SUITE( hermite_2d_simplexproduct_testsuite )

BOOST_AUTO_TEST_CASE( lag21scr64s )
{
    TestHermite<fem::Hermite<2, 1, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag22scr64s )
{
    TestHermite<fem::Hermite<2, 2, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag25scr64s )
{
    TestHermite<fem::Hermite<2, 5, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag210scr64s )
{
    TestHermite<fem::Hermite<2, 10, Scalar, real64_type, SimplexProduct> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // hermite_2d_simplexproduct_testsuite

BOOST_AUTO_TEST_SUITE( hermite_3d_simplexproduct_testsuite )

BOOST_AUTO_TEST_CASE( lag31scr64s )
{
    TestHermite<fem::Hermite<3, 1, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag32scr64s )
{
    TestHermite<fem::Hermite<3, 2, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag35scr64s )
{
    TestHermite<fem::Hermite<3, 5, Scalar, real64_type, SimplexProduct> > t;
    t();
}
BOOST_AUTO_TEST_CASE( lag310scr64s )
{
    TestHermite<fem::Hermite<3, 7, Scalar, real64_type, SimplexProduct> > t;
    t();
}

BOOST_AUTO_TEST_SUITE_END() // hermite_3d_simplexproduct_testsuite

#endif// 0
BOOST_AUTO_TEST_SUITE_END()
#else


#if defined( USE_TEST )
void test_scalar_continuous_simplex( test_suite* test )
{
#if USE_TEST
    // 1D simplex
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<1, 1, Scalar, real64_type, Simplex> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<1, 2, Scalar, real64_type, Simplex> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<1, 3, Scalar, real64_type, Simplex> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<1, 4, Scalar, real64_type, Simplex> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<1, 10, Scalar, real64_type, Simplex> >() ) ) );

    // 2D simplex
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<2, 1, Scalar, real64_type, Simplex> >())  ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<2, 2, Scalar, real64_type, Simplex> >())  ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<2,10, Scalar, real64_type, Simplex> >())  ) );

    // 3D simplex
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<3, 1, Scalar, real64_type, Simplex> >())  ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<3, 2, Scalar, real64_type, Simplex> >())  ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<3, 3, Scalar, real64_type, Simplex> >())  ) );
#else

    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<3, 3, Scalar, qd_real, Simplex> >())  ) );

#endif
}

void test_scalar_continuous_simplex_product( test_suite* test )
{
#if USE_TEST
    // 1D simplex
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<1, 1, Scalar, real64_type, SimplexProduct> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<1, 2, Scalar, real64_type, SimplexProduct> >() ) ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<1,10, Scalar, real64_type, SimplexProduct> >() ) ) );

    // 2D simplex
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<2, 1, Scalar, real64_type, SimplexProduct> >())  ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<2, 2, Scalar, real64_type, SimplexProduct> >())  ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<2,10, Scalar, real64_type, SimplexProduct> >())  ) );

    // 3D simplex
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<3, 1, Scalar, real64_type, SimplexProduct> >())  ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<3, 2, Scalar, real64_type, SimplexProduct> >())  ) );
    test->add( BOOST_TEST_CASE( ( TestHermite<fem::Hermite<3, 3, Scalar, real64_type, SimplexProduct> >())  ) );


#endif
}

test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
#if defined( HAVE_QD_QD_H )
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
#endif

    test_suite* test = BOOST_TEST_SUITE( "Finite element test suite" );


    test_scalar_continuous_simplex( test );
    //test_scalar_continuous_simplex_product( test );



    return test;
}
#else
int main( int argc, char* argv[] )
{
#if 1
    TestHermite<fem::Hermite<1, 3, Scalar, real64_type, Simplex> > t; t();
#if 0
    TestHermite<fem::Hermite<2, 2, Scalar, real64_type, Simplex> >();
    TestHermite<fem::Hermite<2, 10, Scalar, real64_type, Simplex> >();

    TestHermite<fem::Hermite<3, 1, Scalar, real64_type, Simplex> >();
    TestHermite<fem::Hermite<3, 2, Scalar, real64_type, Simplex> >();
    TestHermite<fem::Hermite<3, 3, Scalar, real64_type, Simplex> >();
#endif // 0
#endif

#if 0
    // 1D simplex
    ( ( ( TestHermite<fem::Hermite<1, 1, Scalar, real64_type, SimplexProduct> >() ) ) );
    ( ( ( TestHermite<fem::Hermite<1, 2, Scalar, real64_type, SimplexProduct> >() ) ) );
    ( ( ( TestHermite<fem::Hermite<1,10, Scalar, real64_type, SimplexProduct> >() ) ) );

    // 2D simplex
    ( ( ( TestHermite<fem::Hermite<2, 1, Scalar, real64_type, SimplexProduct> >())  ) );
    ( ( ( TestHermite<fem::Hermite<2, 2, Scalar, real64_type, SimplexProduct> >())  ) );
    ( ( ( TestHermite<fem::Hermite<2,10, Scalar, real64_type, SimplexProduct> >())  ) );

    // 3D simplex
    ( ( ( TestHermite<fem::Hermite<3, 1, Scalar, real64_type, SimplexProduct> >())  ) );
    ( ( ( TestHermite<fem::Hermite<3, 2, Scalar, real64_type, SimplexProduct> >())  ) );
    ( ( ( TestHermite<fem::Hermite<3, 3, Scalar, real64_type, SimplexProduct> >())  ) );
#endif
}

#endif
#endif // 0

