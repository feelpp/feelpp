/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-07-28

  Copyright (C) 2005,2006 EPFL

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
   \file test_jacobi.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-07-28
 */
#define CHAR_BIT 8
#include <boost/timer.hpp>
// Boost.Test
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;


#include <life/lifecore/traits.hpp>
#include <life/lifepoly/jacobi.hpp>
#include <life/lifepoly/expansions.hpp>
#include <life/lifepoly/dubiner.hpp>
#include <life/lifepoly/polynomial.hpp>
#include <life/lifepoly/geomap.hpp>
#include <life/lifepoly/lagrange.hpp>

uint fact( uint  n )
{
    if ( n == 0 )
        return 1;
    return n*fact( n-1 );
}
template<typename T>
T f( T const& t )
{
  return sin(t);
}
template<Life::uint16_type N, typename T>
class TestJacobi
{
public:
    typedef T value_type;
    value_type _M_factor_eps;
    TestJacobi( value_type factor_eps = value_type( 1 ) )
        :
        _M_factor_eps( factor_eps )
    {}
    void operator()() const
    {
        using namespace Life;

        boost::timer __timer;

        for(int n = 0; n < 10000; ++n )
        {
            Jacobi<N, T> p( 0.0, 0.0 );
            value_type v = value_type( double(fact( N+0 ))/value_type( double(fact( N )*fact( 0 )) ) );
            BOOST_CHECK( Life::math::abs( p( 1.0 ) - v ) < _M_factor_eps*Life::type_traits<T>::epsilon() );
            //BOOST_CHECK( Life::math::abs( Jacobi<N, T>( 1.0, 0.0, 1.0 ) -
            //( fact( n-0 )/( fact( n )*fact(0 ) ) ) ) < Life::type_traits<T>::epsilon() );

            if ( N > 0 )
                {
                    // check derivation with derivation relation in KS book appendix A p 351
                    Jacobi<N-1, T> dp( 1.0, 1.0 );
                    value_type dp2 = 0.5 * ( 0.0 + 0.0 + value_type(N) + 1.0 ) * dp( 0.6 );
                    BOOST_CHECK( Life::math::abs( p.derivate( 0.6 ) - dp2 ) <  _M_factor_eps*Life::type_traits<T>::epsilon() );
                }
            BOOST_CHECK( Life::math::abs( Life::details::integrate<N,T>( f<T> ) - 0.0 ) < _M_factor_eps*Life::type_traits<T>::epsilon() );

        }
        std::cout << "meta::jacobi test done in " << __timer.elapsed()/1000 << "\n";
        __timer.restart();
        for(int n = 0; n < 10000; ++n )
        {
            dyna::Jacobi<T> p( N, 0.0, 0.0 );
            value_type v = value_type( double(fact( N+0 ))/value_type( double(fact( N )*fact( 0 )) ) );
            BOOST_CHECK( Life::math::abs( p( 1.0 ) - v ) < _M_factor_eps*Life::type_traits<T>::epsilon() );
            //BOOST_CHECK( Life::math::abs( Jacobi<N, T>( 1.0, 0.0, 1.0 ) -
            //( fact( n-0 )/( fact( n )*fact(0 ) ) ) ) < Life::type_traits<T>::epsilon() );

            if ( N > 0 )
                {
                    // check derivation with derivation relation in KS book appendix A p 351
                    dyna::Jacobi<T> dp( uint16_type(N-1), 1.0, 1.0 );
                    value_type dp2 = 0.5 * ( 0.0 + 0.0 + value_type(N) + 1.0 ) * dp( 0.6 );
                    BOOST_CHECK( Life::math::abs( p.derivate( 0.6 ) - dp2 ) <  _M_factor_eps*Life::type_traits<T>::epsilon() );
                }
            BOOST_CHECK( Life::math::abs( Life::details::integrate<N,T>( f<T> ) - 0.0 ) < _M_factor_eps*Life::type_traits<T>::epsilon() );
        }
        std::cout << "dyna::jacobi test done in " << __timer.elapsed()/1000 << "\n";

        {
            Jacobi<N, T> p( 0.0, 0.0 );
            ublas::vector<T> pts( 2 );
            pts[0] = 1.0;
            pts[1] = 1.0;
            ublas::matrix<T> P ( JacobiBatchEvaluation<N, T>( 0, 0, pts ) );
            BOOST_CHECK( P.size1() == N+1 );
            BOOST_CHECK( P.size2() == 2 );
            value_type v = value_type( double(fact( N+0 ))/value_type( double(fact( N )*fact( 0 )) ) );
            BOOST_CHECK( Life::math::abs( p( 1.0 ) - v ) < _M_factor_eps*Life::type_traits<T>::epsilon() );
            BOOST_CHECK( Life::math::abs( P( N, 0 ) - v ) < _M_factor_eps*Life::type_traits<T>::epsilon() );
            if ( pts.size() > 0 )
                BOOST_CHECK( Life::math::abs( P( N, 1 ) - v ) < _M_factor_eps*Life::type_traits<T>::epsilon() );

            if ( N > 0 )
                {
                    // check derivation with derivation relation in KS book appendix A p 351
                    Jacobi<N-1, T> dp( 1.0, 1.0 );
                    value_type dp2 = 0.5 * ( 0.0 + 0.0 + value_type(N) + 1.0 ) * dp( 0.6 );
                    BOOST_CHECK( Life::math::abs( p.derivate( 0.6 ) - dp2 ) <  _M_factor_eps*Life::type_traits<T>::epsilon() );
                    ublas::vector<T> pts( 2 );
                    pts[0] = 0.6;
                    pts[1] = 0.6;
                    ublas::matrix<T> dP ( JacobiBatchDerivation<N, T>( 0.0, 0.0, pts ) );
                    if ( Life::math::abs( dP( N, 0 ) - dp2 ) > Life::type_traits<T>::epsilon() )
                        {
                            std::cout << "dP( N, 0 ) = " <<  dP( N, 0 ) << "\n"
                                      << "dp2 = " << dp2 << "\n";
                        }
                    BOOST_CHECK( Life::math::abs( dP( N, 0 ) - dp2 ) < _M_factor_eps*Life::type_traits<T>::epsilon() );
                }
        }
    }
};

template<Life::uint16_type N, typename T>
class TestDubiner
{
  public:
    typedef T value_type;
    TestDubiner()
    {}
    void operator()() const
    {
      using namespace Life;
      Dubiner<2,N,Normalized<true>,T> dubexp;
      ublas::matrix<T,ublas::column_major> pts( 2,2  );
      ublas::column( pts, 0 )  = ublas::scalar_vector<value_type>( pts.size2(), 1.0 );
      ublas::column( pts, 1 )  = ublas::scalar_vector<value_type>( pts.size2(), 1.0 );
      std::cout << "[batch] dubexp at pts =  " << dubexp.evaluate( pts ) << "\n";
      std::cout << "[batch] dubexp derivation at pts =  " << dubexp.derivate( pts ) << "\n";
      Dubiner<2,N> dub;
    }
};

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{
    test_suite* test = BOOST_TEST_SUITE( "Master test suite" );

#if 1
    //
    // TestJacobi
    //
    test->add( BOOST_TEST_CASE( ( TestJacobi<0, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( TestJacobi<1, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( TestJacobi<2, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( TestJacobi<10, double>() )  ) );
    //test->add( BOOST_TEST_CASE( ( TestJacobi<15, Life::real96_type>(50.0) )  ) );

#if defined( LIFE_HAVE_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start(&old_cw);

    test->add( BOOST_TEST_CASE( ( TestJacobi<0, qd_real>() )  ) );
    test->add( BOOST_TEST_CASE( ( TestJacobi<1, qd_real>() )  ) );
    test->add( BOOST_TEST_CASE( ( TestJacobi<2, qd_real>() )  ) );
    test->add( BOOST_TEST_CASE( ( TestJacobi<4, qd_real>() )  ) );
    test->add( BOOST_TEST_CASE( ( TestJacobi<10, qd_real>() )  ) );
    //test->add( BOOST_TEST_CASE( ( TestJacobi<15, qd_real>() )  ) );
#endif /* LIFE_HAVE_QD_REAL */

    //
    // TestDubiner
    //
    std::cout << "o- TestDubiner<0, double>()\n";
    test->add( BOOST_TEST_CASE( ( TestDubiner<0, double>() )  ) );
    std::cout << "o- TestDubiner<1, double>()\n";
    test->add( BOOST_TEST_CASE( ( TestDubiner<1, double>() )  ) );
#else
    //test->add( BOOST_TEST_CASE( ( TestJacobi<1, double>() )  ) );
    //test->add( BOOST_TEST_CASE( ( TestJacobi<10, double>() )  ) );

#if defined( LIFE_HAVE_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start(&old_cw);

    test->add( BOOST_TEST_CASE( ( TestJacobi<0, qd_real>() )  ) );
    //test->add( BOOST_TEST_CASE( ( TestJacobi<1, qd_real>() )  ) );
    //test->add( BOOST_TEST_CASE( ( TestJacobi<2, qd_real>() )  ) );
    //test->add( BOOST_TEST_CASE( ( TestJacobi<4, qd_real>() )  ) );
    test->add( BOOST_TEST_CASE( ( TestJacobi<10, qd_real>() )  ) );
    //test->add( BOOST_TEST_CASE( ( TestJacobi<15, qd_real>() )  ) );
#endif /* LIFE_HAVE_QD_REAL */

#endif
    return test;
}


