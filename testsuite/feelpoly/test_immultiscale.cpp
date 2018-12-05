/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Thomas Lantz 
       Date: 2015-04-27

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
   \file test_immultiscale.cpp
   \author Thomas Lantz
   \date 2005-07-29
 */
/*#define BOOST_TEST_MODULE integration methods test
// Boost.Test
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
*/
#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_integrateQuadra
#include <feel/feelcore/testsuite.hpp>
#endif

#include <feel/feelcore/feel.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/imsimplex.hpp>
#include <feel/feelpoly/multiscalequadrature.hpp>



using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "test_immultiscale" ,
                     "test_immultiscale" ,
                     "8.9e-3",
                     "test integrate Quadra",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Thomas Lantz", "student", "", "" );
    return about;
}

template<typename T>
T one( typename Feel::node<T>::type const& /*t*/ )
{
    return 1.0;
}

template<typename T>
T sint( typename Feel::node<T>::type const& t )
{
    return sin( t[0] );
}
template<typename T>
T cost( typename Feel::node<T>::type const& t )
{
    return cos( t[0] );
}

template<typename T>
T xp2( typename Feel::node<T>::type const& t )
{
    //   return pow(t[0]*t[1],T(P));
    return t[0]*t[0];
}
template<typename T>
T xp3( typename Feel::node<T>::type const& t )
{
    //   return pow(t[0]*t[1],T(P));
    return t[0]*t[0]*t[0];
}
template<int P, typename T>
T exp2( typename Feel::node<T>::type const& t )
{
    //   return pow(t[0]*t[1],T(P));
    return exp( t[0]+t[1] );
}
template<typename T>
T coscos( typename Feel::node<T>::type const& t )
{
    return cos( t[0] )*cos( t[1] );
}


template<int D, int N, typename T>
class TestImQK
{
public:
    typedef T value_type;
    typedef typename Feel::node<T>::type node_type;
    
    TestImQK( value_type res,
              boost::function<value_type( node_type const& )> const& func,
              value_type tol = Feel::type_traits<T>::epsilon() )
        :
        M_res( res ),
        M_tol( tol ),
        M_func( func )
    {}
    value_type operator()() const
    {
        using namespace Feel;

        //IMGeneral<D, T, Hypercube, MultiScaleQuadrature> im( N );
        IMGeneral<D, T, Hypercube> im( N );

        value_type res = math::abs( im.integrateAtPoints( M_func ) - M_res );

        if ( !( res < M_tol ) )
        {
            std::cout << "TestImQK:: int   = " << im.integrateAtPoints( M_func ) << "\n";
            std::cout << "TestImQK:: exact = " << M_res << "\n";
            std::cout << "TestImQK:: res   = " << res << "\n";
            std::cout << "TestImQK:: tol   = " << M_tol << "\n";
        }

        BOOST_CHECK( res < M_tol );

        std::cout << "" << std::endl;

        return res;
    }
    value_type M_res;
    value_type M_tol;
    boost::function<value_type( node_type const& )> M_func;
};


// automatically registered test cases could be organized in test suites
//BOOST_FIXTURE_TEST_SUITE( im1d_double_suite, F )

//#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( immultiscale_suite )
BOOST_AUTO_TEST_CASE( im1d_test1 )
{
    TestImQK<1,1, double> t1( 2.0, one<double> );
    t1();
}

BOOST_AUTO_TEST_CASE( im1d_test2 )
{
    TestImQK<1,1, double> t2( 2.0*sin( double( 1.0 ) )  , cost<double>, 0.1 );
    t2();
}

BOOST_AUTO_TEST_CASE( im1d_test3 )
{
    TestImQK<1,2, double> t3( 2.0*sin( double( 1.0 ) )  , cost<double>, 0.1);
    t3();
}

BOOST_AUTO_TEST_CASE( im1d_test4 )
{
    TestImQK<1,3, double> t4( 2.0*sin( double( 1.0 ) )  , cost<double>, 0.1 );
    t4();
}

BOOST_AUTO_TEST_CASE( im1d_test5 )
{
    TestImQK<1,50,double> t5( 2.0*sin( double( 1.0 ) )  , cost<double>, 0.1 );
    t5();
}

BOOST_AUTO_TEST_CASE( im2d_test0 )
{
    TestImQK<2,1, double> t1( 0 , sint<double>, 0.1 );
    t1();
}

BOOST_AUTO_TEST_CASE( im2d_test1 )
{
    TestImQK<2,1, double> t1( 4.*sin(1.) , cost<double> , 0.1 );
    t1();
}

BOOST_AUTO_TEST_CASE( im2d_test2 )
{
    TestImQK<2,1, double> t1( -2.+exp(2.)+(1/exp(2.)) , exp2<2,double> , 0.1);
    t1();
}

BOOST_AUTO_TEST_CASE( im2d_test3 )
{
    TestImQK<2,1, double> t1( 4./3. , xp2<double> , 0.1 );
    t1();
}

BOOST_AUTO_TEST_CASE( im2d_test4 )
{
    TestImQK<2,1, double> t1( 0 , xp3<double> , 0.1 );
    t1();
}

BOOST_AUTO_TEST_CASE( im2d_test5 )
{
    TestImQK<2,1, double> t1( 4.*sin(1.)*sin(1.) , coscos<double> , 0.1);
    t1();
}

BOOST_AUTO_TEST_CASE( im2d_test6 )
{
    TestImQK<2,50,double> t6( 4.0*sin( double( 1.0 ) )  , cost<double> , 0.1);
    t6();
}


BOOST_AUTO_TEST_SUITE_END()

