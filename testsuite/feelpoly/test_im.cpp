/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-29

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007, 2008 Université Joseph Fourier Grenoble 1

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
   \file test_im.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-29
 */
#define BOOST_TEST_MODULE integration methods test
// Boost.Test
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <feel/feelcore/feel.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/imsimplex.hpp>



using namespace Feel;
struct F
{
    F() : i( 0 )
    {
        BOOST_TEST_MESSAGE( "setup fixture" );
    }
    ~F()
    {
        BOOST_TEST_MESSAGE( "teardown fixture" );
    }

    int i;
};

template<typename T>
T two( typename Feel::node<T>::type const& /*t*/ )
{
    return 2.0;
}

template<typename T>
T one( typename Feel::node<T>::type const& /*t*/ )
{
    return 1.0;
}

template<int P, typename T>
T xp( typename Feel::node<T>::type const& t )
{
    if ( P ==0 ) return 1.0;

    T res = t[0];

    for ( int i = 1; i < P; ++ i )
        res *= t[0];

    //std::cout << "res xp <" << P << ">: " << res << "\n";
    return res;
    //return pow( t, T(P) );
}
template<typename T>
struct P2N
{
    typedef T value_type;
    typedef typename Feel::node<T>::type node_type;

    P2N( int alpha1, int alpha2 )
        :
        M_alpha1( alpha1 ),
        M_alpha2( alpha2 )
    {}

    value_type operator()( node_type const& t )
    {
        value_type res = Feel::math::pow( t[0], M_alpha1 )* Feel::math::pow( t[1], M_alpha2 );
        return res;
    }
    value_type integral() const
    {
        value_type a1 = Feel::math::pow( -1.0, M_alpha2+1 )/( M_alpha2 +1 );
        value_type a2 = 1./( M_alpha1+M_alpha2+2 )-1./( M_alpha1+1 );
        value_type a3 = Feel::math::pow( -1.0, M_alpha1+M_alpha2+2 )/( M_alpha1 + M_alpha2+2 );
        value_type a4 = Feel::math::pow( -1.0, M_alpha1+1 )/( M_alpha1 + 1 );

        return a1*( a2-( a3-a4 ) );
    }
    int M_alpha1;
    int M_alpha2;
};

template<typename T>
struct P3N
{
    typedef T value_type;
    typedef typename Feel::node<T>::type node_type;

    P3N( int alpha1, int alpha2 = 1, int alpha3 = 1 )
        :
        M_alpha1( alpha1 ),
        M_alpha2( alpha2 ),
        M_alpha3( alpha3 )
    {}

    value_type operator()( node_type const& t )
    {
        value_type res = Feel::math::pow( t[0], M_alpha1 )* Feel::math::pow( t[1], M_alpha2 )* Feel::math::pow( t[2], M_alpha3 );
        return res;
    }
    value_type integral( int face )
    {
        value_type a = M_alpha1;
        value_type b = M_alpha2;
        value_type c = M_alpha3;

        switch ( face )
        {
        default:
        case FACE_0:
            return 0;

        case FACE_1:
            return 0;

        case FACE_2:
            return 0;

        case FACE_3:
            return -( Feel::math::pow( -1.,a )+1 )*Feel::math::pow( -1.,c )/( ( a+1 )*( a+3 ) );
        }



    }
    value_type integral() const
    {
        value_type a = M_alpha1;
        return ( 2*Feel::math::pow( a,3. )*Feel::math::pow( ( -1. ),a )+21*Feel::math::pow( a,2. )*Feel::math::pow( ( -1. ),a )+
                 67*a*Feel::math::pow( -1.,a )+63*Feel::math::pow( -1.,a )+3*Feel::math::pow( a,2. )+21*a+33 )/( 3*( a+1 )*( a+2 )*( a+3 )*( a+4 )*( a+5 ) );
    }
    int M_alpha1;
    int M_alpha2;
    int M_alpha3;
};

template<int P, typename T>
T xpdim( typename Feel::node<T>::type const& t )
{
    if ( P ==0 ) return 1.0;

    T res = t[0];

    for ( int i = 1; i < P; ++ i )
        res *= t[i];

    //std::cout << "res xp <" << P << ">: " << res << "\n";
    return res;
    //return pow( t, T(P) );
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
T xp1( typename Feel::node<T>::type const& t )
{
    //   return pow(t[0]*t[1],T(P));
    return t[0];
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
template<typename T>
T xp4( typename Feel::node<T>::type const& t )
{
    //   return pow(t[0]*t[1],T(P));
    return t[0]*t[0]*t[0]*t[0];
}
template<int P, typename T>
T exp2( typename Feel::node<T>::type const& t )
{
    //   return pow(t[0]*t[1],T(P));
    return exp( t[0]+t[1] );
}
template<typename T>
T sint2( typename Feel::node<T>::type const& t )
{
    return sin( t[0] )*sin( t[1] );
}
template<typename T>
T coscos( typename Feel::node<T>::type const& t )
{
    return cos( t[0] )*cos( t[1] );
}

template<Feel::uint16_type D,
         Feel::uint16_type N,
         typename T,
         template<class Convex, Feel::uint16_type O, typename T2> class QPS = Feel::Gauss>
class TestImPK
{
public:
    typedef T value_type;
    typedef typename Feel::node<value_type>::type node_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<D>,mpl::int_<1> >,
            mpl::identity<node_type>,
            mpl::identity<node_type> >::type::type arg_type;
    //IM<D, N, T, Simplex, QPS> im;
    template<int IMORDER> struct MyIm
            :
        public IM<D,N,T,Simplex>
    {};
    TestImPK( value_type res,
              boost::function<value_type( arg_type const& )> const& func,
              value_type tol = Feel::type_traits<T>::epsilon() )
        :
        M_res( res ),
        M_tol( tol ),
        M_face( -2 ),
        M_func( func )
    {

    }
    TestImPK( int __face,
              value_type res,
              boost::function<value_type( arg_type const& )> const& func,
              value_type tol = Feel::type_traits<T>::epsilon() )
        :
        M_res( res ),
        M_tol( tol ),
        M_face( __face ),
        M_func( func )
    {

    }
    value_type operator()() const
    {


        MyIm<N> im;
        value_type res = 0.0;
        value_type integral = 0.0;

        if ( M_face == -2 )
        {
            integral = im.integrateAtPoints( M_func );
        }

        if ( M_face != -2 )
        {
            integral = im.integrateAtPoints( Feel::IntegrationFaceEnum( M_face ), M_func );
        }

        res = math::abs( integral - M_res );

        if ( ! ( res < M_tol ) )
        {
            std::cout << "IM: " << im << "\n";
            std::cout << "TestIm<"<<D << "," << typeid( T ).name() << ">::  integral   = " << integral << "\n";
            std::cout << "TestIm<"<<D << "," << typeid( T ).name()<< ">:: exact value = " << M_res << "\n";
            std::cout << "TestIm<"<<D << "," << typeid( T ).name()<< ">::     error   = " << res << "\n";
            std::cout << "TestIm<"<<D << "," << typeid( T ).name()<< ">:: tolerance   = " << M_tol << "\n";
        }

        BOOST_CHECK( res < M_tol );
        return res;
    }
    value_type M_res;
    value_type M_tol;
    int M_face;
    boost::function<value_type( arg_type const& )> M_func;
};
#define TESTS_XP( T )                                                          \
 test->add( BOOST_TEST_CASE( ( TestImPK<1,1, T>( 0.0         , xp<1,T> ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<1,10,T>( T(2.0)/11.0 , xp<10,T> ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<1,2, T>( T(2.0)/3.0  , xp<2,T> ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<1,3, T>( 0.0         , xp<3,T> ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<1,20,T>( T(2.0)/21.0 , xp<20,T>) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<1,21,T>( 0.0         , xp<21,T>) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<1,15,T>( 0.0         , xp<21,T>) )  ) ); \
 /**/
#define TESTS_SIN( T )                                                    \
 test->add( BOOST_TEST_CASE( ( TestImPK<1,1, T>( 0.0         , sint<T> ) )  ) );\
 test->add( BOOST_TEST_CASE( ( TestImPK<1,2, T>( 0.0         , sint<T> ) )  ) );\
 test->add( BOOST_TEST_CASE( ( TestImPK<1,3, T>( 0.0         , sint<T> ) )  ) );\
 test->add( BOOST_TEST_CASE( ( TestImPK<1,20,T>( 0.0         , sint<T> ) ) )  );\
 /**/

#define TESTS_XP2( T )                                                          \
 test->add( BOOST_TEST_CASE( ( TestImPK<2,5, T>( -T(2.0)/3.0 , xp2<1,T> ) )  ) ); \
    test->add( BOOST_TEST_CASE( ( TestImPK<2,17, T>( 1.0+T(1.0)/(exp(T(1.0))*exp(T(1.0))) , exp2<1,T>, 1e-15 ) )  ) ); \
 /**/
#define TESTS_SIN2( T )                                                    \
 test->add( BOOST_TEST_CASE( ( TestImPK<2,1, T>( 0.0         , sint2<T>, 1e-1 ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<2,2, T>( 0.0         , sint2<T>, 1.5e-2 ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<2,3, T>( 0.0         , sint2<T>, 5e-4 ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<2,4, T>( 0.0         , sint2<T>, 5e-4 ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<2,6, T>( 0.0         , sint2<T>, 5e-4 ) )  ) ); \
 test->add( BOOST_TEST_CASE( ( TestImPK<2,20,T>( 0.0         , sint2<T>, 1.5e-16 ) ) )  ); \
 /**/
#define TESTS_COSCOS( T )                                                    \
 test->add( BOOST_TEST_CASE( ( TestImPK<2,40,T>( 2.0*sin(T(1))*sin(T(1)) , coscos<T> ) ) )  );\
 /**/


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

        IM<D, N, T, Hypercube> im;

        value_type res = math::abs( im.integrateAtPoints( M_func ) - M_res );

        if ( !( res < M_tol ) )
        {
            std::cout << "TestImQK:: int   = " << im.integrateAtPoints( M_func ) << "\n";
            std::cout << "TestImQK:: exact = " << M_res << "\n";
            std::cout << "TestImQK:: res   = " << res << "\n";
            std::cout << "TestImQK:: tol   = " << M_tol << "\n";
        }

        BOOST_CHECK( res < M_tol );
        return res;
    }
    value_type M_res;
    value_type M_tol;
    boost::function<value_type( node_type const& )> M_func;
};
#define TESTS_QK_COSCOS( T )                                                    \
 test->add( BOOST_TEST_CASE( ( TestImQK2D<40, T>( 4.0*sin(T(1.0))*sin(T(1.0)) , coscos<T>, 1.0E-14 ) )  ) );\
 /**/


#if 1

// automatically registered test cases could be organized in test suites
BOOST_FIXTURE_TEST_SUITE( im1d_double_suite, F )
BOOST_AUTO_TEST_CASE( im1d_test1 )
{
    TestImPK<1,1, double> t1( 2.0, one<double> );
    t1();
}
BOOST_AUTO_TEST_CASE( im1d_test2 )
{
    TestImPK<1,1, double> t2( 2.0*sin( double( 1.0 ) )  , cost<double>, 1.1E-2 );
    t2();
}
BOOST_AUTO_TEST_CASE( im1d_test3 )
{
    TestImPK<1,2, double> t3( 2.0*sin( double( 1.0 ) )  , cost<double>, 1.1E-2 );
    t3();
}
BOOST_AUTO_TEST_CASE( im1d_test4 )
{
    TestImPK<1,3, double> t4( 2.0*sin( double( 1.0 ) )  , cost<double>, 1.0E-4 );
    t4();
}
BOOST_AUTO_TEST_CASE( im1d_test5 )
{
    TestImPK<1,50,double> t5( 2.0*sin( double( 1.0 ) )  , cost<double> );
    t5();
}
BOOST_AUTO_TEST_CASE( im1d_test6 )
{
    TestImPK<1,1, double> t6( Feel::FACE_0, -1.0 , xp<1,double> );
    t6();
}
BOOST_AUTO_TEST_CASE( im1d_test7 )
{
    TestImPK<1,1, double> t7( Feel::FACE_1, 1.0 , xp<1,double> );
    t7();
}
BOOST_AUTO_TEST_CASE( im1d_face_test1 )
{
    TestImPK<1,1, double> t0( Feel::FACE_0, -1.0 , xp<1,double> );
    t0();
    TestImPK<1,1, double> t1( Feel::FACE_1, 1.0 , xp<1,double> );
    t1();
}
BOOST_AUTO_TEST_CASE( im1d_face_test2 )
{
    TestImPK<1,2, double> t0( Feel::FACE_0, -1.0 , xp<1,double> );
    t0();
    TestImPK<1,2, double> t1( Feel::FACE_1, 1.0 , xp<1,double> );
    t1();
}
BOOST_AUTO_TEST_CASE( im1d_face_test3 )
{
    TestImPK<1,3, double> t0( Feel::FACE_0, -1.0 , xp<1,double> );
    t0();
    TestImPK<1,3, double> t1( Feel::FACE_1, 1.0 , xp<1,double> );
    t1();
}

BOOST_AUTO_TEST_CASE( im2d_test0 )
{

    TestImPK<2,1, double> t1( P2N<double>( 0, 0 ).integral(), P2N<double>( 0, 0 ) );
    t1();
}
BOOST_AUTO_TEST_CASE( im2d_test1 )
{
    const int N = 1;

    for ( int alpha1 = 0; alpha1 < N; ++ alpha1 )
    {
        int alpha2 = N-alpha1;
        TestImPK<2,N, double> t10( P2N<double>( alpha1, alpha2 ).integral(), P2N<double>( alpha1, alpha2 ) );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im2d_test2 )
{
    const int N = 2;

    for ( int alpha1 = 0; alpha1 < N; ++ alpha1 )
    {
        int alpha2 = N-alpha1;
        TestImPK<2,N, double> t10( P2N<double>( alpha1, alpha2 ).integral(), P2N<double>( alpha1, alpha2 ) );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im2d_test3 )
{
    const int N = 3;

    for ( int alpha1 = 0; alpha1 < N; ++ alpha1 )
    {
        int alpha2 = N-alpha1;
        TestImPK<2,N, double> t10( P2N<double>( alpha1, alpha2 ).integral(), P2N<double>( alpha1, alpha2 ) );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im2d_test4 )
{
    const int N = 4;

    for ( int alpha1 = 0; alpha1 < N; ++ alpha1 )
    {
        int alpha2 = N-alpha1;
        TestImPK<2,N, double> t10( P2N<double>( alpha1, alpha2 ).integral(), P2N<double>( alpha1, alpha2 ) );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im2d_test5 )
{
    const int N = 5;

    for ( int alpha1 = 0; alpha1 < N; ++ alpha1 )
    {
        int alpha2 = N-alpha1;
        TestImPK<2,N, double> t10( P2N<double>( alpha1, alpha2 ).integral(), P2N<double>( alpha1, alpha2 ) );
        t10();
    }

}
BOOST_AUTO_TEST_CASE( im2d_test6 )
{
    const int N = 6;

    for ( int alpha1 = 0; alpha1 < N; ++ alpha1 )
    {
        int alpha2 = N-alpha1;
        TestImPK<2,N, double> t10( P2N<double>( alpha1, alpha2 ).integral(), P2N<double>( alpha1, alpha2 ) );
        t10();
    }


}
BOOST_AUTO_TEST_CASE( im2d_face_test1 )
{
    TestImPK<2,1, double> t0( Feel::FACE_1, 2.0 , one<double> );
    t0();
    TestImPK<2,1, double> t1( Feel::FACE_0, Feel::math::sqrt( double( 8.0 ) ) , one<double> );
    t1();
    TestImPK<2,1, double> t2( Feel::FACE_2, 2.0 , one<double> );
    t2();
}
BOOST_AUTO_TEST_CASE( im2d_face_test2 )
{
    TestImPK<2,2, double> t0( Feel::FACE_1, 2.0 , one<double> );
    t0();
    TestImPK<2,2, double> t1( Feel::FACE_0, Feel::math::sqrt( double( 8.0 ) ) , one<double> );
    t1();
    TestImPK<2,2, double> t2( Feel::FACE_2, 2.0 , one<double> );
    t2();
}
BOOST_AUTO_TEST_CASE( im2d_face_test3 )
{
    TestImPK<2,3, double> t0( Feel::FACE_1, 2.0 , one<double> );
    t0();
    TestImPK<2,3, double> t1( Feel::FACE_0, Feel::math::sqrt( double( 8.0 ) ) , one<double> );
    t1();
    TestImPK<2,3, double> t2( Feel::FACE_2, 2.0 , one<double> );
    t2();
}
BOOST_AUTO_TEST_CASE( im2d_face_test4 )
{
    TestImPK<2,4, double> t0( Feel::FACE_1, 2.0 , one<double> );
    t0();
    TestImPK<2,4, double> t1( Feel::FACE_0, Feel::math::sqrt( double( 8.0 ) ) , one<double> );
    t1();
    TestImPK<2,4, double> t2( Feel::FACE_2, 2.0 , one<double> );
    t2();
}
BOOST_AUTO_TEST_CASE( im2d_face_test5 )
{
    TestImPK<2,5, double> t0( Feel::FACE_1, 2.0 , one<double> );
    t0();
    TestImPK<2,5, double> t1( Feel::FACE_0, Feel::math::sqrt( double( 8.0 ) ) , one<double> );
    t1();
    TestImPK<2,5, double> t2( Feel::FACE_2, 2.0 , one<double> );
    t2();
}
BOOST_AUTO_TEST_CASE( im2d_face_test6 )
{
    TestImPK<2,6, double> t0( Feel::FACE_1, 2.0 , one<double> );
    t0();
    TestImPK<2,6, double> t1( Feel::FACE_0, Feel::math::sqrt( double( 8.0 ) ) , one<double> );
    t1();
    TestImPK<2,6, double> t2( Feel::FACE_2, 2.0 , one<double> );
    t2();
}

BOOST_AUTO_TEST_CASE( im3d_test1 )
{
    const int N = 2;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test2 )
{
    const int N = 3;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test3 )
{
    const int N = 4;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test5 )
{
    const int N = 5;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test6 )
{
    const int N = 6;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test7 )
{
    const int N = 7;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test11 )
{
    const int N = 11;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test13 )
{
    const int N = 13;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test15 )
{
    const int N = 15;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test17 )
{
    const int N = 17;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test19 )
{
    const int N = 19;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_test20 )
{
    const int N = 20;

    for ( int alpha1 = 0; alpha1 <= N-2; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = 1;
        TestImPK<3,N, double> t10( P3N<double>( alpha1, alpha2, alpha3 ).integral(), P3N<double>( alpha1, alpha2, alpha3 ), 5e-13 );
        t10();
    }
}

BOOST_AUTO_TEST_CASE( im3d_face_test1 )
{
    TestImPK<3,1, double> t1( Feel::FACE_2, 2.0 , one<double> );
    t1();
    TestImPK<3,1, double> t2( Feel::FACE_1, 2.0 , one<double> );
    t2();
    TestImPK<3,1, double> t3( Feel::FACE_0, 2.0*Feel::math::sqrt( double( 3.0 ) ), one<double> );
    t3();
    TestImPK<3,1, double> t4( Feel::FACE_3, 2.0 , one<double> );
    t4();

    const int N = 1;

    for ( int alpha1 = 0; alpha1 < N-1; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = N-1-alpha1;
        TestImPK<3,N, double> t10( Feel::FACE_3, P3N<double>( alpha1, alpha2, alpha3 ).integral( FACE_3 ), P3N<double>( alpha1, alpha2, alpha3 ) );
        t10();
    }

}

BOOST_AUTO_TEST_CASE( im3d_face_test2 )
{
    const int N = 2;
    TestImPK<3,N, double> t1( Feel::FACE_2, 2.0 , one<double> );
    t1();
    TestImPK<3,N, double> t2( Feel::FACE_1, 2.0 , one<double> );
    t2();
    TestImPK<3,N, double> t3( Feel::FACE_0, 2.0*Feel::math::sqrt( double( 3.0 ) ), one<double> );
    t3();
    TestImPK<3,N, double> t4( Feel::FACE_3, 2.0 , one<double> );
    t4();

    for ( int alpha1 = 0; alpha1 < N-1; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = N-1-alpha1;
        TestImPK<3,N, double> t10( Feel::FACE_3, P3N<double>( alpha1, alpha2, alpha3 ).integral( FACE_3 ), P3N<double>( alpha1, alpha2, alpha3 ) );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_face_test3 )
{
    const int N = 3;
    TestImPK<3,N, double> t1( Feel::FACE_2, 2.0 , one<double> );
    t1();
    TestImPK<3,N, double> t2( Feel::FACE_1, 2.0 , one<double> );
    t2();
    TestImPK<3,N, double> t3( Feel::FACE_0, 2.0*Feel::math::sqrt( double( 3.0 ) ), one<double> );
    t3();
    TestImPK<3,N, double> t4( Feel::FACE_3, 2.0 , one<double> );
    t4();

    for ( int alpha1 = 0; alpha1 < N-1; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = N-1-alpha1;
        TestImPK<3,N, double> t10( Feel::FACE_3, P3N<double>( alpha1, alpha2, alpha3 ).integral( FACE_3 ), P3N<double>( alpha1, alpha2, alpha3 ) );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_face_test4 )
{
    const int N = 4;
    TestImPK<3,N, double> t1( Feel::FACE_2, 2.0 , one<double> );
    t1();
    TestImPK<3,N, double> t2( Feel::FACE_1, 2.0 , one<double> );
    t2();
    TestImPK<3,N, double> t3( Feel::FACE_0, 2.0*Feel::math::sqrt( double( 3.0 ) ), one<double> );
    t3();
    TestImPK<3,N, double> t4( Feel::FACE_3, 2.0 , one<double> );
    t4();

    for ( int alpha1 = 0; alpha1 < N-1; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = N-1-alpha1;
        TestImPK<3,N, double> t10( Feel::FACE_3, P3N<double>( alpha1, alpha2, alpha3 ).integral( FACE_3 ), P3N<double>( alpha1, alpha2, alpha3 ) );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_face_test5 )
{
    const int N = 5;
    TestImPK<3,N, double> t1( Feel::FACE_2, 2.0 , one<double> );
    t1();
    TestImPK<3,N, double> t2( Feel::FACE_1, 2.0 , one<double> );
    t2();
    TestImPK<3,N, double> t3( Feel::FACE_0, 2.0*Feel::math::sqrt( double( 3.0 ) ), one<double> );
    t3();
    TestImPK<3,N, double> t4( Feel::FACE_3, 2.0 , one<double> );
    t4();

    for ( int alpha1 = 0; alpha1 < N-1; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = N-1-alpha1;
        TestImPK<3,N, double> t10( Feel::FACE_3, P3N<double>( alpha1, alpha2, alpha3 ).integral( FACE_3 ), P3N<double>( alpha1, alpha2, alpha3 ) );
        t10();
    }
}
BOOST_AUTO_TEST_CASE( im3d_face_test6 )
{
    const int N = 6;
    TestImPK<3,N, double> t1( Feel::FACE_2, 2.0 , one<double> );
    t1();
    TestImPK<3,N, double> t2( Feel::FACE_1, 2.0 , one<double> );
    t2();
    TestImPK<3,N, double> t3( Feel::FACE_0, 2.0*Feel::math::sqrt( double( 3.0 ) ), one<double> );
    t3();
    TestImPK<3,N, double> t4( Feel::FACE_3, 2.0 , one<double> );
    t4();

    for ( int alpha1 = 0; alpha1 < N-1; ++ alpha1 )
    {
        int alpha2 = 1;
        int alpha3 = N-1-alpha1;
        TestImPK<3,N, double> t10( Feel::FACE_3, P3N<double>( alpha1, alpha2, alpha3 ).integral( FACE_3 ), P3N<double>( alpha1, alpha2, alpha3 ) );
        t10();
    }
}

BOOST_AUTO_TEST_SUITE_END()

#else
template<typename T>
void add_tests( test_suite* test )
{
    using namespace Feel;
    typedef T double;

    // 1D quadrature on segment [-1;1]
    test->add( BOOST_TEST_CASE( ( TestImPK<1,1, value_type>( 2.0, one<value_type> ) )  ) );
    TESTS_XP( value_type );
    TESTS_SIN( value_type );
    test->add( BOOST_TEST_CASE( ( TestImPK<1,1, value_type>( 2.0*sin( value_type( 1.0 ) )  , cost<value_type>, 1.1E-2 ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<1,2, value_type>( 2.0*sin( value_type( 1.0 ) )  , cost<value_type>, 1.1E-2 ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<1,3, value_type>( 2.0*sin( value_type( 1.0 ) )  , cost<value_type>, 1.0E-4 ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<1,50,value_type>( 2.0*sin( value_type( 1.0 ) )  , cost<value_type> ) ) )  );
    test->add( BOOST_TEST_CASE( ( TestImPK<1,1, value_type>( Feel::FACE_0, -1.0 , xp<1,value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<1,1, value_type>( Feel::FACE_1, 1.0 , xp<1,value_type> ) )  ) );

    // 2D quadrature on triangle [-1;1]x[-1;-x]
    test->add( BOOST_TEST_CASE( ( TestImPK<2,1, value_type>( 2.0 , one<value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<2,1, value_type>( Feel::FACE_1, 2.0 , one<value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<2,1, value_type>( Feel::FACE_0, Feel::math::sqrt( value_type( 8.0 ) ) , one<value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<2,1, value_type>( Feel::FACE_2, 2.0 , one<value_type> ) )  ) );

    test->add( BOOST_TEST_CASE( ( TestImPK<2,5, value_type>( -2.0/value_type( 3.0 ) , xp2<1,value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<2,50, value_type>( value_type( 1.0 )+1.0/( exp( value_type( 1.0 ) )*exp( value_type( 1.0 ) ) ) , exp2<1,value_type> ) )  ) );

    // 2D Quadratures on faces of the reference triangle
    test->add( BOOST_TEST_CASE( ( TestImPK<2,1, value_type>( Feel::ALL_FACES, value_type( 4.0 )+sqrt( value_type( 8.0 ) ), one<value_type> ) )  ) );

    // 3D quadrature [-1;1]x-[-1;-x]x[-1;-x]
    test->add( BOOST_TEST_CASE( ( TestImPK<3,1, value_type>( 8.0/value_type( 6.0 ) , one<value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<3,1, value_type>( Feel::FACE_2, 2.0 , one<value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<3,1, value_type>( Feel::FACE_1, 2.0 , one<value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<3,1, value_type>( Feel::FACE_0, 2.0*Feel::math::sqrt( value_type( 3.0 ) ), one<value_type> ) )  ) );
    test->add( BOOST_TEST_CASE( ( TestImPK<3,1, value_type>( Feel::FACE_3, 2.0 , one<value_type> ) )  ) );

    TESTS_XP2( value_type );
    TESTS_SIN2( value_type );
    TESTS_COSCOS( value_type );

    // 2D quadrature
    test->add( BOOST_TEST_CASE( ( TestImQK<2,50, value_type>( 4.0*sin( value_type( 1.0 ) )*sin( value_type( 1.0 ) ) , coscos<value_type>, 1.8e-15 ) )  ) );
    // 3D quadrature
    test->add( BOOST_TEST_CASE( ( TestImQK<3,1, value_type>( 0.0, xpdim<3,value_type>, 1.8e-15 ) )  ) );
    // 4D quadrature
    test->add( BOOST_TEST_CASE( ( TestImQK<4,1, value_type>( 0.0, xpdim<4,value_type>, 1.8e-15 ) )  ) );
}

#if 1
test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{

    test_suite* test = BOOST_TEST_SUITE( "Integration methods test suite" );

    add_tests<double>( test );

#if defined( FEELPP_HAS_QD_QD_H )

    unsigned int old_cw;
    fpu_fix_start( &old_cw );

    add_tests<dd_real>( test );
    add_tests<qd_real>( test );

#endif /* FEELPP_HAS_QD_QD_H  */

#if defined( FEELPP_HAS_MPFR )
    std::cout << "[feelpoly::testsuite::test_im] test with arbitrary precision\n";
    mp::mp_init( 200 );
    add_tests<mp_type>( test );
#endif
    return test;
}
#else
int main()
{
    //unsigned int old_cw;
    //fpu_fix_start(&old_cw);

    TestImPK<1,1, double> im0( Feel::FACE_0, 1.0 , one<double> );
    im0();

    TestImPK<1,1, double> im1( Feel::FACE_1, 1.0 , one<double> );
    im1();
    TestImPK<1,1, double> im11( Feel::FACE_0,-1.0 , xp<1,double> );
    im11();
    TestImPK<1,1, double> im2( Feel::FACE_0, 1.0 , xp<2,double> );
    im2();
    TestImPK<1,1, double> im3( Feel::FACE_1, 1.0 , xp<2,double> );
    im3();

#if 0
#if 0
    //test_suite* test = BOOST_TEST_SUITE( "Integration methods test suite" );
    //add_tests<double>( test );
    TestImPK<2,1, double> im( 2.0 , one<double> );
    im();
    TestImPK<2,1, dd_real> imdd( 2.0 , one<dd_real> );
    imdd();
#endif

    TestImQK<2,1, double> im1( 0.0 , xp<2,dd_real> );
    im1();
    TestImQK<3,1, double> im2( 0.0 , xp<3,dd_real> );
    im2();

    Feel::Tesseract tess;
    std::cout << "Tess.Dim = " << tess.nDim << "\n";
    TestImQK<4,1, double> im3( 0.0 , xp<4,dd_real> );
    im3();
#endif
}
#endif

#endif // 1
