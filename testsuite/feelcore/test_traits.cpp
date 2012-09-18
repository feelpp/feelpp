/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-28

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2009 Université de Grenoble 1 (Joseph Fourier)

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
   \file test_traits.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-28
 */
// Boost.Test
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;


#include <feel/feelcore/constants.hpp>
#include <feel/feelcore/traits.hpp>


template<typename value_type>
void check( value_type x )
{
    BOOST_CHECK( Feel::math::abs( x ) < Feel::type_traits<value_type>::epsilon() );
}

BOOST_AUTO_TEST_CASE( test_functions )
{
    using namespace Feel;

# define FEELPP_CHECK_UNARY_FUNCS_OP_STD(_,TF) \
  FEELPP_CHECK_UNARY_FUNCS_OP_CODE_STD TF      \
  /**/
#
# define FEELPP_CHECK_UNARY_FUNCS_OP_CODE_STD(T,F)                                                                                     \
  check( Feel::math::FEELPP_FUNC_NAME( F )( FEELPP_TRAITS_TYPE( T )( 1.0 ) ) - std::FEELPP_FUNC_NONS( F )( FEELPP_TRAITS_TYPE( T )( 1.0 ) ) ); \
  /**/
#
    BOOST_PP_LIST_FOR_EACH_PRODUCT( FEELPP_CHECK_UNARY_FUNCS_OP_STD, 2, ( BOOST_PP_LIST_FILTER( FEELPP_PRED_FUNC, 0, FEELPP_TRAITS_TYPES ), FEELPP_STD_FUNCS ) );

# define FEELPP_CHECK_UNARY_FUNCS_OP_GLOBAL(_,TF) \
  FEELPP_CHECK_UNARY_FUNCS_OP_CODE_GLOBAL TF      \
  /**/
#
# define FEELPP_CHECK_UNARY_FUNCS_OP_CODE_GLOBAL(T,F)                                                                                  \
  check( Feel::math::FEELPP_FUNC_NAME( F )( FEELPP_TRAITS_TYPE( T )( 1.0 ) ) - ::FEELPP_FUNC_NONS( F )( FEELPP_TRAITS_TYPE( T )( 1.0 ) ) ); \
  /**/
#
    BOOST_PP_LIST_FOR_EACH_PRODUCT( FEELPP_CHECK_UNARY_FUNCS_OP_GLOBAL, 2, ( BOOST_PP_LIST_FILTER( FEELPP_PRED_FUNC, 1, FEELPP_TRAITS_TYPES ), FEELPP_GLOBAL_FUNCS ) );
}

BOOST_AUTO_TEST_CASE( test_promote )
{
    using namespace Feel;
    // integral types
#if !defined(__APPLE__)
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<int, uint>::type, int>::type::value ) );
#endif
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<int, uint64_type>::type, uint64_type>::type::value ) );

    // floating types
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<double, float>::type, double>::type::value ) );
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<float, double>::type, double>::type::value ) );
    //BOOST_CHECK( ( boost::is_same<strongest_numeric_type<double, long double>::type, long double>::type::value ) );
    //BOOST_CHECK( ( boost::is_same<strongest_numeric_type<long double, double>::type, long double>::type::value ) );
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<float,double, int>::type, double>::type::value ) );

    // mix float/integral type
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<double, int>::type, double>::type::value ) );
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<double, uint64_type>::type, double>::type::value ) );

    // multiples types: not working properly yet
    //BOOST_MPL_ASSERT(( boost::is_same< strongest_numeric_type<float,double,long double, std::complex<double> >::type, std::complex<long double> > ));

    BOOST_MPL_ASSERT( ( boost::is_same< strongest_numeric_type<double, std::complex<double> >::type, std::complex<double> > ) );
    BOOST_MPL_ASSERT( ( boost::is_same< strongest_numeric_type<std::complex<float>, double  >::type, std::complex<double> > ) );

    BOOST_MPL_ASSERT( ( boost::is_same< strongest_numeric_type<std::complex<float>, std::complex<double>  >::type, std::complex<double> > ) );



#if defined(FEELPP_HAS_QD_REAL)
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<dd_real, float>::type, dd_real>::type::value ) );
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<dd_real, double>::type, dd_real>::type::value ) );
    //BOOST_CHECK( ( boost::is_same<strongest_numeric_type<dd_real, long double>::type, dd_real>::type::value ) );
    //BOOST_CHECK( ( boost::is_same<strongest_numeric_type<qd_real, long double>::type, qd_real>::type::value ) );
    BOOST_CHECK( ( boost::is_same<strongest_numeric_type<qd_real, dd_real>::type, qd_real>::type::value ) );
    BOOST_CHECK( ( !boost::is_same<strongest_numeric_type<qd_real, dd_real>::type, dd_real>::type::value ) );

    // multiples types
    BOOST_MPL_ASSERT( ( boost::is_same< strongest_numeric_type<float,double,dd_real, qd_real>::type, qd_real> ) );
#endif
}


BOOST_AUTO_TEST_CASE( test_constants )
{
    using namespace Feel;
    BOOST_MESSAGE( "check pi value with double" );
    check( math::Constant<math::pi_tag, double >() - double( M_PI ) );
    BOOST_MESSAGE( "check pi value with float" );
    check( math::Constant<math::pi_tag, float >() - float( M_PI ) );

#if defined(FEELPP_HAS_QD_REAL)
    BOOST_MESSAGE( "check pi value with dd/qd_real" );
    check( math::Constant<math::pi_tag, dd_real >() - dd_real::_pi );
    check( math::Constant<math::pi_tag, qd_real >() - qd_real::_pi );
#endif

}
