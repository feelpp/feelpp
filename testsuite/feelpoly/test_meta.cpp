/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-01-04

  Copyright (C) 2026 Feel++ Consortium

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
   \file test_meta.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-04
   \brief Tests for mp11-based metaprogramming utilities in feelpoly
 */

#define BOOST_TEST_MODULE test_meta
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/meta.hpp>

#include <array>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( meta_suite )

//
// Test type_list (mp::mp_list wrapper)
//
BOOST_AUTO_TEST_CASE( test_type_list )
{
    using list1 = type_list<int, double, float>;
    using list2 = type_list<>;
    using list3 = type_list<std::string>;

    // Check sizes using mp11
    static_assert( mp::mp_size<list1>::value == 3, "list1 should have 3 elements" );
    static_assert( mp::mp_size<list2>::value == 0, "list2 should be empty" );
    static_assert( mp::mp_size<list3>::value == 1, "list3 should have 1 element" );

    // Check that type_list is indeed mp_list
    static_assert( std::is_same_v<list1, mp::mp_list<int, double, float>>,
                   "type_list should be mp_list" );

    BOOST_CHECK( true );
}

//
// Test type_at (element access)
//
BOOST_AUTO_TEST_CASE( test_type_at )
{
    using list = type_list<int, double, float, char>;

    static_assert( std::is_same_v<type_at<list, 0>, int>, "Element 0 should be int" );
    static_assert( std::is_same_v<type_at<list, 1>, double>, "Element 1 should be double" );
    static_assert( std::is_same_v<type_at<list, 2>, float>, "Element 2 should be float" );
    static_assert( std::is_same_v<type_at<list, 3>, char>, "Element 3 should be char" );

    BOOST_CHECK( true );
}

//
// Test type_find (find type index)
//
BOOST_AUTO_TEST_CASE( test_type_find )
{
    using list = type_list<int, double, float, char>;

    static_assert( type_find<list, int>::value == 0, "int should be at index 0" );
    static_assert( type_find<list, double>::value == 1, "double should be at index 1" );
    static_assert( type_find<list, float>::value == 2, "float should be at index 2" );
    static_assert( type_find<list, char>::value == 3, "char should be at index 3" );

    // Not found returns size (mp11 behavior)
    static_assert( type_find<list, std::string>::value == 4,
                   "Not found should return list size" );

    BOOST_CHECK( true );
}

//
// Test type_transform (apply metafunction to list)
//
BOOST_AUTO_TEST_CASE( test_type_transform )
{
    using list = type_list<int, double, float>;

    // Transform to add pointer
    using ptr_list = type_transform<list, std::add_pointer_t>;

    static_assert( std::is_same_v<type_at<ptr_list, 0>, int*>, "Should be int*" );
    static_assert( std::is_same_v<type_at<ptr_list, 1>, double*>, "Should be double*" );
    static_assert( std::is_same_v<type_at<ptr_list, 2>, float*>, "Should be float*" );

    // Transform to add const
    using const_list = type_transform<list, std::add_const_t>;

    static_assert( std::is_same_v<type_at<const_list, 0>, const int>, "Should be const int" );
    static_assert( std::is_same_v<type_at<const_list, 1>, const double>, "Should be const double" );

    BOOST_CHECK( true );
}

//
// Test if_t (conditional type selection)
//
BOOST_AUTO_TEST_CASE( test_if_t )
{
    static_assert( std::is_same_v<if_t<true, int, double>, int>,
                   "if_t<true, int, double> should be int" );
    static_assert( std::is_same_v<if_t<false, int, double>, double>,
                   "if_t<false, int, double> should be double" );

    // More complex conditions
    static_assert( std::is_same_v<if_t<(sizeof(int) == 4), int, long>, int>,
                   "Conditional on sizeof should work" );

    static_assert( std::is_same_v<if_t<std::is_integral_v<double>, int, double>, double>,
                   "Conditional on type trait should work" );

    BOOST_CHECK( true );
}

//
// Test constant, int_c, uint16_c, bool_c
//
BOOST_AUTO_TEST_CASE( test_constants )
{
    // int_c
    static_assert( int_c<42>::value == 42, "int_c<42> should have value 42" );
    static_assert( int_c<-1>::value == -1, "int_c<-1> should have value -1" );
    static_assert( int_c<0>::value == 0, "int_c<0> should have value 0" );

    // uint16_c
    static_assert( uint16_c<100>::value == 100, "uint16_c<100> should have value 100" );
    static_assert( uint16_c<0>::value == 0, "uint16_c<0> should have value 0" );

    // bool_c
    static_assert( bool_c<true>::value == true, "bool_c<true> should be true" );
    static_assert( bool_c<false>::value == false, "bool_c<false> should be false" );

    // constant (generic)
    static_assert( meta::constant<42>::value == 42, "constant<42> should have value 42" );
    static_assert( meta::constant<'A'>::value == 'A', "constant<'A'> should have value 'A'" );

    // Check types
    static_assert( std::is_same_v<typename int_c<5>::value_type, int>,
                   "int_c value_type should be int" );
    static_assert( std::is_same_v<typename uint16_c<5>::value_type, uint16_type>,
                   "uint16_c value_type should be uint16_type" );
    static_assert( std::is_same_v<typename bool_c<true>::value_type, bool>,
                   "bool_c value_type should be bool" );

    BOOST_CHECK( true );
}

//
// Test std::type_identity (C++20 standard)
// Note: Feel++ uses std::type_identity directly from C++20
//
BOOST_AUTO_TEST_CASE( test_type_identity )
{
    static_assert( std::is_same_v<typename std::type_identity<int>::type, int>,
                   "std::type_identity<int>::type should be int" );
    static_assert( std::is_same_v<typename std::type_identity<double>::type, double>,
                   "std::type_identity<double>::type should be double" );

    // Works with complex types
    using complex_type = std::vector<std::pair<int, double>>;
    static_assert( std::is_same_v<typename std::type_identity<complex_type>::type, complex_type>,
                   "std::type_identity should work with complex types" );

    BOOST_CHECK( true );
}

//
// Test mp11 integration - mp::mp_push_back, mp::mp_append, etc.
//
BOOST_AUTO_TEST_CASE( test_mp11_integration )
{
    using list1 = type_list<int, double>;
    using list2 = type_list<float, char>;

    // Append element
    using list1_extended = mp::mp_push_back<list1, std::string>;
    static_assert( mp::mp_size<list1_extended>::value == 3, "Extended list should have 3 elements" );
    static_assert( std::is_same_v<type_at<list1_extended, 2>, std::string>,
                   "Last element should be std::string" );

    // Concatenate lists
    using combined = mp::mp_append<list1, list2>;
    static_assert( mp::mp_size<combined>::value == 4, "Combined list should have 4 elements" );
    static_assert( std::is_same_v<type_at<combined, 0>, int>, "First should be int" );
    static_assert( std::is_same_v<type_at<combined, 3>, char>, "Last should be char" );

    // mp_contains
    static_assert( mp::mp_contains<list1, int>::value, "list1 should contain int" );
    static_assert( !mp::mp_contains<list1, char>::value, "list1 should not contain char" );

    BOOST_CHECK( true );
}

// Helper struct for dimension-based type selection test
template<int Dim>
struct DimensionStorage
{
    using type = if_t<(Dim == 1), std::array<double, 2>,
                 if_t<(Dim == 2), std::array<double, 3>,
                      std::array<double, 4>>>;
};

//
// Test practical use case: type selection based on dimension
//
BOOST_AUTO_TEST_CASE( test_practical_type_selection )
{
    static_assert( std::is_same_v<typename DimensionStorage<1>::type, std::array<double, 2>>,
                   "Dim 1 should use array<2>" );
    static_assert( std::is_same_v<typename DimensionStorage<2>::type, std::array<double, 3>>,
                   "Dim 2 should use array<3>" );
    static_assert( std::is_same_v<typename DimensionStorage<3>::type, std::array<double, 4>>,
                   "Dim 3 should use array<4>" );

    BOOST_CHECK( true );
}

// Helper template alias for floating point filter
template<typename T>
using is_floating_point_type = std::bool_constant<std::is_floating_point_v<T>>;

//
// Test type list iteration pattern
//
BOOST_AUTO_TEST_CASE( test_type_list_iteration )
{
    using types = type_list<int, double, float>;

    // Filter to get floating point types only
    using floating_types = mp::mp_filter<is_floating_point_type, types>;
    static_assert( mp::mp_size<floating_types>::value == 2,
                   "Should have 2 floating point types" );

    // mp_for_each can be used at runtime
    int count = 0;
    mp::mp_for_each<types>( [&count]( auto t ) {
        using T = std::decay_t<decltype(t)>;
        if constexpr ( std::is_floating_point_v<T> )
            ++count;
    } );

    BOOST_CHECK_EQUAL( count, 2 );
}

BOOST_AUTO_TEST_SUITE_END()
