/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE test_tensorformat
#include <feel/feelcore/testsuite.hpp>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <feel/feeldiscr/tensorformat.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( tensorformat )

static_assert( SymmetricTensorFormatLike<SymmetricTensorFormat> );

struct LocalFormat
{
    SymmetricTensorOrder order = SymmetricTensorOrder::DiagonalFirst;
    SymmetricTensorScaling scaling = SymmetricTensorScaling::Mandel;
};

static_assert( SymmetricTensorFormatLike<LocalFormat> );

template<std::size_t N>
std::vector<double>
reordered( std::array<double,N> const& storage, uint16_type n, SymmetricTensorOrder order )
{
    std::vector<double> out( symmetricTensorComponentCount( n, order ), 0. );
    for ( uint16_type storageSlot = 0; storageSlot < symmetricTensorStorageComponentCount( n ); ++storageSlot )
        out[symmetricTensorOutputSlotFromStorage( storageSlot, n, order )] = storage[storageSlot];
    return out;
}

void
check_values( std::vector<double> const& values, std::initializer_list<double> expected )
{
    BOOST_REQUIRE_EQUAL( values.size(), expected.size() );
    auto it = expected.begin();
    for ( std::size_t k = 0; k < values.size(); ++k, ++it )
        BOOST_CHECK_CLOSE( values[k], *it, 1e-12 );
}

void
check_map( std::array<uint16_type,6> const& values, std::array<uint16_type,6> const& expected )
{
    for ( std::size_t k = 0; k < values.size(); ++k )
        BOOST_CHECK_EQUAL( values[k], expected[k] );
}

BOOST_AUTO_TEST_CASE( storage_indices_are_legacy_symmetric_indices )
{
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 0, 0, 2 ), 0 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 0, 1, 2 ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 1, 0, 2 ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 1, 1, 2 ), 2 );

    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 0, 0, 3 ), 0 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 0, 1, 3 ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 1, 0, 3 ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 0, 2, 3 ), 2 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 2, 0, 3 ), 2 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 1, 1, 3 ), 3 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 1, 2, 3 ), 4 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 2, 1, 3 ), 4 );
    BOOST_CHECK_EQUAL( symmetricTensorStorageIndex( 2, 2, 3 ), 5 );
}

BOOST_AUTO_TEST_CASE( mappings_3d )
{
    std::array<double,6> storage = {{ 11., 12., 13., 22., 23., 33. }};

    check_values( reordered( storage, 3, SymmetricTensorOrder::Storage ),
                  { 11., 12., 13., 22., 23., 33. } );
    check_values( reordered( storage, 3, SymmetricTensorOrder::DiagonalFirst ),
                  { 11., 22., 33., 12., 13., 23. } );
    check_values( reordered( storage, 3, SymmetricTensorOrder::VtkTensor6 ),
                  { 11., 22., 33., 12., 23., 13. } );
    check_values( reordered( storage, 3, SymmetricTensorOrder::EnsightTensor6 ),
                  { 11., 22., 33., 12., 13., 23. } );
    check_values( reordered( storage, 3, SymmetricTensorOrder::XdmfTensor6 ),
                  { 11., 12., 13., 22., 23., 33. } );

    auto vtkMap = symmetricTensorStorageToOutputMap( 3, SymmetricTensorOrder::VtkTensor6 );
    BOOST_CHECK_EQUAL( vtkMap[0], 0 );
    BOOST_CHECK_EQUAL( vtkMap[1], 3 );
    BOOST_CHECK_EQUAL( vtkMap[2], 5 );
    BOOST_CHECK_EQUAL( vtkMap[3], 1 );
    BOOST_CHECK_EQUAL( vtkMap[4], 4 );
    BOOST_CHECK_EQUAL( vtkMap[5], 2 );

    auto ensightMap = symmetricTensorStorageToOutputMap( 3, SymmetricTensorOrder::EnsightTensor6 );
    BOOST_CHECK_EQUAL( ensightMap[0], 0 );
    BOOST_CHECK_EQUAL( ensightMap[1], 3 );
    BOOST_CHECK_EQUAL( ensightMap[2], 4 );
    BOOST_CHECK_EQUAL( ensightMap[3], 1 );
    BOOST_CHECK_EQUAL( ensightMap[4], 5 );
    BOOST_CHECK_EQUAL( ensightMap[5], 2 );
}

BOOST_AUTO_TEST_CASE( exporter_tensor6_maps_match_legacy_arrays )
{
    check_map( symmetricTensorStorageToOutputMap( 3, SymmetricTensorOrder::VtkTensor6 ),
               std::array<uint16_type,6>{{ 0, 3, 5, 1, 4, 2 }} );
    check_map( symmetricTensorStorageToOutputMap( 3, SymmetricTensorOrder::EnsightTensor6 ),
               std::array<uint16_type,6>{{ 0, 3, 4, 1, 5, 2 }} );
    check_map( symmetricTensorStorageToOutputMap( 3, SymmetricTensorOrder::XdmfTensor6 ),
               std::array<uint16_type,6>{{ 0, 1, 2, 3, 4, 5 }} );

    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 0, 0, 3, SymmetricTensorOrder::VtkTensor6 ), 0 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 1, 1, 3, SymmetricTensorOrder::VtkTensor6 ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 0, 1, 3, SymmetricTensorOrder::VtkTensor6 ), 3 );

    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 0, 0, 3, SymmetricTensorOrder::EnsightTensor6 ), 0 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 1, 1, 3, SymmetricTensorOrder::EnsightTensor6 ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 0, 1, 3, SymmetricTensorOrder::EnsightTensor6 ), 3 );

    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 0, 0, 3, SymmetricTensorOrder::XdmfTensor6 ), 0 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 0, 1, 3, SymmetricTensorOrder::XdmfTensor6 ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 1, 1, 3, SymmetricTensorOrder::XdmfTensor6 ), 3 );
}

BOOST_AUTO_TEST_CASE( mappings_2d )
{
    std::array<double,3> storage = {{ 11., 12., 22. }};

    check_values( reordered( storage, 2, SymmetricTensorOrder::Storage ),
                  { 11., 12., 22. } );
    check_values( reordered( storage, 2, SymmetricTensorOrder::DiagonalFirst ),
                  { 11., 22., 12. } );
    check_values( reordered( storage, 2, SymmetricTensorOrder::VtkTensor6 ),
                  { 11., 22., 0., 12., 0., 0. } );
    check_values( reordered( storage, 2, SymmetricTensorOrder::EnsightTensor6 ),
                  { 11., 22., 0., 12., 0., 0. } );
    check_values( reordered( storage, 2, SymmetricTensorOrder::XdmfTensor6 ),
                  { 11., 12., 0., 22., 0., 0. } );
}

BOOST_AUTO_TEST_CASE( output_slots_from_components )
{
    LocalFormat format;

    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 2, 0, 3, SymmetricTensorOrder::VtkTensor6 ), 5 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 0, 2, 3, SymmetricTensorOrder::EnsightTensor6 ), 4 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 1, 1, 3, SymmetricTensorOrder::DiagonalFirst ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 1, 1, 3, format ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 1, 1, 2, SymmetricTensorOrder::VtkTensor6 ), 1 );
    BOOST_CHECK_EQUAL( symmetricTensorOutputSlot( 0, 1, 2, SymmetricTensorOrder::XdmfTensor6 ), 1 );
}

BOOST_AUTO_TEST_CASE( component_labels )
{
    BOOST_CHECK_EQUAL( std::string( symmetricTensorOrderName( SymmetricTensorOrder::Storage ) ), "Storage" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorOrderName( SymmetricTensorOrder::DiagonalFirst ) ), "DiagonalFirst" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorOrderName( SymmetricTensorOrder::VtkTensor6 ) ), "VtkTensor6" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorScalingName( SymmetricTensorScaling::Tensor ) ), "Tensor" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorScalingName( SymmetricTensorScaling::EngineeringShear ) ), "EngineeringShear" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorScalingName( SymmetricTensorScaling::Mandel ) ), "Mandel" );

    BOOST_CHECK_EQUAL( std::string( symmetricTensorComponentLabel( 0, 3, SymmetricTensorOrder::Storage ) ), "xx" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorComponentLabel( 2, 3, SymmetricTensorOrder::Storage ) ), "xz" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorComponentLabel( 3, 3, SymmetricTensorOrder::Storage ) ), "yy" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorComponentLabel( 2, 2, SymmetricTensorOrder::DiagonalFirst ) ), "xy" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorComponentLabel( 4, 3, SymmetricTensorOrder::VtkTensor6 ) ), "yz" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorComponentLabel( 4, 3, SymmetricTensorOrder::EnsightTensor6 ) ), "xz" );
    BOOST_CHECK_EQUAL( std::string( symmetricTensorComponentLabel( 3, 2, SymmetricTensorOrder::XdmfTensor6 ) ), "yy" );
}

BOOST_AUTO_TEST_CASE( scaling )
{
    BOOST_CHECK_CLOSE( symmetricTensorComponentScale( 0, 0, 3, SymmetricTensorScaling::Tensor ), 1., 1e-12 );
    BOOST_CHECK_CLOSE( symmetricTensorComponentScale( 0, 1, 3, SymmetricTensorScaling::Tensor ), 1., 1e-12 );
    BOOST_CHECK_CLOSE( symmetricTensorComponentScale( 0, 1, 3, SymmetricTensorScaling::EngineeringShear ), 2., 1e-12 );
    BOOST_CHECK_CLOSE( symmetricTensorComponentScale( 0, 1, 3, SymmetricTensorScaling::Mandel ), std::sqrt( 2. ), 1e-12 );

    BOOST_CHECK_CLOSE( symmetricTensorStorageScale( 0, 3, SymmetricTensorScaling::EngineeringShear ), 1., 1e-12 );
    BOOST_CHECK_CLOSE( symmetricTensorStorageScale( 1, 3, SymmetricTensorScaling::EngineeringShear ), 2., 1e-12 );
    BOOST_CHECK_CLOSE( symmetricTensorStorageScale( 2, 3, SymmetricTensorScaling::Mandel ), std::sqrt( 2. ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( invalid_inputs_fail_clearly )
{
    BOOST_CHECK_THROW( symmetricTensorComponentCount( 4, SymmetricTensorOrder::Storage ), std::invalid_argument );
    BOOST_CHECK_THROW( symmetricTensorStorageIndex( 2, 0, 2 ), std::invalid_argument );
    BOOST_CHECK_THROW( symmetricTensorOutputSlotFromStorage( 6, 3, SymmetricTensorOrder::Storage ), std::invalid_argument );
    BOOST_CHECK_THROW( symmetricTensorComponentLabel( 3, 2, SymmetricTensorOrder::DiagonalFirst ), std::invalid_argument );
    BOOST_CHECK_THROW( symmetricTensorStorageScale( 3, 2, SymmetricTensorScaling::Tensor ), std::invalid_argument );
}

BOOST_AUTO_TEST_SUITE_END()
