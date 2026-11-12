/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_eigenmap.cpp
 * @brief Zero-copy, constness, layout, and serialization tests for Eigen views.
 */
#define BOOST_TEST_MODULE test_eigenmap
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/eigenmap.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <sstream>
#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( eigenmap_suite )

/** @test Verify concepts, aliasing, constness, and column-major strides. */
BOOST_AUTO_TEST_CASE( column_major_matrix_view )
{
    using matrix_type =
        ublas::matrix<double, ublas::column_major>;
    static_assert( UblasDenseMatrix<matrix_type> );
    static_assert( !UblasDenseVector<matrix_type> );

    matrix_type matrix( 2, 3 );
    static_assert( noexcept( eigenMap( matrix ) ) );
    for ( std::size_t j = 0; j < matrix.size2(); ++j )
        for ( std::size_t i = 0; i < matrix.size1(); ++i )
            matrix( i, j ) = 10.0 * j + i;

    auto view = eigenMap( matrix );
    static_assert( !decltype( view )::IsRowMajor );
    static_assert( std::same_as<decltype( view.data() ), double*> );
    BOOST_CHECK_EQUAL( view.data(),
                       std::to_address( matrix.data().begin() ) );
    BOOST_CHECK_EQUAL( view.innerStride(), 1 );
    BOOST_CHECK_EQUAL( view.outerStride(), 2 );
    view( 1, 2 ) = 42.0;
    BOOST_CHECK_EQUAL( matrix( 1, 2 ), 42.0 );

    matrix_type const& constMatrix = matrix;
    auto constView = eigenMap( constMatrix );
    static_assert(
        std::same_as<decltype( constView.data() ), double const*> );
    BOOST_CHECK_EQUAL( constView( 1, 2 ), 42.0 );
}

/** @test Verify row-major layout and bounded contiguous storage. */
BOOST_AUTO_TEST_CASE( row_major_and_bounded_storage )
{
    using row_matrix_type =
        ublas::matrix<double, ublas::row_major>;
    row_matrix_type rowMatrix( 2, 3 );
    auto rowView = eigenMap( rowMatrix );
    static_assert( decltype( rowView )::IsRowMajor );
    BOOST_CHECK_EQUAL( rowView.innerStride(), 1 );
    BOOST_CHECK_EQUAL( rowView.outerStride(), 3 );
    rowView( 1, 2 ) = 7.5;
    BOOST_CHECK_EQUAL( rowMatrix( 1, 2 ), 7.5 );

    using bounded_matrix_type = ublas::matrix<
        double, ublas::column_major, ublas::bounded_array<double, 6>>;
    bounded_matrix_type boundedMatrix( 2, 3 );
    auto boundedView = eigenMap( boundedMatrix );
    boundedView.setConstant( 3.25 );
    BOOST_CHECK_EQUAL( boundedMatrix( 0, 0 ), 3.25 );
    BOOST_CHECK_EQUAL( boundedMatrix( 1, 2 ), 3.25 );
}

/** @test Verify vector semantics, constness, and empty-storage handling. */
BOOST_AUTO_TEST_CASE( vector_and_empty_views )
{
    using vector_type = ublas::vector<double>;
    static_assert( UblasDenseVector<vector_type> );
    static_assert( !UblasDenseMatrix<vector_type> );

    vector_type vector( 3 );
    auto view = eigenMap( vector );
    static_assert( decltype( view )::IsVectorAtCompileTime );
    view << 1.0, 2.0, 3.0;
    BOOST_CHECK_EQUAL( vector( 2 ), 3.0 );

    vector_type const& constVector = vector;
    auto constView = eigenMap( constVector );
    static_assert(
        std::same_as<decltype( constView.data() ), double const*> );
    BOOST_CHECK_EQUAL( constView.sum(), 6.0 );

    vector_type emptyVector;
    BOOST_CHECK_EQUAL( eigenMap( emptyVector ).size(), 0 );
    ublas::matrix<double, ublas::column_major> emptyMatrix( 0, 3 );
    auto emptyMatrixView = eigenMap( emptyMatrix );
    BOOST_CHECK_EQUAL( emptyMatrixView.rows(), 0 );
    BOOST_CHECK_EQUAL( emptyMatrixView.cols(), 3 );
}

/** @test Verify serialization remains attached to the owning uBLAS boundary. */
BOOST_AUTO_TEST_CASE( owning_storage_serialization_round_trip )
{
    using matrix_type =
        ublas::matrix<double, ublas::column_major>;
    matrix_type source( 2, 2 );
    eigenMap( source ) << 1.0, 2.0, 3.0, 4.0;

    std::stringstream stream;
    {
        boost::archive::text_oarchive output( stream );
        output << source;
    }

    matrix_type restored;
    {
        boost::archive::text_iarchive input( stream );
        input >> restored;
    }

    BOOST_REQUIRE_EQUAL( restored.size1(), 2 );
    BOOST_REQUIRE_EQUAL( restored.size2(), 2 );
    BOOST_CHECK_SMALL( ( eigenMap( source ) - eigenMap( restored ) ).norm(),
                       1e-14 );
}

BOOST_AUTO_TEST_SUITE_END()
