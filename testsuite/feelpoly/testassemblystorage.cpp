/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file testassemblystorage.cpp
 * @brief Independent mixed extent, mortar, two-cell and orientation storage contracts.
 */
#define BOOST_TEST_MODULE assemblystorage
#if defined(FEELPP_STORAGE_NDEBUG) && !defined(NDEBUG)
#error "storage release checks require NDEBUG to be defined"
#endif
#if defined(FEELPP_STORAGE_ASSERTIONS) && defined(NDEBUG)
#error "storage assertion checks require NDEBUG to be undefined"
#endif
#include <boost/mp11/list.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelvf/assemblystorage.hpp>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief One static/runtime test and trial axis combination. */
template <int Rows, int Cols>
struct Shape
{
    /** Fixed extent or runtime sentinel for the test side. */
    static constexpr int rows = Rows;
    /** Fixed extent or runtime sentinel for the trial side. */
    static constexpr int cols = Cols;
};
/** All combinations plus P0 and an empty mortar-compatible extent. */
using Shapes = boost::mp11::mp_list<Shape<6, 4>, Shape<6, Dynamic>, Shape<Dynamic, 4>, Shape<Dynamic, Dynamic>,
                                    Shape<1, 1>, Shape<1, 2>, Shape<2, 1>, Shape<0, 1>>;

/** @brief Assemble independent scalar products and verify shape-preserving storage reuse. */
template <typename Matrix>
void checkMatrix( int rows, int cols )
{
    Matrix matrix;
    detail::resizeAssemblyMatrix( matrix, rows, cols );
    Eigen::MatrixXd left( rows, 5 ), right( cols, 5 );
    for ( int i = 0; i < rows; ++i )
        for ( int q = 0; q < 5; ++q )
            left( i, q ) = .2 + i * .1 + q * .03;
    for ( int j = 0; j < cols; ++j )
        for ( int q = 0; q < 5; ++q )
            right( j, q ) = .3 + j * .07 - q * .02;
    // Eigen's GEMM implementation cannot instantiate a statically empty destination.
    if constexpr ( Matrix::SizeAtCompileTime != 0 )
        if ( rows != 0 && cols != 0 ) matrix.noalias() = left * right.transpose();
    for ( int i = 0; i < rows; ++i )
        for ( int j = 0; j < cols; ++j )
        {
            double expected = 0;
            for ( int q = 0; q < 5; ++q )
                expected += left( i, q ) * right( j, q );
            BOOST_CHECK_SMALL( matrix( i, j ) - expected, 1e-12 );
        }
    auto const* address = matrix.data();
    detail::resizeAssemblyMatrix( matrix, rows, cols );
    BOOST_CHECK( matrix.data() == address );
}
} // namespace

BOOST_AUTO_TEST_CASE_TEMPLATE( mixedAssemblyShapes, Case, Shapes )
{
    using Storage = detail::AssemblyStorage<double, Case::rows, Case::cols>;
    constexpr int rows = Case::rows == Dynamic ? 6 : Case::rows;
    constexpr int cols = Case::cols == Dynamic ? 4 : Case::cols;
    static_assert( Storage::Ordinary::RowsAtCompileTime == Case::rows );
    static_assert( Storage::Ordinary::ColsAtCompileTime == Case::cols );
    checkMatrix<typename Storage::Ordinary>( rows, cols );
    checkMatrix<typename Storage::MortarTest>( std::max( rows - 1, 0 ), cols );
    checkMatrix<typename Storage::MortarTrial>( rows, std::max( cols - 1, 0 ) );
    checkMatrix<typename Storage::TwoCell>( 2 * rows, 2 * cols );
    typename Storage::RowSigns rowSigns;
    typename Storage::ColSigns colSigns;
    detail::resizeAssemblyMatrix( rowSigns, rows );
    detail::resizeAssemblyMatrix( colSigns, cols );
    rowSigns.setOnes();
    colSigns.setConstant( -1 );
    BOOST_CHECK_EQUAL( rowSigns.sum(), rows );
    BOOST_CHECK_EQUAL( colSigns.sum(), -cols );
}

BOOST_AUTO_TEST_CASE( allocationAndOverflowContracts )
{
    using Large = detail::MatrixStoragePolicy<double, 256, 256, detail::assemblyInlineMatrixBytes>;
    static_assert( !Large::usesInlineStorage );
    static_assert( Large::rowsAtCompileTime == 256 && Large::colsAtCompileTime == 256 );
    static_assert( sizeof( Large::matrix_type ) < 128 );
    static_assert( detail::AssemblyStorage<double, Dynamic, 4>::Ordinary::ColsAtCompileTime == 4 );
    static_assert( detail::AssemblyStorage<double, 6, Dynamic>::Ordinary::RowsAtCompileTime == 6 );
    static_assert( detail::assemblyExtent( Dynamic, 1 ) == Dynamic );
    static_assert( detail::assemblyExtent( 0, 1 ) == 0 );
    static_assert( detail::assemblyExtent( 1, 1 ) == 0 );
    static_assert( detail::assemblyExtent( 3, 0, 2 ) == 6 );
    BOOST_CHECK_THROW( static_cast<void>( detail::assemblyExtent( std::numeric_limits<int>::max(), 0, 2 ) ), std::length_error );
    BOOST_CHECK_THROW( static_cast<void>( detail::assemblyDofCount( std::numeric_limits<std::size_t>::max() ) ), std::length_error );
    BOOST_CHECK_THROW( static_cast<void>( detail::assemblyExtent( -2 ) ), std::invalid_argument );
    Eigen::Matrix<double, 2, Dynamic> mixed;
    BOOST_CHECK_THROW( detail::resizeAssemblyMatrix( mixed, 3, 4 ), std::invalid_argument );
    BOOST_CHECK_THROW( detail::resizeAssemblyMatrix( mixed, 2, -1 ), std::invalid_argument );
    checkMatrix<Large::matrix_type>( 256, 256 );
}

BOOST_AUTO_TEST_CASE( independent_static_and_runtime_local_axes )
{
    using Storage = detail::AssemblyStorage<double, 6, Dynamic>;
    static_assert( Storage::Ordinary::RowsAtCompileTime == 6 );
    static_assert( Storage::Ordinary::ColsAtCompileTime == Dynamic );
    Storage::Ordinary matrix;
    detail::resizeAssemblyMatrix( matrix, 6, 10 );
    matrix.setConstant( 2.0 );
    BOOST_CHECK_EQUAL( matrix.rows(), 6 );
    BOOST_CHECK_EQUAL( matrix.cols(), 10 );
    BOOST_CHECK_THROW( detail::resizeAssemblyMatrix( matrix, 7, 10 ), std::invalid_argument );
    BOOST_CHECK_EQUAL( matrix.rows(), 6 );
    BOOST_CHECK_EQUAL( matrix.cols(), 10 );
    BOOST_CHECK_EQUAL( matrix.sum(), 120.0 );
    checkMatrix<Storage::Ordinary>( 6, 10 );

    using TransposedStorage = detail::AssemblyStorage<double, Dynamic, 10>;
    static_assert( TransposedStorage::Ordinary::ColsAtCompileTime == 10 );
    TransposedStorage::Ordinary transposed;
    detail::resizeAssemblyMatrix( transposed, 6, 10 );
    BOOST_CHECK_THROW( detail::resizeAssemblyMatrix( transposed, 6, 11 ), std::invalid_argument );

    detail::EigenMatrixBuffer<double, 6, Dynamic> bound( 6, 10, 6, 10 );
    BOOST_CHECK_THROW( bound.resize( 7, 10 ), std::invalid_argument );
    BOOST_CHECK_THROW( bound.resize( 6, 11 ), std::invalid_argument );
    detail::EigenMatrixBuffer<double, 256, 256> large;
    static_assert( !decltype( large )::usesInlineStorage );
    BOOST_CHECK_THROW( large.resize( 256, 257 ), std::invalid_argument );
}
