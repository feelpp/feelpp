/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 Feel++ Consortium
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/**
   \file test_policy.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-04
   \brief Tests for feelpoly policy helpers and storage types
 */

#define BOOST_TEST_MODULE test_policy
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/policy.hpp>

#include <Eigen/Core>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <cstddef>
#include <boost/mp11.hpp>
#include <array>
#include <type_traits>

namespace ublas = boost::numeric::ublas;
using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief Borrow either numerical representation without copying coefficients. */
template<typename Matrix>
decltype(auto) matrixView( Matrix& matrix )
{
    if constexpr ( EigenMatrix<Matrix> ) return ( matrix );
    else return eigenMap( matrix );
}

/** @brief Fill either representation in the same mathematical row/column order. */
template<typename Matrix>
void fillMatrix( Matrix& matrix )
{
    auto&& view = matrixView( matrix );
    typename std::remove_cvref_t<decltype( view )>::Scalar value = 1;
    for ( Eigen::Index i = 0; i < view.rows(); ++i )
        for ( Eigen::Index j = 0; j < view.cols(); ++j ) view( i, j ) = value++;
}

/** @brief Share exact reshape/roundtrip assertions across both representations. */
template<typename MatrixA, typename MatrixB>
void checkMatricesEqual( MatrixA const& a, MatrixB const& b )
{
    auto&& left = matrixView( a );
    auto&& right = matrixView( b );
    BOOST_REQUIRE_EQUAL( left.rows(), right.rows() );
    BOOST_REQUIRE_EQUAL( left.cols(), right.cols() );
    for ( Eigen::Index i = 0; i < left.rows(); ++i )
        for ( Eigen::Index j = 0; j < left.cols(); ++j ) BOOST_CHECK_EQUAL( left( i, j ), right( i, j ) );
}
} // namespace

BOOST_AUTO_TEST_SUITE( policyContracts )

BOOST_AUTO_TEST_CASE( storageEigenTypes )
{
    using storage = StorageEigen<double>;
    static_assert( std::is_same_v<storage::matrix_type,
                                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> );
    static_assert( std::is_same_v<storage::vector_type,
                                  Eigen::Matrix<double, Eigen::Dynamic, 1>> );
    static_assert( std::is_same_v<storage::vector_matrix_type,
                                  std::vector<storage::matrix_type>> );
    static_assert( std::is_same_v<storage::node_type,
                                  Eigen::Matrix<double, Eigen::Dynamic, 1>> );
    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( rankHelpers )
{
    static_assert( std::is_same_v<typename RankUp<Scalar<2>>::type, Vectorial<2>> );
    static_assert( std::is_same_v<typename RankUp2<Scalar<2>>::type, Tensor2<2>> );
    static_assert( std::is_same_v<typename RankDown<Tensor2<2>>::type, Vectorial<2>> );
    static_assert( std::is_same_v<typename RankDown<Tensor2Symm<2>>::type, Vectorial<2>> );
    static_assert( std::is_same_v<typename RankDown2<Tensor3<2>>::type, Vectorial<2>> );
    static_assert( std::is_same_v<typename RankDown<Scalar<2>>::type, Scalar<2>> );
    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( scalarCoefficientIdentity )
{
    ublas::matrix<double> u( 2, 2 );
    fillMatrix( u );
    auto const& u_mat = Scalar<2>::toMatrix( u );
    checkMatricesEqual( u, u_mat );
    auto const& u_type = Scalar<2>::toType( u );
    checkMatricesEqual( u, u_type );

    Eigen::Matrix<double, 2, 2> e;
    fillMatrix( e );
    auto e_mat = Scalar<2>::toMatrix( e );
    checkMatricesEqual( e, e_mat );
    auto e_type = Scalar<2>::toType( e );
    checkMatricesEqual( e, e_type );
}

/** @test The same packing rule supports arbitrary polynomial counts in both representations. */
using ComponentShapes = boost::mp11::mp_list<Vectorial<2>, Tensor2<2>, Tensor2Symm<2>>;
BOOST_AUTO_TEST_CASE_TEMPLATE( componentCoefficientPacking, Shape, ComponentShapes )
{
    constexpr Eigen::Index components = Shape::nComponents;
    for ( Eigen::Index functions : std::array<Eigen::Index, 4>{ 0, 1, 2, 7 } )
        for ( Eigen::Index modes : std::array<Eigen::Index, 3>{ 0, 1, 3 } )
        {
            Eigen::MatrixXd native( components * functions, modes );
            fillMatrix( native );
            ublas::matrix<double> legacy( native.rows(), native.cols() );
            eigenMap( legacy ) = native;
            auto packed = Shape::toMatrix( native );
            auto legacyPacked = Shape::toMatrix( legacy );
            BOOST_REQUIRE_EQUAL( packed.rows(), functions );
            BOOST_REQUIRE_EQUAL( packed.cols(), components * modes );
            checkMatricesEqual( packed, legacyPacked );
            for ( Eigen::Index function = 0; function < functions; ++function )
                for ( Eigen::Index component = 0; component < components; ++component )
                    for ( Eigen::Index mode = 0; mode < modes; ++mode )
                        BOOST_CHECK_EQUAL( packed( function, component * modes + mode ),
                                           native( function * components + component, mode ) );
            checkMatricesEqual( Shape::toType( packed ), native );
            checkMatricesEqual( Shape::toType( legacyPacked ), legacy );
        }
    Eigen::MatrixXd invalidRows = Eigen::MatrixXd::Zero( components + 1, 2 );
    Eigen::MatrixXd invalidCols = Eigen::MatrixXd::Zero( 2, components + 1 );
    BOOST_CHECK_THROW( (void)Shape::toMatrix( invalidRows ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)Shape::toType( invalidCols ), std::invalid_argument );
}

/** @test Scalar identities borrow lvalues and own evaluated temporary expressions. */
BOOST_AUTO_TEST_CASE( coefficientIdentityLifetimes )
{
    Eigen::Matrix2d input = Eigen::Matrix2d::Identity();
    static_assert( std::is_lvalue_reference_v<decltype( Scalar<2>::toMatrix( input ) )> );
    auto temporary = Scalar<2>::toMatrix( input + input );
    static_assert( decltype( temporary )::RowsAtCompileTime == 2 );
    input.setZero();
    BOOST_CHECK_EQUAL( temporary( 0, 0 ), 2. );
    auto nested = Scalar<2>::toType( Eigen::Matrix2d::Identity() + Eigen::Matrix2d::Identity() );
    BOOST_CHECK_EQUAL( nested( 1, 1 ), 2. );
    static_assert( !std::is_reference_v<decltype( Scalar<2>::toMatrix( Eigen::Matrix2d{} ) )> );
    auto const& ownedTemporary = Scalar<2>::toMatrix( Eigen::Matrix2d::Identity().eval() );
    BOOST_CHECK_EQUAL( ownedTemporary( 0, 0 ), 1. );
}

/** @test Component packing retains known axes independently of layout and allocation. */
BOOST_AUTO_TEST_CASE( fixedComponentExtents )
{
    Eigen::Matrix<double, 8, 3> input;
    fillMatrix( input );
    auto vector = Vectorial<2>::toMatrix( input );
    auto tensor = Tensor2<2>::toMatrix( input );
    static_assert( decltype( vector )::RowsAtCompileTime == 4 );
    static_assert( decltype( vector )::ColsAtCompileTime == 6 );
    static_assert( decltype( tensor )::RowsAtCompileTime == 2 );
    static_assert( decltype( tensor )::ColsAtCompileTime == 12 );
    checkMatricesEqual( Tensor2<2>::toType( tensor ), input );
}

BOOST_AUTO_TEST_SUITE_END()
