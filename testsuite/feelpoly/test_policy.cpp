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
#include <type_traits>

namespace ublas = boost::numeric::ublas;
using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template<typename Matrix>
void fill_ublas( Matrix& m )
{
    double v = 1.0;
    for ( std::size_t i = 0; i < m.size1(); ++i )
        for ( std::size_t j = 0; j < m.size2(); ++j )
            m( i, j ) = v++;
}

template<typename Matrix>
void check_ublas_equal( Matrix const& a, Matrix const& b )
{
    BOOST_REQUIRE_EQUAL( a.size1(), b.size1() );
    BOOST_REQUIRE_EQUAL( a.size2(), b.size2() );
    for ( std::size_t i = 0; i < a.size1(); ++i )
        for ( std::size_t j = 0; j < a.size2(); ++j )
            BOOST_CHECK_EQUAL( a( i, j ), b( i, j ) );
}

template<typename Matrix>
void fill_eigen( Matrix& m )
{
    typename Matrix::Scalar v = 1;
    for ( Eigen::Index i = 0; i < m.rows(); ++i )
        for ( Eigen::Index j = 0; j < m.cols(); ++j )
            m( i, j ) = v++;
}

template<typename MatrixA, typename MatrixB>
void check_eigen_equal( MatrixA const& a, MatrixB const& b )
{
    BOOST_REQUIRE_EQUAL( a.rows(), b.rows() );
    BOOST_REQUIRE_EQUAL( a.cols(), b.cols() );
    for ( Eigen::Index i = 0; i < a.rows(); ++i )
        for ( Eigen::Index j = 0; j < a.cols(); ++j )
            BOOST_CHECK_EQUAL( a( i, j ), b( i, j ) );
}
} // namespace

BOOST_AUTO_TEST_SUITE( policy_suite )

BOOST_AUTO_TEST_CASE( test_storage_eigen_types )
{
    using storage = StorageEigen<double>;
    static_assert( std::is_same_v<storage::matrix_type,
                                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> );
    static_assert( std::is_same_v<storage::vector_type,
                                  Eigen::Matrix<double, Eigen::Dynamic, 1>> );
    static_assert( std::is_same_v<storage::vector_matrix_type,
                                  ublas::vector<storage::matrix_type>> );
    static_assert( std::is_same_v<storage::node_type,
                                  Eigen::Matrix<double, Eigen::Dynamic, 1>> );
    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_rank_helpers )
{
    static_assert( std::is_same_v<typename RankUp<Scalar<2>>::type, Vectorial<2>> );
    static_assert( std::is_same_v<typename RankUp2<Scalar<2>>::type, Tensor2<2>> );
    static_assert( std::is_same_v<typename RankDown<Tensor2<2>>::type, Vectorial<2>> );
    static_assert( std::is_same_v<typename RankDown<Tensor2Symm<2>>::type, Vectorial<2>> );
    static_assert( std::is_same_v<typename RankDown2<Tensor3<2>>::type, Vectorial<2>> );
    static_assert( std::is_same_v<typename RankDown<Scalar<2>>::type, Scalar<2>> );
    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_scalar_to_matrix )
{
    ublas::matrix<double> u( 2, 2 );
    fill_ublas( u );
    auto const& u_mat = Scalar<2>::toMatrix( u );
    check_ublas_equal( u, u_mat );
    auto const& u_type = Scalar<2>::toType( u );
    check_ublas_equal( u, u_type );

    Eigen::Matrix<double, 2, 2> e;
    fill_eigen( e );
    auto e_mat = Scalar<2>::toMatrix( e );
    check_eigen_equal( e, e_mat );
    auto e_type = Scalar<2>::toType( e );
    check_eigen_equal( e, e_type );
}

BOOST_AUTO_TEST_CASE( test_vectorial_roundtrip_ublas )
{
    const std::size_t comp = Vectorial<2>::nComponents;
    ublas::matrix<double> input( comp * 2, 3 );
    fill_ublas( input );

    auto reshaped = Vectorial<2>::toMatrix( input );
    BOOST_CHECK_EQUAL( reshaped.size1(), input.size1() / comp );
    BOOST_CHECK_EQUAL( reshaped.size2(), input.size2() * comp );

    auto roundtrip = Vectorial<2>::toType( reshaped );
    check_ublas_equal( input, roundtrip );
}

BOOST_AUTO_TEST_CASE( test_vectorial_roundtrip_eigen )
{
    const Eigen::Index comp = static_cast<Eigen::Index>( Vectorial<2>::nComponents );
    Eigen::MatrixXd input( comp * 2, 3 );
    fill_eigen( input );

    auto reshaped = Vectorial<2>::toMatrix( input );
    BOOST_CHECK_EQUAL( reshaped.rows(), input.rows() / comp );
    BOOST_CHECK_EQUAL( reshaped.cols(), input.cols() * comp );

    auto roundtrip = Vectorial<2>::toType( reshaped );
    check_eigen_equal( input, roundtrip );
}

BOOST_AUTO_TEST_CASE( test_tensor2_roundtrip_ublas )
{
    const std::size_t comp = Tensor2<2>::nComponents;
    const std::size_t nCols = 1;
    ublas::matrix<double> input( nCols * comp * comp, nCols );
    fill_ublas( input );

    auto reshaped = Tensor2<2>::toMatrix( input );
    BOOST_CHECK_EQUAL( reshaped.size1(), input.size1() / comp );
    BOOST_CHECK_EQUAL( reshaped.size2(), input.size2() * comp );

    auto roundtrip = Tensor2<2>::toType( reshaped );
    check_ublas_equal( input, roundtrip );
}

BOOST_AUTO_TEST_CASE( test_tensor2_roundtrip_eigen )
{
    const Eigen::Index comp = static_cast<Eigen::Index>( Tensor2<2>::nComponents );
    const Eigen::Index nCols = 1;
    Eigen::MatrixXd input( nCols * comp * comp, nCols );
    fill_eigen( input );

    auto reshaped = Tensor2<2>::toMatrix( input );
    BOOST_CHECK_EQUAL( reshaped.rows(), input.rows() / comp );
    BOOST_CHECK_EQUAL( reshaped.cols(), input.cols() * comp );

    auto roundtrip = Tensor2<2>::toType( reshaped );
    check_eigen_equal( input, roundtrip );
}

BOOST_AUTO_TEST_SUITE_END()
