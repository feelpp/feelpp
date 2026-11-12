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
   \file test_concepts.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-04
   \brief Tests for C++20 concepts in feelpoly
 */

#define BOOST_TEST_MODULE test_concepts
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>

#include <Eigen/Core>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;
using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
struct WithStaticOrder
{
    static constexpr int nOrder = 5;
};

struct WithOrderMethod
{
    int order() const { return 3; }
};
} // namespace

BOOST_AUTO_TEST_SUITE( concepts_suite )

//
// Test Dynamic constant
//
BOOST_AUTO_TEST_CASE( test_dynamic_constant )
{
    BOOST_CHECK_EQUAL( Feel::Dynamic, -1 );
    static_assert( Feel::Dynamic == -1, "Dynamic should be -1" );
}

//
// Test EigenMatrix and EigenVector concepts
//
BOOST_AUTO_TEST_CASE( test_eigen_concepts )
{
    // Static matrices
    static_assert( EigenMatrix<Eigen::Matrix3d>, "Matrix3d should satisfy EigenMatrix" );
    static_assert( EigenMatrix<Eigen::MatrixXd>, "MatrixXd should satisfy EigenMatrix" );
    static_assert( EigenMatrix<Eigen::Matrix<double, 2, 3>>, "Fixed matrix should satisfy EigenMatrix" );

    // Static vectors
    static_assert( EigenVector<Eigen::Vector3d>, "Vector3d should satisfy EigenVector" );
    static_assert( EigenVector<Eigen::VectorXd>, "VectorXd should satisfy EigenVector" );
    static_assert( EigenVector<Eigen::RowVector3d>, "RowVector3d should satisfy EigenVector" );

    // Vectors are also matrices
    static_assert( EigenMatrix<Eigen::Vector3d>, "Vector3d should also satisfy EigenMatrix" );

    // But matrices with multiple rows and cols are not vectors
    static_assert( !EigenVector<Eigen::Matrix3d>, "Matrix3d should NOT satisfy EigenVector" );
    static_assert( !EigenVector<Eigen::MatrixXd>, "MatrixXd should NOT satisfy EigenVector" );

    // Test with references and const
    static_assert( EigenMatrix<Eigen::Matrix3d const&>, "const ref Matrix3d should satisfy EigenMatrix" );
    static_assert( EigenVector<Eigen::Vector3d&&>, "rvalue Vector3d should satisfy EigenVector" );

    BOOST_CHECK( true );
}

//
// Test UBlasMatrix and UBlasVector concepts
//
BOOST_AUTO_TEST_CASE( test_ublas_concepts )
{
    using matrix_t = ublas::matrix<double>;
    using vector_t = ublas::vector<double>;

    static_assert( UBlasMatrix<matrix_t>, "ublas::matrix should satisfy UBlasMatrix" );
    static_assert( UBlasVector<vector_t>, "ublas::vector should satisfy UBlasVector" );

    // Test with references
    static_assert( UBlasMatrix<matrix_t const&>, "const ref ublas::matrix should satisfy UBlasMatrix" );
    static_assert( UBlasVector<vector_t&>, "ref ublas::vector should satisfy UBlasVector" );

    // Cross-check: ublas types should not satisfy Eigen concepts
    static_assert( !EigenMatrix<matrix_t>, "ublas::matrix should NOT satisfy EigenMatrix" );
    static_assert( !EigenVector<vector_t>, "ublas::vector should NOT satisfy EigenVector" );

    BOOST_CHECK( true );
}

//
// Test MatrixLike and VectorLike concepts
//
BOOST_AUTO_TEST_CASE( test_matrix_vector_like_concepts )
{
    // MatrixLike accepts both Eigen and ublas
    static_assert( MatrixLike<Eigen::MatrixXd>, "Eigen::MatrixXd should satisfy MatrixLike" );
    static_assert( MatrixLike<ublas::matrix<double>>, "ublas::matrix should satisfy MatrixLike" );

    // VectorLike accepts both Eigen and ublas vectors
    static_assert( VectorLike<Eigen::VectorXd>, "Eigen::VectorXd should satisfy VectorLike" );
    static_assert( VectorLike<ublas::vector<double>>, "ublas::vector should satisfy VectorLike" );

    BOOST_CHECK( true );
}

//
// Test KokkosCompatibleStorage concept
//
BOOST_AUTO_TEST_CASE( test_kokkos_compatible_storage )
{
    // POD types should be Kokkos compatible
    static_assert( KokkosCompatibleStorage<double>, "double should be KokkosCompatibleStorage" );
    static_assert( KokkosCompatibleStorage<int>, "int should be KokkosCompatibleStorage" );

    // Arrays of POD should be compatible
    struct PodStruct { double x; double y; double z; };
    static_assert( KokkosCompatibleStorage<PodStruct>, "POD struct should be KokkosCompatibleStorage" );

    // std::vector is NOT Kokkos compatible (not trivially copyable)
    static_assert( !KokkosCompatibleStorage<std::vector<double>>, "std::vector should NOT be KokkosCompatibleStorage" );

    BOOST_CHECK( true );
}

//
// Test StaticOrder and DynamicOrder concepts
//
BOOST_AUTO_TEST_CASE( test_order_concepts )
{
    // Create order tags
    using static_order_3 = std::integral_constant<int, 3>;
    using static_order_0 = std::integral_constant<int, 0>;
    using dynamic_order = std::integral_constant<int, Dynamic>;

    static_assert( StaticOrder<static_order_3>, "Order 3 should be StaticOrder" );
    static_assert( StaticOrder<static_order_0>, "Order 0 should be StaticOrder" );
    static_assert( !StaticOrder<dynamic_order>, "Dynamic should NOT be StaticOrder" );

    static_assert( DynamicOrder<dynamic_order>, "Dynamic should be DynamicOrder" );
    static_assert( !DynamicOrder<static_order_3>, "Order 3 should NOT be DynamicOrder" );

    BOOST_CHECK( true );
}

//
// Test ScalarField, VectorField, Tensor2Field concepts
//
BOOST_AUTO_TEST_CASE( test_field_concepts )
{
    static_assert( ScalarFieldConcept<Scalar<2>>, "Scalar<2> should satisfy ScalarFieldConcept" );
    static_assert( ScalarFieldConcept<Scalar<3>>, "Scalar<3> should satisfy ScalarFieldConcept" );
    static_assert( !ScalarFieldConcept<Vectorial<3>>, "Vectorial<3> should NOT satisfy ScalarFieldConcept" );

    static_assert( VectorFieldConcept<Vectorial<2>>, "Vectorial<2> should satisfy VectorFieldConcept" );
    static_assert( VectorFieldConcept<Vectorial<3>>, "Vectorial<3> should satisfy VectorFieldConcept" );
    static_assert( !VectorFieldConcept<Scalar<3>>, "Scalar<3> should NOT satisfy VectorFieldConcept" );

    static_assert( Tensor2FieldConcept<Tensor2<2>>, "Tensor2<2> should satisfy Tensor2FieldConcept" );
    static_assert( Tensor2FieldConcept<Tensor2<3>>, "Tensor2<3> should satisfy Tensor2FieldConcept" );
    static_assert( Tensor2FieldConcept<Tensor2Symm<3>>, "Tensor2Symm<3> should satisfy Tensor2FieldConcept" );
    static_assert( !Tensor2FieldConcept<Vectorial<3>>, "Vectorial<3> should NOT satisfy Tensor2FieldConcept" );

    BOOST_CHECK( true );
}

//
// Test Convex concepts
//
BOOST_AUTO_TEST_CASE( test_convex_concepts )
{
    // Simplex convexes
    static_assert( ConvexConcept<Simplex<1>>, "Simplex<1> should satisfy ConvexConcept" );
    static_assert( ConvexConcept<Simplex<2>>, "Simplex<2> should satisfy ConvexConcept" );
    static_assert( ConvexConcept<Simplex<3>>, "Simplex<3> should satisfy ConvexConcept" );

    static_assert( SimplexConvex<Simplex<1>>, "Simplex<1> should satisfy SimplexConvex" );
    static_assert( SimplexConvex<Simplex<2>>, "Simplex<2> should satisfy SimplexConvex" );
    static_assert( SimplexConvex<Simplex<3>>, "Simplex<3> should satisfy SimplexConvex" );

    // Hypercube convexes
    static_assert( ConvexConcept<Hypercube<1>>, "Hypercube<1> should satisfy ConvexConcept" );
    static_assert( ConvexConcept<Hypercube<2>>, "Hypercube<2> should satisfy ConvexConcept" );
    static_assert( ConvexConcept<Hypercube<3>>, "Hypercube<3> should satisfy ConvexConcept" );

    static_assert( HypercubeConvex<Hypercube<1>>, "Hypercube<1> should satisfy HypercubeConvex" );
    static_assert( HypercubeConvex<Hypercube<2>>, "Hypercube<2> should satisfy HypercubeConvex" );
    static_assert( HypercubeConvex<Hypercube<3>>, "Hypercube<3> should satisfy HypercubeConvex" );

    // Cross-checks
    static_assert( !SimplexConvex<Hypercube<2>>, "Hypercube<2> should NOT satisfy SimplexConvex" );
    static_assert( !HypercubeConvex<Simplex<2>>, "Simplex<2> should NOT satisfy HypercubeConvex" );

    BOOST_CHECK( true );
}

//
// Test HasOrder concept
//
BOOST_AUTO_TEST_CASE( test_has_order_concept )
{
    // Types with static nOrder member
    static_assert( HasOrder<WithStaticOrder>, "Type with nOrder should satisfy HasOrder" );

    // Types with order() method
    static_assert( HasOrder<WithOrderMethod>, "Type with order() method should satisfy HasOrder" );

    // integral_constant types
    using order_tag = std::integral_constant<int, 4>;
    static_assert( HasOrder<order_tag>, "integral_constant should satisfy HasOrder" );

    BOOST_CHECK( true );
}

//
// Test PolynomialOrder concepts
//
BOOST_AUTO_TEST_CASE( test_polynomial_order_concepts )
{
    static_assert( PolynomialOrder<0>, "Order 0 should satisfy PolynomialOrder" );
    static_assert( PolynomialOrder<1>, "Order 1 should satisfy PolynomialOrder" );
    static_assert( PolynomialOrder<10>, "Order 10 should satisfy PolynomialOrder" );

    static_assert( LowOrder<0>, "Order 0 should satisfy LowOrder" );
    static_assert( LowOrder<1>, "Order 1 should satisfy LowOrder" );
    static_assert( LowOrder<2>, "Order 2 should satisfy LowOrder" );
    static_assert( !LowOrder<3>, "Order 3 should NOT satisfy LowOrder" );

    static_assert( !HighOrder<2>, "Order 2 should NOT satisfy HighOrder" );
    static_assert( HighOrder<3>, "Order 3 should satisfy HighOrder" );
    static_assert( HighOrder<10>, "Order 10 should satisfy HighOrder" );

    BOOST_CHECK( true );
}

//
// Runtime test to verify concepts work at runtime too
//
BOOST_AUTO_TEST_CASE( test_concepts_runtime_usage )
{
    // Function constrained by EigenMatrix concept
    auto process_eigen = []<EigenMatrix M>( M const& m ) {
        return m.rows() * m.cols();
    };

    Eigen::MatrixXd mat( 3, 4 );
    auto result = process_eigen( mat );
    BOOST_CHECK_EQUAL( result, 12 );

    // Function constrained by MatrixLike concept (works with both)
    auto get_size = []<MatrixLike M>( M const& m ) {
        if constexpr ( EigenMatrix<M> )
            return m.rows() * m.cols();
        else
            return m.size1() * m.size2();
    };

    ublas::matrix<double> ublas_mat( 2, 5 );
    BOOST_CHECK_EQUAL( get_size( mat ), 12 );
    BOOST_CHECK_EQUAL( get_size( ublas_mat ), 10 );
}

BOOST_AUTO_TEST_SUITE_END()
