/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
 * @file test_nedelec_static.cpp
 * @brief Unisolvence, reproduction, support-envelope, and dimension tests for static Nedelec elements.
 */

#define BOOST_TEST_MODULE test_nedelec_static
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/nedelec.hpp>

#include <cmath>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/**
 * @brief Verify unisolvence and deterministic polynomial reproduction.
 * @tparam FE concrete static Nedelec finite element
 * @param label Boost.Test context label
 * @param expectedDof expected dimension of both primal and dual spaces
 */
template<typename FE>
void
checkLowestOrderNedelec( std::string const& label, uint16_type expectedDof )
{
    FE fe;

    BOOST_TEST_CONTEXT( label )
    {
        BOOST_REQUIRE_EQUAL( fe.primal().polynomialDimension(), expectedDof );
        BOOST_REQUIRE_EQUAL( FE::nLocalDof, expectedDof );

        auto const vandermonde = fe.dual()( fe.primal() );
        BOOST_REQUIRE_EQUAL( vandermonde.size1(), expectedDof );
        BOOST_REQUIRE_EQUAL( vandermonde.size2(), expectedDof );

        using eigen_matrix_type = Eigen::Matrix<typename FE::value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using eigen_vector_type = Eigen::Matrix<typename FE::value_type, Eigen::Dynamic, 1>;
        Eigen::Map<eigen_matrix_type const> vandermondeMap(
            vandermonde.data().begin(), vandermonde.size1(), vandermonde.size2() );
        Eigen::FullPivLU<eigen_matrix_type> lu( vandermondeMap );
        BOOST_CHECK_EQUAL( lu.rank(), expectedDof );

        IMGeneral<FE::nDim, typename FE::value_type, Simplex> im( 6 );
        auto const primalValues = fe.primal().evaluate( im.points() );
        auto const nodalValues = fe.evaluate( im.points() );
        Eigen::Map<eigen_matrix_type const> primalValuesMap(
            primalValues.data().begin(), primalValues.size1(), primalValues.size2() );
        Eigen::Map<eigen_matrix_type const> nodalValuesMap(
            nodalValues.data().begin(), nodalValues.size1(), nodalValues.size2() );

        eigen_vector_type primalCoefficients( expectedDof );
        for ( uint16_type j = 0; j < expectedDof; ++j )
            primalCoefficients( j ) = typename FE::value_type( j + 1 ) / typename FE::value_type( expectedDof + 1 );

        eigen_vector_type const dofValues = vandermondeMap * primalCoefficients;
        for ( int q = 0; q < im.nPoints(); ++q )
        {
            for ( uint16_type c = 0; c < FE::nComponents; ++c )
            {
                typename FE::value_type exact = 0;
                typename FE::value_type reconstructed = 0;
                for ( uint16_type j = 0; j < expectedDof; ++j )
                {
                    exact += primalCoefficients( j ) * primalValuesMap( FE::nComponents*j + c, q );
                    reconstructed += dofValues( j ) * nodalValuesMap( FE::nComponents*j + c, q );
                }
                BOOST_CHECK_SMALL( reconstructed - exact, typename FE::value_type( 1e-10 ) );
            }
        }
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( nedelec_static_suite )

/** @test Validate the supported lowest-order first-kind triangle and tetrahedron. */
BOOST_AUTO_TEST_CASE( lowest_order_ned1_is_unisolvent_and_reproducing )
{
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 0>::polynomialDegree == 1 );
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 0>::totalDof == 3 );
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 2>::edgeDof == 3 );
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 2>::cellDof == 6 );
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 2>::totalDof == 15 );
    static_assert( NedelecTraits<NedelecKind::NED1, 3, 2>::faceDof == 6 );
    static_assert( NedelecTraits<NedelecKind::NED1, 3, 2>::cellDof == 3 );
    static_assert( NedelecTraits<NedelecKind::NED1, 3, 2>::totalDof == 45 );
    static_assert( nedelecSimplexIsSupported<2, 0, NedelecKind::NED1> );
    static_assert( nedelecSimplexIsSupported<3, 0, NedelecKind::NED1> );
    static_assert( !nedelecSimplexIsSupported<3, 1, NedelecKind::NED1> );

    using ned1_2d = typename Nedelec<0, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using ned1_3d = typename Nedelec<0, NedelecKind::NED1>::template apply<3, 3, double, Simplex<3>>::type;

    checkLowestOrderNedelec<ned1_2d>( "NED1 triangle order 0", 3 );
    checkLowestOrderNedelec<ned1_3d>( "NED1 tetrahedron order 0", 6 );
}

/** @test Validate the supported lowest-order second-kind triangle. */
BOOST_AUTO_TEST_CASE( lowest_order_ned2_is_unisolvent_and_reproducing )
{
    static_assert( NedelecTraits<NedelecKind::NED2, 2, 0>::polynomialDegree == 1 );
    static_assert( NedelecTraits<NedelecKind::NED2, 2, 0>::edgeDof == 2 );
    static_assert( NedelecTraits<NedelecKind::NED2, 2, 0>::totalDof == 6 );
    static_assert( NedelecTraits<NedelecKind::NED2, 3, 2>::edgeDof == 4 );
    static_assert( NedelecTraits<NedelecKind::NED2, 3, 2>::faceDof == 8 );
    static_assert( NedelecTraits<NedelecKind::NED2, 3, 2>::cellDof == 4 );
    static_assert( NedelecTraits<NedelecKind::NED2, 3, 2>::totalDof == 60 );
    static_assert( nedelecSimplexIsSupported<2, 0, NedelecKind::NED2> );
    static_assert( !nedelecSimplexIsSupported<2, 1, NedelecKind::NED2> );
    static_assert( !nedelecSimplexIsSupported<3, 0, NedelecKind::NED2> );

    using ned2_2d = typename Nedelec<0, NedelecKind::NED2>::template apply<2, 2, double, Simplex<2>>::type;

    checkLowestOrderNedelec<ned2_2d>( "NED2 triangle order 0", 6 );
}

/** @test Validate extracted moment kernels for higher-order first-kind triangles. */
BOOST_AUTO_TEST_CASE( higher_order_ned1_triangle_uses_extracted_moment_kernels )
{
    using ned1_order1 = typename Nedelec<1, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using ned1_order2 = typename Nedelec<2, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( NedelecTraits<NedelecKind::NED1, 2, 1>::isSupported );
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 2>::isSupported );
    checkLowestOrderNedelec<ned1_order1>( "NED1 triangle order 1", 8 );
    checkLowestOrderNedelec<ned1_order2>( "NED1 triangle order 2", 15 );
}

BOOST_AUTO_TEST_SUITE_END()
