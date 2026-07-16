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
#include <feel/feelmesh/entities.hpp>

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

        std::array<uint16_type,4> attachedByDimension{};
        for ( uint16_type localDof = 0; localDof < FE::nLocalDof; ++localDof )
        {
            auto const attachment = fe.dofAttachment( localDof );
            BOOST_REQUIRE( attachment.isValid() );
            ++attachedByDimension[static_cast<uint16_type>( attachment.entityDim )];
        }
        BOOST_CHECK_EQUAL( attachedByDimension[0], FE::reference_convex_type::numVertices*FE::nDofPerVertex );
        BOOST_CHECK_EQUAL( attachedByDimension[1], FE::reference_convex_type::numEdges*FE::nDofPerEdge );
        BOOST_CHECK_EQUAL( attachedByDimension[2], FE::reference_convex_type::numFaces*FE::nDofPerFace );
        BOOST_CHECK_EQUAL( attachedByDimension[3], FE::nDofPerVolume );

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
                const auto tolerance = typename FE::value_type(
                    expectedDof >= 100 ? 5e-10 : 1e-10 );
                BOOST_CHECK_SMALL( reconstructed - exact, tolerance );
            }
        }
    }
}

/**
 * @brief Verify that the shared face interpolation kernel applies the same
 * reference functionals as the NED1 dual builder.
 *
 * Applying the face moments to the nodal basis must recover the face rows of
 * the identity interpolation matrix.  This isolates the sampling layout and
 * moment kernel from geometric mapping and global orientation handling.
 */
template<typename FE>
void
checkReferenceFaceInterpolationKernel()
{
    static_assert( FE::nDim == 3 );
    FE fe;
    auto const values = fe.evaluate( fe.points() );
    auto const& sampling = fe.dual().samplingLayout();

    struct ReferenceExpression
    {
        typename FE::matrix_type const& values;
        uint16_type basis = 0;

        [[nodiscard]] uint32_type nPoints() const noexcept
        {
            return static_cast<uint32_type>( values.size2() );
        }

        [[nodiscard]] typename FE::value_type
        evalq( uint16_type component, uint16_type, uint32_type point ) const
        {
            return values( FE::nComponents*basis+component, point );
        }
    };

    for ( uint16_type basis = 0; basis < FE::nLocalDof; ++basis )
    {
        ReferenceExpression expression{ values, basis };
        auto interpolant = fe.localInterpolant();
        if constexpr ( FE::kind == NedelecKind::NED1 )
        {
            for ( auto const& block : sampling.faces )
                Feel::detail::interpolateHCurlTriangularFaceMoments(
                    expression, interpolant, fe.referenceConvex(), block.face,
                    static_cast<int>( FE::nOrder )-2, block.firstDof,
                    block.elementPointOffset, block.momentCount );
        }
        else
        {
            RaviartThomasDynamicPolynomialSet<2> moments(
                RuntimeOrder{ static_cast<uint16_type>( FE::kind == NedelecKind::NED2
                                                            ? FE::nOrder-2 : 0 ) } );
            for ( auto const& block : sampling.faces )
                Feel::detail::interpolateHCurlTriangularFaceRTMoments(
                    expression, interpolant, fe.referenceConvex(), moments,
                    block.face, block.firstDof, block.elementPointOffset );

            if ( sampling.cellMomentCount > 0 )
            {
                RaviartThomasDynamicPolynomialSet<3> cellMoments(
                    RuntimeOrder{ static_cast<uint16_type>( FE::nOrder-3 ) } );
                Feel::detail::interpolateHCurlCellMoments<3>(
                    expression, interpolant, cellMoments,
                    static_cast<uint16_type>( FE::nOrder+1 ),
                    sampling.cellFirstDof, sampling.cellPointOffset,
                    sampling.cellMomentCount, FE::nComponents );
            }
        }

        for ( auto const& block : sampling.faces )
            for ( uint16_type ordinal = 0; ordinal < block.momentCount; ++ordinal )
            {
                const uint16_type dof = static_cast<uint16_type>( block.firstDof+ordinal );
                const double expected = dof == basis ? 1.0 : 0.0;
                BOOST_CHECK_SMALL( interpolant( dof )-expected, 2e-10 );
            }
        if constexpr ( FE::kind == NedelecKind::NED2 )
            for ( uint16_type ordinal = 0;
                  ordinal < sampling.cellMomentCount; ++ordinal )
            {
                const uint16_type dof = static_cast<uint16_type>(
                    sampling.cellFirstDof+ordinal );
                const double expected = dof == basis ? 1.0 : 0.0;
                BOOST_CHECK_SMALL( interpolant( dof )-expected, 2e-10 );
            }
    }
}

/**
 * @brief Verify that the 3D first-kind primal contains every vector affine
 * polynomial.
 *
 * The public orders tested here must contain `[P_1]^3`.  Checking this
 * directly complements unisolvence: a square primal/dual pairing can be
 * nonsingular while spanning the wrong polynomial subspace.
 */
template<typename FE>
void
checkVectorAffineContainment()
{
    static_assert( FE::nDim == 3 );
    FE fe;
    IMGeneral<3, typename FE::value_type, Simplex> quadrature( 8 );
    auto const values = fe.evaluate( quadrature.points() );
    const Eigen::Index sampleCount = 3*quadrature.nPoints();
    Eigen::MatrixXd basisValues( sampleCount, FE::nLocalDof );
    for ( int q = 0; q < quadrature.nPoints(); ++q )
        for ( int component = 0; component < 3; ++component )
            for ( uint16_type basis = 0; basis < FE::nLocalDof; ++basis )
                basisValues( 3*q+component, basis ) = values( 3*basis+component, q );

    for ( int targetComponent = 0; targetComponent < 3; ++targetComponent )
        for ( int monomial = 0; monomial < 4; ++monomial )
        {
            Eigen::VectorXd target = Eigen::VectorXd::Zero( sampleCount );
            for ( int q = 0; q < quadrature.nPoints(); ++q )
                target( 3*q+targetComponent ) = monomial == 0
                    ? 1.0 : quadrature.points()( monomial-1, q );
            Eigen::VectorXd const coefficients =
                basisValues.colPivHouseholderQr().solve( target );
            BOOST_CHECK_SMALL( ( basisValues*coefficients-target ).norm(), 2e-10 );
        }
}

/** Check a guarded NED2 primal/dual pair before exposing its FE shell. */
template<int Dim, int Order>
void
checkCandidateNed2CiarletPair( uint16_type expectedDof )
{
    using primal_type = Feel::detail::OrthonormalPolynomialSet<
        Dim, Order+1, Dim, Vectorial, double, 0, Simplex>;
    using dual_type = fem::detail::Nedelec2Dual<primal_type, PointSetEquiSpaced>;
    primal_type primal;
    dual_type dual( primal );
    auto const matrix = dual( primal );
    BOOST_REQUIRE_EQUAL( primal.polynomialDimension(), expectedDof );
    BOOST_REQUIRE_EQUAL( matrix.size1(), expectedDof );
    BOOST_REQUIRE_EQUAL( matrix.size2(), expectedDof );
    Eigen::MatrixXd dense( expectedDof, expectedDof );
    for ( uint16_type row = 0; row < expectedDof; ++row )
        for ( uint16_type column = 0; column < expectedDof; ++column )
            dense( row, column ) = matrix( row, column );
    Eigen::FullPivLU<Eigen::MatrixXd> rank( dense );
    rank.setThreshold( 1e-11 );
    BOOST_CHECK_EQUAL( rank.rank(), expectedDof );
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
    static_assert( nedelecSimplexIsSupported<3, 1, NedelecKind::NED1> );
    static_assert( nedelecSimplexIsSupported<3, 2, NedelecKind::NED1> );
    static_assert( nedelecSimplexIsSupported<3, 3, NedelecKind::NED1> );

    using ned1_2d = typename Nedelec<0, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using ned1_3d = typename Nedelec<0, NedelecKind::NED1>::template apply<3, 3, double, Simplex<3>>::type;

    checkLowestOrderNedelec<ned1_2d>( "NED1 triangle order 0", 3 );
    checkLowestOrderNedelec<ned1_3d>( "NED1 tetrahedron order 0", 6 );
}

/** @test Validate second-kind FE shells on triangles/tetrahedra through order three. */
BOOST_AUTO_TEST_CASE( ned2_is_unisolvent_and_reproducing_through_order_three )
{
    static_assert( NedelecTraits<NedelecKind::NED2, 2, 0>::polynomialDegree == 1 );
    static_assert( NedelecTraits<NedelecKind::NED2, 2, 0>::edgeDof == 2 );
    static_assert( NedelecTraits<NedelecKind::NED2, 2, 0>::totalDof == 6 );
    static_assert( NedelecTraits<NedelecKind::NED2, 3, 2>::edgeDof == 4 );
    static_assert( NedelecTraits<NedelecKind::NED2, 3, 2>::faceDof == 8 );
    static_assert( NedelecTraits<NedelecKind::NED2, 3, 2>::cellDof == 4 );
    static_assert( NedelecTraits<NedelecKind::NED2, 3, 2>::totalDof == 60 );
    static_assert( nedelecSimplexIsSupported<2, 0, NedelecKind::NED2> );
    static_assert( nedelecSimplexIsSupported<2, 1, NedelecKind::NED2> );
    static_assert( nedelecSimplexIsSupported<3, 0, NedelecKind::NED2> );
    static_assert( nedelecSimplexIsSupported<3, 3, NedelecKind::NED2> );

    using ned2_2d_0 = typename Nedelec<0, NedelecKind::NED2>::template apply<2, 2, double, Simplex<2>>::type;
    using ned2_2d_1 = typename Nedelec<1, NedelecKind::NED2>::template apply<2, 2, double, Simplex<2>>::type;
    using ned2_2d_2 = typename Nedelec<2, NedelecKind::NED2>::template apply<2, 2, double, Simplex<2>>::type;
    using ned2_2d_3 = typename Nedelec<3, NedelecKind::NED2>::template apply<2, 2, double, Simplex<2>>::type;
    using ned2_3d_0 = typename Nedelec<0, NedelecKind::NED2>::template apply<3, 3, double, Simplex<3>>::type;
    using ned2_3d_1 = typename Nedelec<1, NedelecKind::NED2>::template apply<3, 3, double, Simplex<3>>::type;
    using ned2_3d_2 = typename Nedelec<2, NedelecKind::NED2>::template apply<3, 3, double, Simplex<3>>::type;
    using ned2_3d_3 = typename Nedelec<3, NedelecKind::NED2>::template apply<3, 3, double, Simplex<3>>::type;

    checkLowestOrderNedelec<ned2_2d_0>( "NED2 triangle order 0", 6 );
    checkLowestOrderNedelec<ned2_2d_1>( "NED2 triangle order 1", 12 );
    checkLowestOrderNedelec<ned2_2d_2>( "NED2 triangle order 2", 20 );
    checkLowestOrderNedelec<ned2_2d_3>( "NED2 triangle order 3", 30 );
    checkLowestOrderNedelec<ned2_3d_0>( "NED2 tetrahedron order 0", 12 );
    checkLowestOrderNedelec<ned2_3d_1>( "NED2 tetrahedron order 1", 30 );
    checkLowestOrderNedelec<ned2_3d_2>( "NED2 tetrahedron order 2", 60 );
    checkLowestOrderNedelec<ned2_3d_3>( "NED2 tetrahedron order 3", 105 );
}

/**
 * @test Audit candidate higher-order NED2 moments before lifting public guards.
 */
BOOST_AUTO_TEST_CASE( candidate_higher_order_ned2_moments_are_unisolvent )
{
    checkCandidateNed2CiarletPair<2, 1>( 12 );
    checkCandidateNed2CiarletPair<2, 2>( 20 );
    checkCandidateNed2CiarletPair<2, 3>( 30 );
    checkCandidateNed2CiarletPair<3, 0>( 12 );
    checkCandidateNed2CiarletPair<3, 1>( 30 );
    checkCandidateNed2CiarletPair<3, 2>( 60 );
    checkCandidateNed2CiarletPair<3, 3>( 105 );
}

/** @test Validate extracted moment kernels for higher-order first-kind triangles. */
BOOST_AUTO_TEST_CASE( higher_order_ned1_triangle_uses_extracted_moment_kernels )
{
    using ned1_order1 = typename Nedelec<1, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using ned1_order2 = typename Nedelec<2, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using ned1_order3 = typename Nedelec<3, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using ned1_order4 = typename Nedelec<4, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( NedelecTraits<NedelecKind::NED1, 2, 1>::isSupported );
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 2>::isSupported );
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 3>::isSupported );
    static_assert( NedelecTraits<NedelecKind::NED1, 2, 4>::isSupported );
    checkLowestOrderNedelec<ned1_order1>( "NED1 triangle order 1", 8 );
    checkLowestOrderNedelec<ned1_order2>( "NED1 triangle order 2", 15 );
    checkLowestOrderNedelec<ned1_order3>( "NED1 triangle order 3", 24 );
    checkLowestOrderNedelec<ned1_order4>( "NED1 triangle order 4", 35 );
}

/** @test Validate complete static first-kind tetrahedra through order three. */
BOOST_AUTO_TEST_CASE( higher_order_ned1_tetrahedron_is_unisolvent_and_reproducing )
{
    using ned1_order1 = typename Nedelec<1, NedelecKind::NED1>::template apply<3, 3, double, Simplex<3>>::type;
    using ned1_order2 = typename Nedelec<2, NedelecKind::NED1>::template apply<3, 3, double, Simplex<3>>::type;
    using ned1_order3 = typename Nedelec<3, NedelecKind::NED1>::template apply<3, 3, double, Simplex<3>>::type;

    checkLowestOrderNedelec<ned1_order1>( "NED1 tetrahedron order 1", 20 );
    checkLowestOrderNedelec<ned1_order2>( "NED1 tetrahedron order 2", 45 );
    checkLowestOrderNedelec<ned1_order3>( "NED1 tetrahedron order 3", 84 );
    checkVectorAffineContainment<ned1_order1>();
    checkVectorAffineContainment<ned1_order2>();
    checkVectorAffineContainment<ned1_order3>();
}

/** @test Match the shared tetrahedral face interpolant to the static dual. */
BOOST_AUTO_TEST_CASE( ned1_tetrahedron_face_interpolation_matches_dual )
{
    using ned1_order1 = typename Nedelec<1, NedelecKind::NED1>::template apply<3, 3, double, Simplex<3>>::type;
    using ned1_order2 = typename Nedelec<2, NedelecKind::NED1>::template apply<3, 3, double, Simplex<3>>::type;
    using ned1_order3 = typename Nedelec<3, NedelecKind::NED1>::template apply<3, 3, double, Simplex<3>>::type;

    checkReferenceFaceInterpolationKernel<ned1_order1>();
    checkReferenceFaceInterpolationKernel<ned1_order2>();
    checkReferenceFaceInterpolationKernel<ned1_order3>();
}

/** @test Match NED2 tetrahedral RT face interpolation to the static dual. */
BOOST_AUTO_TEST_CASE( ned2_tetrahedron_face_interpolation_matches_dual )
{
    using ned2_order1 = typename Nedelec<1, NedelecKind::NED2>::template apply<3, 3, double, Simplex<3>>::type;
    using ned2_order2 = typename Nedelec<2, NedelecKind::NED2>::template apply<3, 3, double, Simplex<3>>::type;
    using ned2_order3 = typename Nedelec<3, NedelecKind::NED2>::template apply<3, 3, double, Simplex<3>>::type;

    checkReferenceFaceInterpolationKernel<ned2_order1>();
    checkReferenceFaceInterpolationKernel<ned2_order2>();
    checkReferenceFaceInterpolationKernel<ned2_order3>();
}

/** @test Validate all six triangular-face actions and their group law. */
BOOST_AUTO_TEST_CASE( ned1_triangular_face_orientation_is_a_group_representation )
{
    using permutation_type = triangular_faces_type;
    using matrix_type = Eigen::MatrixXd;
    std::array<std::array<uint16_type,3>,7> const vertices{{
        {{0,1,2}}, {{0,1,2}}, {{2,1,0}}, {{0,2,1}},
        {{1,0,2}}, {{1,2,0}}, {{2,0,1}}
    }};
    auto permutationFromVertices = [&vertices]( std::array<uint16_type,3> const& candidate )
    {
        for ( uint16_type p = permutation_type::IDENTITY; p < permutation_type::N_PERMUTATIONS; ++p )
            if ( vertices[p] == candidate )
                return p;
        return uint16_type( 0 );
    };
    auto dense = []( auto const& storage )
    {
        matrix_type result( storage->rows, storage->cols );
        for ( uint16_type i = 0; i < storage->rows; ++i )
            for ( uint16_type j = 0; j < storage->cols; ++j )
                result( i, j ) = ( *storage )( i, j );
        return result;
    };

    for ( uint16_type degree = 0; degree <= 2; ++degree )
    {
        std::array<matrix_type,7> transforms;
        for ( uint16_type p = permutation_type::IDENTITY; p < permutation_type::N_PERMUTATIONS; ++p )
        {
            transforms[p] = dense( Feel::detail::nedelecTriangularFaceMomentTransform(
                permutation_type( p ), degree ) );
            BOOST_CHECK_SMALL( ( transforms[p].fullPivLu().rank() - transforms[p].rows() ) * 1., 1e-12 );
        }
        BOOST_CHECK_SMALL( ( transforms[permutation_type::IDENTITY] -
                             matrix_type::Identity( transforms[1].rows(), transforms[1].cols() ) ).norm(), 1e-11 );

        for ( uint16_type p = permutation_type::IDENTITY; p < permutation_type::N_PERMUTATIONS; ++p )
            for ( uint16_type q = permutation_type::IDENTITY; q < permutation_type::N_PERMUTATIONS; ++q )
            {
                std::array<uint16_type,3> composition{};
                for ( uint16_type i = 0; i < 3; ++i )
                    composition[i] = vertices[p][vertices[q][i]];
                auto const pq = permutationFromVertices( composition );
                BOOST_REQUIRE_NE( pq, 0 );
                BOOST_CHECK_SMALL( ( transforms[q] * transforms[p] - transforms[pq] ).norm(), 2e-10 );
            }
    }
}

/** @test Validate all six NED2 triangular-face RT moment actions. */
BOOST_AUTO_TEST_CASE( ned2_triangular_face_orientation_is_a_group_representation )
{
    using matrix_type = Eigen::MatrixXd;
    std::array<std::array<uint16_type,3>,6> const vertices{{
        {{0,1,2}}, {{2,1,0}}, {{0,2,1}},
        {{1,0,2}}, {{1,2,0}}, {{2,0,1}}
    }};
    auto permutationFromVertices = [&vertices]( std::array<uint16_type,3> const& candidate )
    {
        for ( uint16_type p = 0; p < vertices.size(); ++p )
            if ( vertices[p] == candidate )
                return p;
        return static_cast<uint16_type>( vertices.size() );
    };
    auto dense = []( auto const& storage )
    {
        matrix_type result( storage->rows, storage->cols );
        for ( uint16_type i = 0; i < storage->rows; ++i )
            for ( uint16_type j = 0; j < storage->cols; ++j )
                result( i, j ) = ( *storage )( i, j );
        return result;
    };

    for ( uint16_type publicOrder = 1; publicOrder <= 3; ++publicOrder )
    {
        std::array<matrix_type,6> transforms;
        for ( uint16_type p = 0; p < vertices.size(); ++p )
        {
            transforms[p] = dense(
                Feel::detail::nedelecTriangularFaceRTMomentTransform(
                    vertices[p], publicOrder ) );
            BOOST_CHECK_EQUAL( transforms[p].fullPivLu().rank(),
                               transforms[p].rows() );
            auto const singularValues = transforms[p].jacobiSvd().singularValues();
            BOOST_CHECK_LT( singularValues( 0 )/
                            singularValues( singularValues.size()-1 ), 3. );

            RaviartThomasDynamicPolynomialSet<2> moments(
                RuntimeOrder{ static_cast<uint16_type>( publicOrder-1 ) } );
            IMGeneral<2, double, Simplex> quadrature(
                static_cast<uint16_type>( 2*( publicOrder+3 ) ) );
            auto const momentValues = moments.evaluate( quadrature.points() );
            Eigen::VectorXd canonicalMoments = Eigen::VectorXd::Zero(
                moments.polynomialDimension() );
            Eigen::VectorXd localMoments = Eigen::VectorXd::Zero(
                moments.polynomialDimension() );
            std::array<Eigen::Vector2d,3> const referenceVertices{{
                Eigen::Vector2d( -1., -1. ), Eigen::Vector2d( 1., -1. ),
                Eigen::Vector2d( -1., 1. ) }};
            Eigen::Vector2d const referenceOrigin = referenceVertices[0];
            Eigen::Vector2d const origin = referenceVertices[vertices[p][0]];
            Eigen::Matrix2d jacobian;
            jacobian.col( 0 ) =
                ( referenceVertices[vertices[p][1]]-origin )/2.;
            jacobian.col( 1 ) =
                ( referenceVertices[vertices[p][2]]-origin )/2.;
            for ( uint16_type q = 0; q < quadrature.nPoints(); ++q )
            {
                Eigen::Vector2d const localPoint(
                    quadrature.points()( 0, q ), quadrature.points()( 1, q ) );
                Eigen::Vector2d const canonicalPoint =
                    origin+jacobian*( localPoint-referenceOrigin );
                auto field = []( Eigen::Vector2d const& x )
                {
                    return Eigen::Vector2d(
                        1.+x[0]+x[0]*x[1], 2.-x[1]+x[0]*x[0] );
                };
                Eigen::Vector2d const canonicalValue = field( localPoint );
                Eigen::Vector2d const localValue =
                    jacobian.transpose()*field( canonicalPoint );
                for ( uint16_type moment = 0;
                      moment < moments.polynomialDimension(); ++moment )
                {
                    Eigen::Vector2d const testValue(
                        momentValues( 2*moment, q ),
                        momentValues( 2*moment+1, q ) );
                    canonicalMoments( moment ) += quadrature.weight( q )*
                        canonicalValue.dot( testValue );
                    localMoments( moment ) += quadrature.weight( q )*
                        localValue.dot( testValue );
                }
            }
            BOOST_CHECK_SMALL(
                ( localMoments-transforms[p]*canonicalMoments ).norm(),
                2e-10 );
        }
        BOOST_CHECK_SMALL(
            ( transforms[0] - matrix_type::Identity(
                  transforms[0].rows(), transforms[0].cols() ) ).norm(),
            1e-11 );

        for ( uint16_type p = 0; p < vertices.size(); ++p )
            for ( uint16_type q = 0; q < vertices.size(); ++q )
            {
                std::array<uint16_type,3> composition{};
                for ( uint16_type i = 0; i < 3; ++i )
                    composition[i] = vertices[p][vertices[q][i]];
                auto const pq = permutationFromVertices( composition );
                BOOST_REQUIRE_LT( pq, vertices.size() );
                BOOST_CHECK_SMALL(
                    ( transforms[q]*transforms[p] - transforms[pq] ).norm(),
                    2e-10 );
            }
    }
}

BOOST_AUTO_TEST_SUITE_END()
