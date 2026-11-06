/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
    SPDX-FileContributor: Hanna Chetouane

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * \file test_sb9pinching.cpp
 * \brief Unit tests for SB9 pinching kinematic operators.
 */

#define BOOST_TEST_MODULE test_sb9pinching
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/sb9_quadrature.hpp>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

namespace
{
/**
 * \brief Check that one field has zero SB9 pinching energy.
 *
 * The helper assembles the SB9 pinching bilinear form and evaluates it on 
 * the interpolated field. It is used for rigid translation and rigid rotation
 * fields, both of which must be in the kernel of the pinching strain operator.
 *
 * \tparam ProductSpacePtrType Feel++ displacement product function space pointer type.
 * \tparam TrialUType Trial basis proxy type for the Q1 vector displacement space.
 * \tparam TestUType Test basis proxy type for the Q1 vector displacement space.
 * \tparam TrialAType Trial basis proxy type for the scalar space carrying the 25th
 *         SB9 degree of freedom.
 * \tparam TestAType Test basis proxy type for the scalar space carrying the 25th
 *         SB9 degree of freedom.
 * \tparam FieldExprType Feel++ vector expression type for the checked field.
 * \param mesh Shell mesh.
 * \param Xh Product function space combining the Q1 vector displacement space 
 *        and the scalar space carrying the 25th SB9 degree of freedom.
 * \param u Trial basis proxy for the displacement space.
 * \param v Test basis proxy for the displacement space.
 * \param alpha Trial basis proxy for the scalar space.
 * \param beta Test basis proxy for the scalar space.
 * \param fieldExpr Field interpolated in the displacement component of \p Xh 
 *        before energy evaluation.
 */
template <typename ProductSpacePtrType, typename TrialUType, typename TestUType, typename TrialAType, typename TestAType, typename FieldExprType>
void
checkSb9PinchingRigidMode( mesh_ptrtype const& mesh,
                           ProductSpacePtrType const& Xh,
                           TrialUType const& u,
                           TestUType const& v,
                           TrialAType const& alpha,
                           TestAType const& beta,
                           FieldExprType const& fieldExpr )
{
    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;
    auto C = isotropic_stiffness<3>( lambda, mu );

    auto b = backend( _rebuild=true );
    auto pinching = blockform2( *Xh, solve::strategy::monolithic, b );
    auto zt = zeta();

    pinching( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=ddot( C, sb9Pinching( u, zt ), sb9Pinching( v, zt ) ) );
    pinching( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=ddot( C, sb9PinchingW9( alpha ), sb9Pinching( v, zt ) ) );
    pinching( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=ddot( C, sb9Pinching( u, zt ), sb9PinchingW9( beta ) ) );
    pinching( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=ddot( C, sb9PinchingW9( alpha ), sb9PinchingW9( beta ) ) );
    pinching.close();

    auto U = Xh->element();
    U( 0_c ).on( _range=elements( mesh ), _expr=fieldExpr, _close=true );

    U.buildVector( b );
    U.updateVectorFromSubVectors();
    BOOST_CHECK_SMALL( formEnergy( pinching, U.vectorMonolithic(), U.vectorMonolithic() ), g_tol );
}

/**
 * \brief Check rigid translations and rotations on one shell patch.
 *
 * \param mesh Shell mesh used to build the SB9 displacement product space.
 */
void
checkSb9PinchingRigidModesOnPatch( mesh_ptrtype const& mesh )
{
    auto Uh = Pchv<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto Xh = productPtr( Uh, Ah );

    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );
    auto alpha = trial( Ah, "alpha" );
    auto beta = test( Ah, "beta" );

    checkSb9PinchingRigidMode( mesh, Xh, u, v, alpha, beta, rigidTranslationField() );
    checkSb9PinchingRigidMode( mesh, Xh, u, v, alpha, beta, rigidRotationField() );
}

/**
 * \brief Compute SB9 pinching energy of a linear displacement field.
 *
 * This assembles the SB9 pinching bilinear form with `lambda = 0` and
 * `mu = 0.5`, then evaluates it on the displacement field `(ux, uy, uz)`.
 *
 * \tparam UxExpr Feel++ expression type for the x displacement component.
 * \tparam UyExpr Feel++ expression type for the y displacement component.
 * \tparam UzExpr Feel++ expression type for the z displacement component.
 * \param mesh Shell mesh used to build the SB9 displacement product space.
 * \param ux X displacement expression.
 * \param uy Y displacement expression.
 * \param uz Z displacement expression.
 * \return Quadratic pinching energy of the interpolated field.
 */
template <typename UxExpr, typename UyExpr, typename UzExpr>
double
sb9PinchingLinearEnergy( mesh_ptrtype const& mesh,
                         UxExpr const& ux,
                         UyExpr const& uy,
                         UzExpr const& uz )
{
    constexpr double lambda = 0.0;
    constexpr double mu = 0.5;
    auto C = isotropic_stiffness<3>( lambda, mu );

    auto Uh = Pchv<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto Xh = productPtr( Uh, Ah );

    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );
    auto alpha = trial( Ah, "alpha" );
    auto beta = test( Ah, "beta" );

    auto b = backend( _rebuild=true, _worldcomm=Uh->worldCommPtr() );
    auto a = blockform2( *Xh, solve::strategy::monolithic, b );
    auto zt = zeta();

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=ddot( C, sb9Pinching( u, zt ), sb9Pinching( v, zt ) ) );
    a( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=ddot( C, sb9PinchingW9( alpha ), sb9Pinching( v, zt ) ) );
    a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=ddot( C, sb9Pinching( u, zt ), sb9PinchingW9( beta ) ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=ddot( C, sb9PinchingW9( alpha ), sb9PinchingW9( beta ) ) );
    a.close();

    auto U = Xh->element();
    U( 0_c ).on( _range=elements( mesh ), _expr=vec(ux, uy, uz), _close=true );

    U.buildVector( b );
    U.updateVectorFromSubVectors();
    return formEnergy( a, U.vectorMonolithic(), U.vectorMonolithic() );
}

/**
 * \brief Create an oblique one-cell patch with a non-axis-aligned thickness.
 *
 * This geometry makes the `Bpz` pinching contribution observable in the
 * displacement/internal-mode coupling block.
 */
mesh_ptrtype
createObliqueThicknessPatch( std::string const& caseName )
{
    return createShellPatch( caseName,
                             test_vector_type( 0.0, 0.0, 0.0 ),
                             test_vector_type( 1.0, 0.0, 0.0 ),
                             test_vector_type( 0.0, 1.0, 0.0 ),
                             test_vector_type( 0.25, -0.15, 1.0 ) );
}

/**
 * \brief Return the norm of the displacement/internal pinching coupling block.
 *
 * With `bpzScale = 0`, the integrand is odd in `zeta` on the oblique patch and
 * the block must vanish. With `bpzScale = 1`, the `zeta * Bpz` term couples to
 * the internal `W9` mode and the block must be nonzero.
 */
double
sb9MixedCouplingBlockNorm( mesh_ptrtype const& mesh,
                           double bpzScale )
{
    auto Uh = Pchv<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto Xh = productPtr( Uh, Ah );
    auto b = backend( _rebuild=true, _worldcomm=Uh->worldCommPtr() );
    auto u = trial( Uh, "u" );
    auto beta = test( Ah, "beta" );

    auto a = blockform2( *Xh, solve::strategy::monolithic, b );
    auto C = isotropic_stiffness<3>( 0.0, 0.5 );
    auto sb9Quad = sb9ThroughThicknessLobatto5();

    a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                _quad=sb9Quad,
                                _expr=ddot( C,
                                            sb9Pinching( u, zeta(), cst( bpzScale ) ),
                                            sb9PinchingW9( beta ) ) );
    a.close();
    return a.l1Norm();
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_pinching_suite )

/**
 * \test Verify that SB9 pinching terms preserve rigid modes on a flat patch.
 */
BOOST_AUTO_TEST_CASE( sb9_pinching_preserves_rigid_modes_on_flat_patch )
{
    checkSb9PinchingRigidModesOnPatch( createFlatShellPatch( "sb9_pinching_flat_patch" ) );
}

/**
 * \test Verify that SB9 pinching terms preserve rigid modes after rotation.
 */
BOOST_AUTO_TEST_CASE( sb9_pinching_preserves_rigid_modes_on_rotated_patch )
{
    checkSb9PinchingRigidModesOnPatch( createRotatedShellPatch( "sb9_pinching_rotated_patch" ) );
}

/**
 * \test Verify simple linear fields on an axis-aligned patch.
 *
 * On the unit patch and with the chosen material constants, unit extension in
 * the z direction should produce energy one.
 * 
 * This test also verifies that the energy is quadratic and vanishes for a field
 * with zero pinching strain.
 */
BOOST_AUTO_TEST_CASE( sb9_pinching_linear_fields_have_expected_energy_on_axis_aligned_patch )
{
    auto mesh = createAxisAlignedUnitPatch( "sb9_pinching_axis_aligned_unit_patch" );

    double const energyZ = sb9PinchingLinearEnergy( mesh, cst( 0.0 ), cst( 0.0 ), Pz() );
    double const energyMZ = sb9PinchingLinearEnergy( mesh, cst( 0.0 ), cst( 0.0 ), -Pz() );
    double const energyX = sb9PinchingLinearEnergy( mesh, cst( 0.0 ), cst( 0.0 ), Px() );

    BOOST_CHECK_CLOSE( energyZ, 1.0, 1e-8 );
    BOOST_CHECK_SMALL( energyZ - energyMZ, g_tol );
    BOOST_CHECK_SMALL( energyX, g_tol );
}

/**
 * \test Verify that the `zeta * Bpz` pinching term couples to the internal
 * scalar pinching mode on an oblique-thickness patch.
 */
BOOST_AUTO_TEST_CASE( sb9_pinching_w9_blockform_couples_through_bpz )
{
    auto mesh = createObliqueThicknessPatch( "sb9_pinching_w9_oblique_patch" );

    double const bpcOnlyCoupling = sb9MixedCouplingBlockNorm( mesh, 0.0 );
    double const fullPinchingCoupling = sb9MixedCouplingBlockNorm( mesh, 1.0 );

    BOOST_TEST_MESSAGE( "sb9 pinching-w9 coupling norms: Bpc-only="
                        << bpcOnlyCoupling << " full=" << fullPinchingCoupling );
    BOOST_CHECK_SMALL( bpcOnlyCoupling, 1e-12 );
    BOOST_CHECK_GT( fullPinchingCoupling, 1e-8 );
}

BOOST_AUTO_TEST_SUITE_END()
