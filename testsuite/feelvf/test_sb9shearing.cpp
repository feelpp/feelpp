/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
 * \file test_sb9shearing.cpp
 * \brief Unit tests for SB9 shearing kinematic operators.
 */

#define BOOST_TEST_MODULE test_sb9shearing
#include <feel/feelcore/testsuite.hpp>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

namespace
{
/**
 * \brief Check that one field has zero SB9 transverse-shearing energy.
 *
 * The helper assembles the SB9 shearing bilinear form and
 * evaluates it on the interpolated field. It is used for rigid translation and
 * rigid rotation fields, both of which must be in the kernel of the shearing
 * strain operator.
 *
 * \tparam SpaceType Feel++ displacement function space type.
 * \tparam TrialType Trial basis proxy type.
 * \tparam TestType Test basis proxy type.
 * \tparam FieldExprType Feel++ vector expression type for the checked field.
 * \param Uh Displacement function space.
 * \param u Trial basis proxy.
 * \param v Test basis proxy.
 * \param fieldExpr Field interpolated in \p Uh before energy evaluation.
 */
template <typename SpaceType, typename TrialType, typename TestType, typename FieldExprType>
void
checkSb9ShearingRigidMode( SpaceType const& Uh,
                           TrialType const& u,
                           TestType const& v,
                           FieldExprType const& fieldExpr )
{
    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;
    auto C = isotropic_stiffness<3>( lambda, mu );

    auto shearing = form2( _trial=Uh, _test=Uh );
    shearing = integrate( _range=elements( Uh->mesh() ),
                          _expr=ddot( C, sb9Shearing( u, zeta() ), sb9Shearing( v, zeta() ) ) );

    auto uh = Uh->element( "uh" );
    auto vh = Uh->element( "vh" );
    uh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );
    vh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );

    BOOST_CHECK_SMALL( formEnergy( shearing, vh, uh ), g_tol );
}

/**
 * \brief Check rigid translations and rotations on one shell patch.
 *
 * \param mesh Shell mesh used to build the Q1 vector displacement space.
 */
void
checkSb9ShearingRigidModesOnPatch( mesh_ptrtype const& mesh )
{
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );

    checkSb9ShearingRigidMode( Uh, u, v, rigidTranslationField() );
    checkSb9ShearingRigidMode( Uh, u, v, rigidRotationField() );
}

/**
 * \brief Compute SB9 shearing energy of a linear displacement field.
 *
 * This assembles the SB9 shearing bilinear form with
 * `lambda = 0` and `mu = 0.5`, then evaluates it on the displacement field
 * `(ux, uy, uz)`.
 *
 * \tparam UxExpr Feel++ expression type for the x displacement component.
 * \tparam UyExpr Feel++ expression type for the y displacement component.
 * \tparam UzExpr Feel++ expression type for the z displacement component.
 * \param mesh Shell mesh used to build the Q1 vector displacement space.
 * \param ux X displacement expression.
 * \param uy Y displacement expression.
 * \param uz Z displacement expression.
 * \return Quadratic shearing energy of the interpolated field.
 */
template <typename UxExpr, typename UyExpr, typename UzExpr>
double
sb9ShearingLinearEnergy( mesh_ptrtype const& mesh,
                         UxExpr const& ux,
                         UyExpr const& uy,
                         UzExpr const& uz )
{
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );
    constexpr double lambda = 0.0;
    constexpr double mu = 0.5;
    auto C = isotropic_stiffness<3>( lambda, mu );

    auto a = form2( _trial=Uh, _test=Uh );
    a = integrate( _range=elements( mesh ),
                   _expr=ddot( C, sb9Shearing( u, zeta() ), sb9Shearing( v, zeta() ) ) );

    auto uh = Uh->element( "uh" );
    uh.on( _range=elements( mesh ), _expr=vec( ux, uy, uz ), _close=true );  
    return formEnergy( a, uh, uh );
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_shearing_suite )

/**
 * \test Verify that SB9 shearing terms preserve rigid modes on a flat patch.
 */
BOOST_AUTO_TEST_CASE( sb9_shearing_preserves_rigid_modes_on_flat_patch )
{
    checkSb9ShearingRigidModesOnPatch( createFlatShellPatch( "sb9_shearing_flat_patch" ) );
}

/**
 * \test Verify that SB9 shearing terms preserve rigid modes after rotation.
 */
BOOST_AUTO_TEST_CASE( sb9_shearing_preserves_rigid_modes_on_rotated_patch )
{
    checkSb9ShearingRigidModesOnPatch( createRotatedShellPatch( "sb9_shearing_rotated_patch" ) );
}

/**
 * \test Verify SB9's shearing term preserve properties (symetry and quadratic homogeneity).
 */
BOOST_AUTO_TEST_CASE( sb9_shearing_properties_on_axis_aligned_patch )  
{
    auto mesh = createAxisAlignedUnitPatch( "sb9_shearing_axis_aligned_unit_patch" );

    double const energyX = sb9ShearingLinearEnergy( mesh, cst( 0.0 ), cst( 0.0 ), Px() );
    double const energyY = sb9ShearingLinearEnergy( mesh, cst( 0.0 ), cst( 0.0 ), Py() );
    double const energy2X = sb9ShearingLinearEnergy( mesh, cst( 0.0 ), cst( 0.0 ), cst( 2.0 ) * Px() );

    BOOST_CHECK_SMALL( energyX - energyY, g_tol );
    BOOST_CHECK_SMALL( 4*energyX - energy2X, g_tol );
}

/**
 * \test Verify simple axis-aligned in-plane extension energies.
 *
 * On the unit patch and with the chosen material constants, unit transverse shear
 * deformation would have unit energy. Due to the shear correction factor embedded 
 * in the SB9 transverse-shearing formulation, the expected reference energy is 5/6.
 *
 * This test also verifies that in-plane membrane modes do not contribute to the SB9 
 * transverse-shearing energy.
 */
BOOST_AUTO_TEST_CASE( sb9_shearing_linear_have_expected_energy_on_axis_aligned_patch )
{
    auto mesh = createAxisAlignedUnitPatch( "sb9_shearing_axis_aligned_unit_patch" );

    double const energyX = sb9ShearingLinearEnergy( mesh, cst( 0.0 ), cst( 0.0 ), Px() );
    double const energyMem = sb9ShearingLinearEnergy( mesh, Px(), Py(), cst( 0.0 ) );

    BOOST_CHECK_CLOSE( energyX, 5.0/6.0, g_tol );
    BOOST_CHECK_SMALL( energyMem, g_tol );
}

BOOST_AUTO_TEST_SUITE_END()
