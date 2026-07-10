/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
    SPDX-FileContributor: Hanna Chetouane

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * \file test_sb9stabilization.cpp
 * \brief Unit tests for SB9 stabilization kinematic operators.
 */

#define BOOST_TEST_MODULE test_sb9stabilization
#include <feel/feelcore/testsuite.hpp>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

namespace
{
/**
 * \brief Check that one field has zero SB9 stabilization energy.
 *
 * These helpers assemble the different SB9 stabilization bilinear forms (mode
 * stabilization, transverse shear stabilization, and the sum of both) and
 * evaluates its energy on a given interpolated field. It is used for rigid
 * translation and rigid rotation fields, both of which must be in the kernel
 * of the stabilization strain operator.
 *
 * \tparam FieldExprType Feel++ vector expression type for the checked field.
 * \param mesh Shell mesh used to build the Q1 vector displacement space.
 * \param fieldExpr Field interpolated in \p Uh before energy evaluation.
 */
template <typename FieldExprType>
void
checkSb9ModeStabilizationRigidMode( mesh_ptrtype const& mesh, FieldExprType const& fieldExpr )
{
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );

    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;

    auto modeStabilization = form2( _trial=Uh, _test=Uh );
    modeStabilization = integrate( _range=elements( mesh ),
                                    _quad=sb9ThroughThicknessLobatto5(),
                                    _expr=sb9ModeStabilization(u, v, lambda, mu ) );

    auto uh = Uh->element( "uh" );
    auto vh = Uh->element( "vh" );
    uh.on( _range=elements( mesh ), _expr=fieldExpr, _close=true );
    vh.on( _range=elements( mesh ), _expr=fieldExpr, _close=true );

    BOOST_CHECK_SMALL( formEnergy( modeStabilization, vh, uh ), g_tol );
}
template <typename FieldExprType>
void
checkSb9ShearingStabilizationRigidMode( mesh_ptrtype const& mesh, FieldExprType const& fieldExpr )
{
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );

    constexpr double mu = 1.1;

    auto shearStabilization = form2( _trial=Uh, _test=Uh );
    shearStabilization = integrate( _range=elements( mesh ),
                                    _quad=sb9ThroughThicknessLobatto5(),
                                    _expr=sb9ShearingStabilization(u, v, mu ) );

    auto uh = Uh->element( "uh" );
    auto vh = Uh->element( "vh" );
    uh.on( _range=elements( mesh ), _expr=fieldExpr, _close=true );
    vh.on( _range=elements( mesh ), _expr=fieldExpr, _close=true );

    BOOST_CHECK_SMALL( formEnergy( shearStabilization, vh, uh ), g_tol );
}
template <typename FieldExprType>
void
checkSb9StabilizationRigidMode( mesh_ptrtype const& mesh, FieldExprType const& fieldExpr )
{
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );

    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;

    auto stabilization = form2( _trial=Uh, _test=Uh );
    stabilization = integrate( _range=elements( mesh ),
                               _quad=sb9ThroughThicknessLobatto5(),
                               _expr=sb9Stabilization(u, v, lambda, mu ) );

    auto uh = Uh->element( "uh" );
    auto vh = Uh->element( "vh" );
    uh.on( _range=elements( mesh ), _expr=fieldExpr, _close=true );
    vh.on( _range=elements( mesh ), _expr=fieldExpr, _close=true );

    BOOST_CHECK_SMALL( formEnergy( stabilization, vh, uh ), g_tol );
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_stabilization_suite )

/**
 * \test Verify that SB9 stabilization terms preserve rigid modes on a flat patch.
 */
BOOST_AUTO_TEST_CASE( sb9_mode_stabilization_preserves_rigid_modes_on_flat_patch )
{
    auto mesh = createFlatShellPatch( "sb9_mode_stabilization_flat_patch" );
    checkSb9ModeStabilizationRigidMode( mesh, rigidTranslationField() );
    checkSb9ModeStabilizationRigidMode( mesh, rigidRotationField() );
}
BOOST_AUTO_TEST_CASE( sb9_shearing_stabilization_preserves_rigid_modes_on_flat_patch )
{
    auto mesh = createFlatShellPatch( "sb9_shearing_stabilization_flat_patch" );
    checkSb9ShearingStabilizationRigidMode( mesh, rigidTranslationField() );
    checkSb9ShearingStabilizationRigidMode( mesh, rigidRotationField() );
}
BOOST_AUTO_TEST_CASE( sb9_stabilization_preserves_rigid_modes_on_flat_patch )
{
    auto mesh = createFlatShellPatch( "sb9_stabilization_flat_patch" );
    checkSb9StabilizationRigidMode( mesh, rigidTranslationField() );
    checkSb9StabilizationRigidMode( mesh, rigidRotationField() );
}

/**
 * \test Verify that SB9 stabilization terms preserve rigid modes after rotation.
 */
BOOST_AUTO_TEST_CASE( sb9_mode_stabilization_preserves_rigid_modes_on_rotated_patch )
{
    auto mesh = createRotatedShellPatch( "sb9_mode_stabilization_rotated_patch" );
    checkSb9ModeStabilizationRigidMode( mesh, rigidTranslationField() );
    checkSb9ModeStabilizationRigidMode( mesh, rigidRotationField() );
}
BOOST_AUTO_TEST_CASE( sb9_shearing_stabilization_preserves_rigid_modes_on_rotated_patch )
{
    auto mesh = createRotatedShellPatch( "sb9_shearing_stabilization_rotated_patch" );
    checkSb9ShearingStabilizationRigidMode( mesh, rigidTranslationField() );
    checkSb9ShearingStabilizationRigidMode( mesh, rigidRotationField() );
}
BOOST_AUTO_TEST_CASE( sb9_stabilization_preserves_rigid_modes_on_rotated_patch )
{
    auto mesh = createRotatedShellPatch( "sb9_stabilization_rotated_patch" );
    checkSb9StabilizationRigidMode( mesh, rigidTranslationField() );
    checkSb9StabilizationRigidMode( mesh, rigidRotationField() );
}

BOOST_AUTO_TEST_SUITE_END()
