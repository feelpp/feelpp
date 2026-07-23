/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
    SPDX-FileContributor: Hanna Chetouane

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
 * \file test_sb9stabilization.cpp
 * \brief Unit tests for SB9 mode-stabilization operators.
 *
 * The tests verify the SB9 mode-stabilization row layout:
 * `Bs1` and `Bs2` are one-row blocks, `Bs3` is a two-row block, and `Bs4` is a
 * three-row block. Each block must preserve rigid translations and rigid
 * rotations on flat and rotated shell patches.
 */

#define BOOST_TEST_MODULE test_sb9stabilization
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelvf/sb9_stabilization.hpp>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

namespace
{
/**
 * \brief Check that all SB9 `Bs*` stabilization blocks annihilate one rigid mode.
 *
 * The helper assembles one quadratic form for each mode-stabilization block and
 * evaluates it on a prescribed field interpolated in the displacement space.
 * The compile-time row checks protect the documented MATLAB row layout:
 * `Bs1,Bs2,Bs3,Bs4 = 1,1,2,3`.
 *
 * \tparam SpaceType Feel++ displacement function space type.
 * \tparam TrialType Trial basis proxy type.
 * \tparam TestType Test basis proxy type.
 * \tparam FieldExprType Feel++ vector expression type for the checked field.
 * \param Uh Displacement function space.
 * \param u Trial basis proxy.
 * \param v Test basis proxy.
 * \param fieldExpr Rigid-mode field interpolated before energy evaluation.
 */
template <typename SpaceType, typename TrialType, typename TestType, typename FieldExprType>
void
checkSb9StabilizationRigidMode( SpaceType const& Uh,
                                TrialType const& u,
                                TestType const& v,
                                FieldExprType const& fieldExpr )
{
    static_assert( Feel::vf::detail::expression_rows_v<decltype( sb9Bs1( u ) )> == 1 );
    static_assert( Feel::vf::detail::expression_rows_v<decltype( sb9Bs2( u ) )> == 1 );
    static_assert( Feel::vf::detail::expression_rows_v<decltype( sb9Bs3( u ) )> == 2 );
    static_assert( Feel::vf::detail::expression_rows_v<decltype( sb9Bs4( u ) )> == 3 );

    auto bs1 = form2( _trial=Uh, _test=Uh );
    bs1 = integrate( _range=elements( Uh->mesh() ),
                     _quad=sb9ThroughThicknessLobatto5(),
                     _expr=component<0, 0>( sb9Bs1( u ) ) *
                           component<0, 0>( sb9Bs1( v ) ) );

    auto bs2 = form2( _trial=Uh, _test=Uh );
    bs2 = integrate( _range=elements( Uh->mesh() ),
                     _quad=sb9ThroughThicknessLobatto5(),
                     _expr=component<0, 0>( sb9Bs2( u ) ) *
                           component<0, 0>( sb9Bs2( v ) ) );

    auto bs3 = form2( _trial=Uh, _test=Uh );
    bs3 = integrate( _range=elements( Uh->mesh() ),
                     _quad=sb9ThroughThicknessLobatto5(),
                     _expr=component<0, 0>( sb9Bs3( u ) ) *
                           component<0, 0>( sb9Bs3( v ) ) +
                           component<1, 0>( sb9Bs3( u ) ) *
                           component<1, 0>( sb9Bs3( v ) ) );

    auto bs4 = form2( _trial=Uh, _test=Uh );
    bs4 = integrate( _range=elements( Uh->mesh() ),
                     _quad=sb9ThroughThicknessLobatto5(),
                     _expr=component<0, 0>( sb9Bs4( u ) ) *
                           component<0, 0>( sb9Bs4( v ) ) +
                           component<1, 0>( sb9Bs4( u ) ) *
                           component<1, 0>( sb9Bs4( v ) ) +
                           component<2, 0>( sb9Bs4( u ) ) *
                           component<2, 0>( sb9Bs4( v ) ) );

    auto uh = Uh->element( "uh" );
    auto vh = Uh->element( "vh" );
    uh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );
    vh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );

    BOOST_CHECK_SMALL( formEnergy( bs1, vh, uh ), g_tol );
    BOOST_CHECK_SMALL( formEnergy( bs2, vh, uh ), g_tol );
    BOOST_CHECK_SMALL( formEnergy( bs3, vh, uh ), g_tol );
    BOOST_CHECK_SMALL( formEnergy( bs4, vh, uh ), g_tol );
}

/**
 * \brief Check rigid translations and rotations on one shell patch.
 *
 * \param mesh Shell mesh used to build the Q1 vector displacement space.
 */
void
checkSb9StabilizationRigidModesOnPatch( mesh_ptrtype const& mesh )
{
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );

    checkSb9StabilizationRigidMode( Uh, u, v, rigidTranslationField() );
    checkSb9StabilizationRigidMode( Uh, u, v, rigidRotationField() );
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_stabilization_suite )

/**
 * \test Verify that `Bs1..Bs4` preserve rigid modes on a flat shell patch.
 */
BOOST_AUTO_TEST_CASE( sb9_stabilization_preserves_rigid_modes_on_flat_patch )
{
    checkSb9StabilizationRigidModesOnPatch( createFlatShellPatch( "sb9_stabilization_flat_patch" ) );
}

/**
 * \test Verify that `Bs1..Bs4` preserve rigid modes after a patch rotation.
 */
BOOST_AUTO_TEST_CASE( sb9_stabilization_preserves_rigid_modes_on_rotated_patch )
{
    checkSb9StabilizationRigidModesOnPatch( createRotatedShellPatch( "sb9_stabilization_rotated_patch" ) );
}

BOOST_AUTO_TEST_SUITE_END()
