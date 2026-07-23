/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE test_sb9shear
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelvf/sb9_shear.hpp>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

namespace
{
template <typename SpaceType, typename TrialType, typename TestType, typename FieldExprType>
void
checkSb9ShearRigidMode( SpaceType const& Uh,
                        TrialType const& u,
                        TestType const& v,
                        FieldExprType const& fieldExpr )
{
    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;
    auto C = isotropic_stiffness<3>( lambda, mu );

    auto shear = form2( _trial=Uh, _test=Uh );
    shear = integrate( _range=elements( Uh->mesh() ),
                       _expr=ddot( C, sb9Shear( u, cst( 1.0 ) ), sb9Shear( v, cst( 1.0 ) ) ) );

    auto stabilizationC1 = form2( _trial=Uh, _test=Uh );
    stabilizationC1 = integrate( _range=elements( Uh->mesh() ),
                                 _expr=inner( sb9Bc1( u ), sb9Bc1( v ) ) );

    auto stabilizationC2 = form2( _trial=Uh, _test=Uh );
    stabilizationC2 = integrate( _range=elements( Uh->mesh() ),
                                 _expr=inner( sb9Bc2( u ), sb9Bc2( v ) ) );

    auto uh = Uh->element( "uh" );
    auto vh = Uh->element( "vh" );
    uh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );
    vh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );

    BOOST_CHECK_SMALL( formEnergy( shear, vh, uh ), g_tol );
    BOOST_CHECK_SMALL( formEnergy( stabilizationC1, vh, uh ), g_tol );
    BOOST_CHECK_SMALL( formEnergy( stabilizationC2, vh, uh ), g_tol );
}

void
checkSb9ShearRigidModesOnPatch( mesh_ptrtype const& mesh )
{
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );

    checkSb9ShearRigidMode( Uh, u, v, rigidTranslationField() );
    checkSb9ShearRigidMode( Uh, u, v, rigidRotationField() );
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_shear_suite )

BOOST_AUTO_TEST_CASE( sb9_shear_preserves_rigid_modes_on_flat_patch )
{
    checkSb9ShearRigidModesOnPatch( createFlatShellPatch( "sb9_shear_flat_patch" ) );
}

BOOST_AUTO_TEST_CASE( sb9_shear_preserves_rigid_modes_on_rotated_patch )
{
    checkSb9ShearRigidModesOnPatch( createRotatedShellPatch( "sb9_shear_rotated_patch" ) );
}

BOOST_AUTO_TEST_SUITE_END()
