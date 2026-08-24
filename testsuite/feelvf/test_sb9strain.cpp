/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE test_sb9strain
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelvf/sb9_strain.hpp>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

namespace
{
template <typename UxExpr, typename UyExpr, typename UzExpr>
double
sb9ShellStrainLinearEnergy( mesh_ptrtype const& mesh,
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
                   _expr=ddot( C,
                               sb9ShellStrain( u, zeta(), cst( 0.0 ) ),
                               sb9ShellStrain( v, zeta(), cst( 0.0 ) ) ) );

    auto uh = Uh->element( "uh" );
    uh.on( _range=elements( mesh ), _expr=vec( ux, uy, uz ), _close=true );
    return formEnergy( a, uh, uh );
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_strain_suite )

BOOST_AUTO_TEST_CASE( sb9_shell_strain_linear_x_and_y_have_same_energy_without_shear_on_axis_aligned_patch )
{
    auto mesh = createAxisAlignedUnitPatch( "sb9_strain_axis_aligned_unit_patch" );

    double const energyX = sb9ShellStrainLinearEnergy( mesh, Px(), cst( 0.0 ), cst( 0.0 ) );
    double const energyY = sb9ShellStrainLinearEnergy( mesh, cst( 0.0 ), Py(), cst( 0.0 ) );
    double const energyZ = sb9ShellStrainLinearEnergy( mesh, cst( 0.0 ), cst( 0.0 ), Pz() + cst( 0.5 ) );

    BOOST_TEST_MESSAGE( "sb9 shell strain extension energies: Ex=" << energyX << " Ey=" << energyY << " Ez=" << energyZ );
    BOOST_CHECK_CLOSE( energyX, energyY, 1e-8 );
}

BOOST_AUTO_TEST_SUITE_END()
