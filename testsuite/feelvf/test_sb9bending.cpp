/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Feel++ Consortium

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * \file test_sb9bending.cpp
 * \brief Unit tests for SB9 membrane and bending kinematic operators.
 */

#define BOOST_TEST_MODULE test_sb9bending
#include <feel/feelcore/testsuite.hpp>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

namespace
{
/**
 * \brief Check that one field has zero SB9 membrane/bending energy.
 *
 * The helper assembles the SB9 membrane-plus-bending bilinear form and
 * evaluates it on the interpolated field. It is used for rigid translation and
 * rigid rotation fields, both of which must be in the kernel of the bending
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
checkSb9BendingRigidMode( SpaceType const& Uh,
                          TrialType const& u,
                          TestType const& v,
                          FieldExprType const& fieldExpr )
{
    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;
    auto C = isotropic_stiffness<3>( lambda, mu );

    auto bending = form2( _trial=Uh, _test=Uh );
    bending = integrate( _range=elements( Uh->mesh() ),
                         _expr=ddot( C,
                                     sb9MembraneBending( u, zeta() ),
                                     sb9MembraneBending( v, zeta() ) ) );

    auto uh = Uh->element( "uh" );
    auto vh = Uh->element( "vh" );
    uh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );
    vh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );

    BOOST_CHECK_SMALL( formEnergy( bending, vh, uh ), g_tol );
}

/**
 * \brief Check rigid translations and rotations on one shell patch.
 *
 * \param mesh Shell mesh used to build the Q1 vector displacement space.
 */
void
checkSb9BendingRigidModesOnPatch( mesh_ptrtype const& mesh )
{
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );

    checkSb9BendingRigidMode( Uh, u, v, rigidTranslationField() );
    checkSb9BendingRigidMode( Uh, u, v, rigidRotationField() );
}

/**
 * \brief Compute SB9 membrane/bending energy of a linear displacement field.
 *
 * This assembles the SB9 membrane-plus-bending bilinear form with
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
 * \return Quadratic membrane/bending energy of the interpolated field.
 */
template <typename UxExpr, typename UyExpr, typename UzExpr>
double
sb9BendingLinearEnergy( mesh_ptrtype const& mesh,
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
                   _expr=ddot( C, sb9MembraneBending( u, zeta() ), sb9MembraneBending( v, zeta() ) ) );

    auto uh = Uh->element( "uh" );
    uh.on( _range=elements( mesh ), _expr=vec( ux, uy, uz ), _close=true );
    return formEnergy( a, uh, uh );
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_bending_suite )

/**
 * \test Verify that SB9 membrane/bending terms preserve rigid modes on a flat patch.
 */
BOOST_AUTO_TEST_CASE( sb9_bending_preserves_rigid_modes_on_flat_patch )
{
    checkSb9BendingRigidModesOnPatch( createFlatShellPatch( "sb9_bending_flat_patch" ) );
}

/**
 * \test Verify that SB9 membrane/bending terms preserve rigid modes after rotation.
 */
BOOST_AUTO_TEST_CASE( sb9_bending_preserves_rigid_modes_on_rotated_patch )
{
    checkSb9BendingRigidModesOnPatch( createRotatedShellPatch( "sb9_bending_rotated_patch" ) );
}

/**
 * \test Verify simple axis-aligned in-plane extension energies.
 *
 * On the unit patch and with the chosen material constants, unit extension in
 * the x or y direction should both produce energy one.
 */
BOOST_AUTO_TEST_CASE( sb9_bending_linear_x_and_y_have_expected_energy_on_axis_aligned_patch )
{
    auto mesh = createAxisAlignedUnitPatch( "sb9_bending_axis_aligned_unit_patch" );

    double const energyX = sb9BendingLinearEnergy( mesh, Px(), cst( 0.0 ), cst( 0.0 ) );
    double const energyY = sb9BendingLinearEnergy( mesh, cst( 0.0 ), Py(), cst( 0.0 ) );

    BOOST_TEST_MESSAGE( "sb9 bending linear extension energies: Ex=" << energyX << " Ey=" << energyY );
    BOOST_CHECK_CLOSE( energyX, 1.0, 1e-8 );
    BOOST_CHECK_CLOSE( energyY, 1.0, 1e-8 );
}

BOOST_AUTO_TEST_SUITE_END()
