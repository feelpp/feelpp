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
 * \file test_sb9pinching.cpp
 * \brief Unit tests for SB9 pinching kinematic operators.
 */

#define BOOST_TEST_MODULE test_sb9pinching
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>

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

BOOST_AUTO_TEST_SUITE_END()
