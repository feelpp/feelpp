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

// TEST FORME BILINEAIRE PAR BLOC (25*25)
/**
 * \brief Check that one field has zero SB9 pinching energy.
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
template <typename ProductSpaceType, typename ProductTrialType, typename ProductTestType, typename FieldExprType>
void
checkSb9PinchingRigidMode( mesh_ptrtype const& mesh,
                          ProductSpaceType &ps,
                          ProductTrialType &u,
                          ProductTestType &v,
                          FieldExprType const& fieldExpr ) // mettre les espaces en plus et uw, vw
{
    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;
    auto C = isotropic_stiffness<3>( lambda, mu );

    auto Uh = Pchv<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    // auto u  = trial( Uh, "u" );
    // auto v  = test( Uh, "v" );
    auto uw = trial( Ah, "uw" );
    auto vw = test( Ah, "vw" );


    // auto frame = shellFrame();
    // std::cout << "Size du frame = " << frame.size() << std::endl;


    backend( _rebuild=true );
    auto pinching = blockform2( ps, solve::strategy::monolithic, backend() );
    pinching( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                      _expr=ddot( C, sb9Pinching( u, zeta() ), sb9Pinching( v, zeta() ) ) );
    pinching( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                      _expr=ddot( C, sb9Pinching( u, zeta() ), sb9PinchingW9( vw, zeta() ) ) );
    pinching( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                      _expr=ddot( C, sb9PinchingW9( uw, zeta() ), sb9Pinching( v, zeta() ) ) );
    pinching( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                      _expr=ddot( C, sb9PinchingW9( uw, zeta() ), sb9PinchingW9( vw, zeta() ) ) );
                                    

    auto uh = ps.element();
    auto vh = ps.element();

    // auto fe = vec( cst(0.0), cst(0.0), cst(1.0) );
    auto fe = vec( cst(0.0), cst(0.0), Px()*Px() );

    // uh( 0_c ).on( _range=elements( mesh ), _expr=fieldExpr, _close=true );
    // vh( 0_c ).on( _range=elements( mesh ), _expr=fieldExpr, _close=true );
    // uh( 1_c ).on( _range=elements( mesh ), _expr=cst( 0.0 ), _close=true );
    // vh( 1_c ).on( _range=elements( mesh ), _expr=cst( 0.0 ), _close=true );
    uh( 0_c ).on( _range=elements( mesh ), _expr=fe, _close=true );
    vh( 0_c ).on( _range=elements( mesh ), _expr=fe, _close=true );
    uh( 1_c ).on( _range=elements( mesh ), _expr= Px()*Px(), _close=true );
    vh( 1_c ).on( _range=elements( mesh ), _expr= Px()*Px(), _close=true );

    pinching( 0_c, 0_c ).matrix().printMatlab( "Pinching25_form2.m" );
    uh( 0_c ).printMatlab( "Pinching24_deplacement.m" );
    uh( 1_c ).printMatlab( "Pinching25eme_deplacement.m" );

    uh.buildVector( backend() );
    vh.buildVector( backend() );
    uh.updateVectorFromSubVectors();
    vh.updateVectorFromSubVectors();

    double energy = formEnergy( pinching, vh.vectorMonolithic(), uh.vectorMonolithic() );
    BOOST_TEST_MESSAGE( "sb9 pinching rigid mode energy = " << energy );
    BOOST_CHECK_SMALL( energy, g_tol );
}

/**
 * \brief Check rigid translations and rotations on one shell patch.
 *
 * \param mesh Shell mesh used to build the Q1 vector displacement space.
 */
void
checkSb9PinchingRigidModesOnPatch( mesh_ptrtype const& mesh )
{
    auto Uh = Pchv<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto ps = product( Uh, Ah );

    auto u  = trial( Uh, "u" );
    auto v  = test( Uh, "v" );

    // auto u  = trial(Uh, "u");
    // auto v  = test(Uh, "v");
    // auto uw = trial(Ah, "uw");
    // auto vw = test(Ah, "vw");

    // auto U = ps.element();
    // auto V = ps.element();
    // auto U = Uh->element();
    // auto V = Ah->element();

    checkSb9PinchingRigidMode( mesh, ps, u, v, rigidTranslationField() );
    // checkSb9PinchingRigidMode( mesh, ps, u, v, rigidRotationField() );
}



// // TEST COUPLAGE Uh,Uh
// template <typename ProductSpaceType, typename ProductTrialType, typename ProductTestType, typename FieldExprType>
// void
// checkSb9PinchingRigidMode( mesh_ptrtype const& mesh,
//                           ProductSpaceType &ps,
//                           ProductTrialType &u,
//                           ProductTestType &v,
//                           FieldExprType const& fieldExpr ) // mettre les espaces en plus et uw, vw
// {
//     constexpr double lambda = 2.3;
//     constexpr double mu = 1.1;
//     auto C = isotropic_stiffness<3>( lambda, mu );

//     auto Uh = Pchv<1>( mesh );
//     // auto Ah = Pdh<0>( mesh );
//     // auto u  = trial( Uh, "u" );
//     // auto v  = test( Uh, "v" );
//     auto uw = trial( Uh, "uw" );
//     auto vw = test( Uh, "vw" );


//     backend( _rebuild=true );
//     auto pinching = blockform2( ps, solve::strategy::monolithic, backend() );
//     pinching( 0_c, 0_c ) += integrate( _range=elements( mesh ),
//                                       _expr=ddot( C, sb9Pinching( u, zeta() ), sb9Pinching( v, zeta() ) ) );
//     pinching( 0_c, 1_c ) += integrate( _range=elements( mesh ),
//                                       _expr=ddot( C, sb9Pinching( u, zeta() ), sb9Pinching( vw, zeta() ) ) );
//     pinching( 1_c, 0_c ) += integrate( _range=elements( mesh ),
//                                       _expr=ddot( C, sb9Pinching( uw, zeta() ), sb9Pinching( v, zeta() ) ) );
//     pinching( 1_c, 1_c ) += integrate( _range=elements( mesh ),
//                                       _expr=ddot( C, sb9Pinching( uw, zeta() ), sb9Pinching( vw, zeta() ) ) );
                                    

//     auto uh = ps.element();
//     auto vh = ps.element();

//     auto fe = vec( cst(0.0), cst(0.0), cst(0.0) );

//     // uh( 0_c ).on( _range=elements( mesh ), _expr=fieldExpr, _close=true );
//     // vh( 0_c ).on( _range=elements( mesh ), _expr=fieldExpr, _close=true );
//     // uh( 1_c ).on( _range=elements( mesh ), _expr=cst( 0.0 ), _close=true );
//     // vh( 1_c ).on( _range=elements( mesh ), _expr=cst( 0.0 ), _close=true );
//     uh( 0_c ).on( _range=elements( mesh ), _expr=fe, _close=true );
//     vh( 0_c ).on( _range=elements( mesh ), _expr=fe, _close=true );
//     uh( 1_c ).on( _range=elements( mesh ), _expr=vec( cst(1.0), cst(1.0), cst(1.0) ), _close=true );
//     vh( 1_c ).on( _range=elements( mesh ), _expr=vec( cst(1.0), cst(1.0), cst(1.0) ), _close=true );

//     pinching( 0_c, 0_c ).matrix().printMatlab( "Pinching25_form2.m" );
//     uh( 0_c ).printMatlab( "Pinching24_deplacement.m" );
//     uh( 1_c ).printMatlab( "Pinching25eme_deplacement.m" );

//     uh.buildVector( backend() );
//     vh.buildVector( backend() );
//     uh.updateVectorFromSubVectors();
//     vh.updateVectorFromSubVectors();

//     double energy = formEnergy( pinching, vh.vectorMonolithic(), uh.vectorMonolithic() );
//     BOOST_TEST_MESSAGE( "sb9 pinching rigid mode energy = " << energy );
//     BOOST_CHECK_SMALL( energy, g_tol );
// }

// /**
//  * \brief Check rigid translations and rotations on one shell patch.
//  *
//  * \param mesh Shell mesh used to build the Q1 vector displacement space.
//  */
// void
// checkSb9PinchingRigidModesOnPatch( mesh_ptrtype const& mesh )
// {
//     auto Uh = Pchv<1>( mesh );
//     // auto Ah = Pdh<0>( mesh );
//     auto ps = product( Uh, Uh );

//     auto u  = trial( Uh, "u" );
//     auto v  = test( Uh, "v" );

//     // auto u  = trial(Uh, "u");
//     // auto v  = test(Uh, "v");
//     // auto uw = trial(Ah, "uw");
//     // auto vw = test(Ah, "vw");

//     // auto U = ps.element();
//     // auto V = ps.element();
//     // auto U = Uh->element();
//     // auto V = Ah->element();

//     checkSb9PinchingRigidMode( mesh, ps, u, v, rigidTranslationField() );
//     // checkSb9PinchingRigidMode( mesh, ps, u, v, rigidRotationField() );
// }





// // ----------------------------------------------------------------------------
// // TEST BLOC PRINCIPAL (24*24)
// template <typename SpaceType, typename TrialType, typename TestType, typename FieldExprType>
// void checkSb9PinchingRigidMode( SpaceType const& Uh,
//                                 TrialType const& u,
//                                 TestType const& v,
//                                 FieldExprType const& fieldExpr )
// {
//     constexpr double lambda = 2.3;
//     constexpr double mu = 1.1;
//     auto C = isotropic_stiffness<3>( lambda, mu );

//     auto pinching = form2( _trial=Uh, _test=Uh );
//     pinching = integrate( _range=elements( Uh->mesh() ),
//                          _expr=ddot( C,
//                                      sb9Pinching( u, zeta() ),
//                                      sb9Pinching( v, zeta() ) ) );

//     auto uh = Uh->element( "uh" );
//     auto vh = Uh->element( "vh" );                             
//     // uh.on( _range=elements( Uh->mesh() ), _expr=vec( cst(1.0), cst(2.0), cst(3.0) ), _close=true );   // fieldExpr
//     // vh.on( _range=elements( Uh->mesh() ), _expr=vec( cst(1.0), cst(2.0), cst(3.0) ), _close=true );
//     uh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );   // fieldExpr
//     vh.on( _range=elements( Uh->mesh() ), _expr=fieldExpr, _close=true );

//     pinching.matrix().printMatlab( "A_testsb9Pinching_form2.m" );
//     uh.printMatlab( "A_testsb9Pinching_deplacement.m" );

//     auto energy = formEnergy( pinching, uh, vh );

//     BOOST_TEST_MESSAGE( "sb9 pinching energies: E = " << energy );

//     BOOST_CHECK_SMALL( energy, g_tol );
// }

// /**
//  * \brief Check rigid translations and rotations on one shell patch.
//  *
//  * \param mesh Shell mesh used to build the Q1 vector displacement space.
//  */
// void
// checkSb9PinchingRigidModesOnPatch( mesh_ptrtype const& mesh )
// {
//     auto Uh = Pchv<1>( mesh );
//     auto u  = trial(Uh, "u");
//     auto v  = test(Uh, "v");

//     checkSb9PinchingRigidMode( Uh, u, v, rigidTranslationField() );
//     // checkSb9PinchingRigidMode( Uh, u, v, rigidRotationField() );
// }
// // ----------------------------------------------------------------------------






// // ----------------------------------------------------------------------------
// // TEST BLOC SECONDAIRE (1*1)
// template <typename SpaceType, typename TrialType, typename TestType, typename FieldExprType>
// void checkSb9PinchingRigidMode( SpaceType const& Ah,
//                                 TrialType const& u,
//                                 TestType const& v,
//                                 FieldExprType const& fieldExpr )
// {
//     constexpr double lambda = 2.3;
//     constexpr double mu = 1.1;
//     auto C = isotropic_stiffness<3>( lambda, mu );


//     // auto Pinchw9_term_u = mandel_vec<3>( cst(0.0), cst(0.0), - cst(4.0)* shellThickness() * id( u ), cst(0.0), cst(0.0), cst(0.0) );
//     // auto Pinchw9_term_v = mandel_vec<3>( cst(0.0), cst(0.0), - cst(4.0)* shellThickness() * id( v ), cst(0.0), cst(0.0), cst(0.0) );


//     auto pinching = form2( _trial=Ah, _test=Ah );
//     pinching = integrate( _range=elements( Ah->mesh() ),
//                          _expr=ddot( C,
//                                      sb9PinchingW9( u, zeta() ),
//                                      sb9PinchingW9( v, zeta() ) ) );
//     // pinching = integrate( _range=elements( Ah->mesh() ),
//     //                      _expr=ddot( C, Pinchw9_term_u, Pinchw9_term_v ) );


//     auto uh = Ah->element( "uh" );
//     auto vh = Ah->element( "vh" );                             
//     uh.on( _range=elements( Ah->mesh() ), _expr=cst(1.0), _close=true );
//     vh.on( _range=elements( Ah->mesh() ), _expr=cst(1.0), _close=true );

//     pinching.matrix().printMatlab( "heho_testsb9Pinchingw9_form2.m" );
//     uh.printMatlab( "heho_testsb9Pinchingw9_deplacement.m" );

//     auto energy = formEnergy( pinching, uh, vh );

//     BOOST_TEST_MESSAGE( "sb9 pinching energies: E = " << energy );

//     BOOST_CHECK_SMALL( energy, g_tol );


// }

// void checkSb9PinchingRigidModesOnPatch( mesh_ptrtype const& mesh )
// {
//     auto Ah = Pdh<0>( mesh );
//     auto u  = trial(Ah, "u");
//     auto v  = test(Ah, "v");

//     checkSb9PinchingRigidMode( Ah, u, v, rigidTranslationField() );
//     checkSb9PinchingRigidMode( Ah, u, v, rigidRotationField() );
// }

// // ----------------------------------------------------------------------------








} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_pinching_suite )

/**
 * \test Verify that SB9 membrane/bending terms preserve rigid modes on a flat patch.
 */
BOOST_AUTO_TEST_CASE( sb9_pinching_preserves_rigid_modes_on_flat_patch )
{
    // checkSb9PinchingRigidModesOnPatch( createFlatShellPatch( "sb9_pinching_flat_patch" ) );
    checkSb9PinchingRigidModesOnPatch( createAxisAlignedUnitPatch( "sb9_pinching_flat_patch" ) );
}

/**
 * \test Verify that SB9 membrane/bending terms preserve rigid modes after rotation.
 */
BOOST_AUTO_TEST_CASE( sb9_pinching_preserves_rigid_modes_on_rotated_patch )
{
    checkSb9PinchingRigidModesOnPatch( createRotatedShellPatch( "sb9_pinching_rotated_patch" ) );
}

BOOST_AUTO_TEST_SUITE_END()
