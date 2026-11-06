/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme
       Date: 2026-03-25

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#define BOOST_TEST_MODULE trial_test_basis testsuite
#include <feel/feelcore/testsuite.hpp>

#include <cmath>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

template <typename FormType, typename TestElement, typename TrialElement>
double
formEnergy( FormType& form,
            TestElement const& testElement,
            TrialElement const& trialElement )
{
    form.close();
    return form.matrixPtr()->energy( testElement, trialElement );
}

template <int Dim>
auto
makeStructuredMesh( int n )
{
    using mesh_t = MeshStructured<Hypercube<Dim>>;

    auto discretisation = nl::json::array();
    for ( int axis = 0; axis < Dim; ++axis )
        discretisation.push_back( n );

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", discretisation } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();
    return mesh;
}

template <typename EpsTrialExprT, typename EpsTestExprT>
auto
elasticityEnergyDensity( EpsTrialExprT const& epsu, EpsTestExprT const& epsv )
{
    constexpr double E = 1.0;
    constexpr double nu = 0.45;
    constexpr double lambda = E*nu/( ( 1.0 + nu )*( 1.0 - 2.0*nu ) );
    constexpr double mu = E/( 2.0*( 1.0 + nu ) );

    return cst( lambda )*trace( epsu )*trace( epsv ) + cst( 2.0*mu )*inner( epsu, epsv );
}

template <typename MeshType>
void
check3DSymmGradProxyMatchesExplicitSyntax( std::shared_ptr<MeshType> const& mesh )
{
    auto Xh = Pchv<1>( mesh );

    auto uExplicit = Xh->element( "u_explicit" );
    auto vExplicit = Xh->element( "v_explicit" );
    auto uFromElement = trial( uExplicit );
    auto vFromElement = test( vExplicit );
    auto uFromSpace = trial( Xh, "u_space" );
    auto vFromSpace = test( Xh, "v_space" );

    auto explicitSymGrad = form2( _test=Xh, _trial=Xh );
    explicitSymGrad = integrate( _range=elements( mesh ),
                                 _expr=elasticityEnergyDensity( sym( gradt( uExplicit ) ), sym( grad( vExplicit ) ) ) );

    auto explicitSymmGrad = form2( _test=Xh, _trial=Xh );
    explicitSymmGrad = integrate( _range=elements( mesh ),
                                  _expr=elasticityEnergyDensity( symm_gradt( uExplicit ), symm_grad( vExplicit ) ) );

    auto proxyFromElement = form2( _test=Xh, _trial=Xh );
    proxyFromElement = integrate( _range=elements( mesh ),
                                  _expr=elasticityEnergyDensity( symm_grad( uFromElement ), symm_grad( vFromElement ) ) );

    auto proxyFromSpace = form2( _test=Xh, _trial=Xh );
    proxyFromSpace = integrate( _range=elements( mesh ),
                                _expr=elasticityEnergyDensity( symm_grad( uFromSpace ), symm_grad( vFromSpace ) ) );

    auto trialProbe = Xh->element( vec( Px()+2.0*Py()-Pz(),
                                        -2.0*Px()+Py()+0.5*Pz(),
                                        3.0*Px()-Py()+2.0*Pz() ),
                                   "trial_probe" );
    auto testProbe = Xh->element( vec( cst( 1.0 )+Px()-Py(),
                                       2.0*Px()+Py()*Pz(),
                                       -Px()+0.5*Py()+Pz() ),
                                  "test_probe" );

    double explicitSymGradEnergy = formEnergy( explicitSymGrad, testProbe, trialProbe );
    double explicitSymmGradEnergy = formEnergy( explicitSymmGrad, testProbe, trialProbe );
    double proxyFromElementEnergy = formEnergy( proxyFromElement, testProbe, trialProbe );
    double proxyFromSpaceEnergy = formEnergy( proxyFromSpace, testProbe, trialProbe );
    double scale = std::max( 1.0, std::abs( explicitSymGradEnergy ) );

    BOOST_CHECK_SMALL( std::abs( explicitSymGradEnergy-explicitSymmGradEnergy )/scale, 1e-12 );
    BOOST_CHECK_SMALL( std::abs( explicitSymGradEnergy-proxyFromElementEnergy )/scale, 1e-12 );
    BOOST_CHECK_SMALL( std::abs( explicitSymGradEnergy-proxyFromSpaceEnergy )/scale, 1e-12 );

    BOOST_CHECK_GT( std::abs( explicitSymGradEnergy ), 1e-8 );
}

} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( trial_test_basis_suite )

BOOST_AUTO_TEST_CASE( scalar_basis_proxy_matches_explicit_syntax )
{
    auto mesh = unitSquare();
    auto Xh = Pch<2>( mesh );

    auto uExplicit = Xh->element( "u_explicit" );
    auto vExplicit = Xh->element( "v_explicit" );
    auto u = trial( Xh, "u" );
    auto v = test( Xh, "v" );

    auto explicitForm = form2( _test=Xh, _trial=Xh );
    explicitForm = integrate( _range=elements( mesh ),
                              _expr=inner( gradt( uExplicit ), grad( vExplicit ) ) + idt( uExplicit )*id( vExplicit ) );

    auto modern = form2( _test=Xh, _trial=Xh );
    modern = integrate( _range=elements( mesh ),
                        _expr=inner( grad( u ), grad( v ) ) + u*v );

    auto trialProbe = Xh->element( Px() + 2.0*Py(), "trial_probe" );
    auto testProbe = Xh->element( cst( 1.0 ) + Px()*Py(), "test_probe" );

    double explicitEnergy = formEnergy( explicitForm, testProbe, trialProbe );
    double modernEnergy = formEnergy( modern, testProbe, trialProbe );

    BOOST_CHECK_SMALL( std::abs( explicitEnergy-modernEnergy ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( vector_basis_proxy_matches_explicit_symm_grad )
{
    auto mesh = unitSquare();
    auto Xh = Pchv<2>( mesh );

    auto uExplicit = Xh->element( "u_explicit" );
    auto vExplicit = Xh->element( "v_explicit" );
    auto u = trial( uExplicit );
    auto v = test( vExplicit );

    auto explicitForm = form2( _test=Xh, _trial=Xh );
    explicitForm = integrate( _range=elements( mesh ),
                              _expr=inner( sym( gradt( uExplicit ) ), sym( grad( vExplicit ) ) ) + inner( idt( uExplicit ), id( vExplicit ) ) );

    auto modern = form2( _test=Xh, _trial=Xh );
    modern = integrate( _range=elements( mesh ),
                        _expr=inner( symm_grad( u ), symm_grad( v ) ) + inner( u, v ) );

    auto trialProbe = Xh->element( vec( Px()+Py(), 2.0*Px()-Py() ), "trial_probe" );
    auto testProbe = Xh->element( vec( cst( 1.0 )+Px(), Px()*Py() ), "test_probe" );

    double explicitEnergy = formEnergy( explicitForm, testProbe, trialProbe );
    double modernEnergy = formEnergy( modern, testProbe, trialProbe );

    BOOST_CHECK_SMALL( std::abs( explicitEnergy-modernEnergy ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( vector_basis_proxy_matches_explicit_symm_grad_in_3d_simplex )
{
    check3DSymmGradProxyMatchesExplicitSyntax( unitCube( 0.7 ) );
}

BOOST_AUTO_TEST_CASE( vector_basis_proxy_matches_explicit_symm_grad_in_3d_hypercube )
{
    check3DSymmGradProxyMatchesExplicitSyntax( makeStructuredMesh<3>( 3 ) );
}

BOOST_AUTO_TEST_SUITE_END()
