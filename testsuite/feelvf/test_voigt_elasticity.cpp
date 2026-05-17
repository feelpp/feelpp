/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

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

#define BOOST_TEST_MODULE voigt_elasticity testsuite
#include <feel/feelcore/testsuite.hpp>

#include <cmath>

#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

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

template <int Dim>
auto
sampleTrialField()
{
    if constexpr ( Dim == 2 )
    {
        return vec( 2.0*Px() + 3.0*Py(),
                    -Px() + 4.0*Py() );
    }
    else
    {
        return vec( 2.0*Px() + 3.0*Py() - Pz(),
                    -Px() + 4.0*Py() + 2.0*Pz(),
                    Px() - Py() + 5.0*Pz() );
    }
}

template <int Dim>
auto
sampleTestField()
{
    if constexpr ( Dim == 2 )
    {
        return vec( -Px() + Py(),
                    5.0*Px() + 2.0*Py() );
    }
    else
    {
        return vec( -Px() + Py() + 2.0*Pz(),
                    5.0*Px() + 2.0*Py() - Pz(),
                    -2.0*Px() + 3.0*Py() + Pz() );
    }
}

template <int Dim, typename TrialExprT, typename TestExprT>
auto
tensorElasticityForm( TrialExprT const& u, TestExprT const& v, double lambda, double mu )
{
    auto epsu = symm_grad( u );
    auto epsv = symm_grad( v );
    return cst( lambda )*trace( epsu )*trace( epsv ) + cst( 2.0*mu )*inner( epsu, epsv );
}

template <int Dim, typename TrialExprT, typename TestExprT>
auto
mandelElasticityForm( TrialExprT const& u, TestExprT const& v, double lambda, double mu )
{
    auto C = isotropic_stiffness<Dim>( lambda, mu );
    return ddot( C, symm_grad( u ), symm_grad( v ) );
}

template <int Dim, typename TrialExprT, typename TestExprT>
auto
voigtElasticityForm( TrialExprT const& u, TestExprT const& v, double lambda, double mu )
{
    auto C = isotropic_stiffness<Dim, SymmetricTensorNotation::Voigt>( lambda, mu );
    return ddot<SymmetricTensorNotation::Voigt>( C, symm_grad( u ), symm_grad( v ) );
}

template <int Dim>
void
checkElasticityEquivalence()
{
    auto mesh = makeStructuredMesh<Dim>( Dim == 2 ? 8 : 4 );
    auto Xh = Pchv<1>( mesh );
    auto u = trial( Xh, "u" );
    auto v = test( Xh, "v" );

    constexpr double E = 1.3e6;
    constexpr double nu = 0.31;
    constexpr double lambda = E*nu/( ( 1.0 + nu )*( 1.0 - 2.0*nu ) );
    constexpr double mu = E/( 2.0*( 1.0 + nu ) );

    auto aTensor = form2( _trial=Xh, _test=Xh );
    aTensor = integrate( _range=elements( mesh ),
                         _expr=tensorElasticityForm<Dim>( u, v, lambda, mu ) );
    aTensor.close();

    auto aMandel = form2( _trial=Xh, _test=Xh );
    aMandel = integrate( _range=elements( mesh ),
                         _expr=mandelElasticityForm<Dim>( u, v, lambda, mu ) );
    aMandel.close();

    auto aVoigt = form2( _trial=Xh, _test=Xh );
    aVoigt = integrate( _range=elements( mesh ),
                        _expr=voigtElasticityForm<Dim>( u, v, lambda, mu ) );
    aVoigt.close();

    auto uh = Xh->element( "uh" );
    auto vh = Xh->element( "vh" );
    uh.on( _range=elements( mesh ), _expr=sampleTrialField<Dim>(), _close=true );
    vh.on( _range=elements( mesh ), _expr=sampleTestField<Dim>(), _close=true );

    auto tensorValue = aTensor( vh, uh );
    auto mandelValue = aMandel( vh, uh );
    auto voigtValue = aVoigt( vh, uh );
    auto scale = std::abs( tensorValue );
    if ( scale < 1.0 )
        scale = 1.0;

    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelValue )/scale, 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - voigtValue )/scale, 1e-12 );
}

} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( voigt_elasticity_suite )

BOOST_AUTO_TEST_CASE( elasticity_forms_are_equivalent_in_2d )
{
    checkElasticityEquivalence<2>();
}

BOOST_AUTO_TEST_CASE( elasticity_forms_are_equivalent_in_3d )
{
    checkElasticityEquivalence<3>();
}

BOOST_AUTO_TEST_SUITE_END()
