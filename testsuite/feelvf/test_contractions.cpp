/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme
       Date: 2026-03-26

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

#define BOOST_TEST_MODULE contractions testsuite
#include <feel/feelcore/testsuite.hpp>

#include <cmath>

#include <feel/feeldiscr/meshstructured.hpp>
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

template <typename MeshT, typename ExprT>
double
integratedSquaredNorm( std::shared_ptr<MeshT> const& mesh, ExprT const& expr )
{
    return integrate( _range=elements( mesh ), _expr=inner( expr, expr ) ).evaluate()( 0, 0 );
}

template <typename MeshT, typename ExprT>
double
integratedValue( std::shared_ptr<MeshT> const& mesh, ExprT const& expr )
{
    return integrate( _range=elements( mesh ), _expr=expr ).evaluate()( 0, 0 );
}

template <int Dim>
auto
sampleSymmetricTensor()
{
    if constexpr ( Dim == 2 )
    {
        return mat<2, 2>( 2.0*Px() + 3.0*Py(), -Px() + 4.0*Py(),
                          -Px() + 4.0*Py(), 5.0*Px() - 2.0*Py() );
    }
    else
    {
        return mat<3, 3>( 2.0*Px() + 3.0*Py() - Pz(), -Px() + 4.0*Py() + 2.0*Pz(), 3.0*Px() - Py() + Pz(),
                          -Px() + 4.0*Py() + 2.0*Pz(), 5.0*Px() - 2.0*Py() + 4.0*Pz(), 2.0*Px() + Py() - 3.0*Pz(),
                          3.0*Px() - Py() + Pz(), 2.0*Px() + Py() - 3.0*Pz(), -Px() + 2.0*Py() + 6.0*Pz() );
    }
}

template <int Dim>
auto
sampleSymmetricTensor2()
{
    if constexpr ( Dim == 2 )
    {
        return mat<2, 2>( -Px() + 2.0*Py(), 3.0*Px() + Py(),
                          3.0*Px() + Py(), 4.0*Px() - Py() );
    }
    else
    {
        return mat<3, 3>( -Px() + 2.0*Py() + Pz(), 3.0*Px() + Py() - Pz(), Px() - 2.0*Py() + 4.0*Pz(),
                          3.0*Px() + Py() - Pz(), 4.0*Px() - Py() + 2.0*Pz(), -2.0*Px() + 3.0*Py() + Pz(),
                          Px() - 2.0*Py() + 4.0*Pz(), -2.0*Px() + 3.0*Py() + Pz(), 5.0*Px() + Py() - 3.0*Pz() );
    }
}

template <int Dim, typename TensorExprT>
auto
tensorStress( double lambda, double mu, TensorExprT const& eps )
{
    return cst( lambda )*trace( eps )*eye<Dim, Dim>() + cst( 2.0*mu )*eps;
}

struct Orthotropic2DCoefficients
{
    static constexpr double c0000 = 8.0;
    static constexpr double c1111 = 6.0;
    static constexpr double c0011 = 1.5;
    static constexpr double c0001 = 0.6;
    static constexpr double c1101 = 0.4;
    static constexpr double c0101 = 2.2;
};

template <typename TensorExprT>
auto
orthotropicTensorStress2D( TensorExprT const& eps )
{
    return mat<2, 2>(
        cst( Orthotropic2DCoefficients::c0000 ) * component<0, 0>( eps ) +
        cst( Orthotropic2DCoefficients::c0011 ) * component<1, 1>( eps ) +
        cst( 2.0 * Orthotropic2DCoefficients::c0001 ) * component<0, 1>( eps ),
        cst( Orthotropic2DCoefficients::c0001 ) * component<0, 0>( eps ) +
        cst( Orthotropic2DCoefficients::c1101 ) * component<1, 1>( eps ) +
        cst( 2.0 * Orthotropic2DCoefficients::c0101 ) * component<0, 1>( eps ),
        cst( Orthotropic2DCoefficients::c0001 ) * component<0, 0>( eps ) +
        cst( Orthotropic2DCoefficients::c1101 ) * component<1, 1>( eps ) +
        cst( 2.0 * Orthotropic2DCoefficients::c0101 ) * component<0, 1>( eps ),
        cst( Orthotropic2DCoefficients::c0011 ) * component<0, 0>( eps ) +
        cst( Orthotropic2DCoefficients::c1111 ) * component<1, 1>( eps ) +
        cst( 2.0 * Orthotropic2DCoefficients::c1101 ) * component<0, 1>( eps ) );
}

auto
orthotropicMandelMatrix2D()
{
    constexpr double sqrt2 = std::numbers::sqrt2_v<double>;

    return mat<3, 3>(
        cst( Orthotropic2DCoefficients::c0000 ),
        cst( sqrt2 * Orthotropic2DCoefficients::c0001 ),
        cst( Orthotropic2DCoefficients::c0011 ),
        cst( sqrt2 * Orthotropic2DCoefficients::c0001 ),
        cst( 2.0 * Orthotropic2DCoefficients::c0101 ),
        cst( sqrt2 * Orthotropic2DCoefficients::c1101 ),
        cst( Orthotropic2DCoefficients::c0011 ),
        cst( sqrt2 * Orthotropic2DCoefficients::c1101 ),
        cst( Orthotropic2DCoefficients::c1111 ) );
}

auto
orthotropicVoigtMatrix2D()
{
    return mat<3, 3>(
        cst( Orthotropic2DCoefficients::c0000 ),
        cst( 2.0 * Orthotropic2DCoefficients::c0001 ),
        cst( Orthotropic2DCoefficients::c0011 ),
        cst( Orthotropic2DCoefficients::c0001 ),
        cst( 2.0 * Orthotropic2DCoefficients::c0101 ),
        cst( Orthotropic2DCoefficients::c1101 ),
        cst( Orthotropic2DCoefficients::c0011 ),
        cst( 2.0 * Orthotropic2DCoefficients::c1101 ),
        cst( Orthotropic2DCoefficients::c1111 ) );
}

template <int Dim>
void
checkConstitutiveTensorAction()
{
    auto mesh = makeStructuredMesh<Dim>( Dim == 2 ? 8 : 4 );
    constexpr double lambda = 2.7;
    constexpr double mu = 1.4;
    auto eps = sampleSymmetricTensor<Dim>();

    auto Cmandel = isotropic_stiffness<Dim>( lambda, mu );
    auto Cvoigt = isotropic_stiffness<Dim, SymmetricTensorNotation::Voigt>( lambda, mu );
    auto expected = tensorStress<Dim>( lambda, mu, eps );
    auto epsMandel = mandel( eps );
    auto epsVoigt = voigt( eps );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, contract( Cmandel, eps ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, ddot( Cmandel, eps ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, double_contract( Cmandel, eps ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, contract<SymmetricTensorNotation::Voigt>( Cvoigt, eps ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, ddot<SymmetricTensorNotation::Voigt>( Cvoigt, eps ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, voigt_contract( Cvoigt, eps ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, unmandel( contract( Cmandel, epsMandel ) ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, unmandel( ddot( Cmandel, epsMandel ) ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, unmandel( double_contract( Cmandel, epsMandel ) ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, unvoigt( contract<SymmetricTensorNotation::Voigt>( Cvoigt, epsVoigt ) ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, unvoigt( ddot<SymmetricTensorNotation::Voigt>( Cvoigt, epsVoigt ) ) - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, unvoigt( voigt_contract( Cvoigt, epsVoigt ) ) - expected ), 1e-12 );
}

template <int Dim>
void
checkConstitutiveBilinearAction()
{
    auto mesh = makeStructuredMesh<Dim>( Dim == 2 ? 8 : 4 );
    constexpr double lambda = 2.1;
    constexpr double mu = 0.9;
    auto eps = sampleSymmetricTensor<Dim>();
    auto eta = sampleSymmetricTensor2<Dim>();

    auto Cmandel = isotropic_stiffness<Dim>( lambda, mu );
    auto Cvoigt = isotropic_stiffness<Dim, SymmetricTensorNotation::Voigt>( lambda, mu );
    auto expected = tensorStress<Dim>( lambda, mu, eps );
    auto epsMandel = mandel( eps );
    auto etaMandel = mandel( eta );
    auto epsVoigt = voigt( eps );
    auto etaVoigt = voigt( eta );
    auto tensorValue = integratedValue( mesh, inner( expected, eta ) );
    auto mandelContractValue = integratedValue( mesh, contract( Cmandel, eps, eta ) );
    auto mandelDdotValue = integratedValue( mesh, ddot( Cmandel, eps, eta ) );
    auto mandelDoubleContractValue = integratedValue( mesh, double_contract( Cmandel, eps, eta ) );
    auto voigtValue = integratedValue( mesh, contract<SymmetricTensorNotation::Voigt>( Cvoigt, eps, eta ) );
    auto voigtDdotValue = integratedValue( mesh, ddot<SymmetricTensorNotation::Voigt>( Cvoigt, eps, eta ) );
    auto voigtAliasValue = integratedValue( mesh, voigt_contract( Cvoigt, eps, eta ) );
    auto mandelStorageValue = integratedValue( mesh, contract( Cmandel, epsMandel, etaMandel ) );
    auto mandelStorageDdotValue = integratedValue( mesh, ddot( Cmandel, epsMandel, etaMandel ) );
    auto mandelStorageDoubleContractValue = integratedValue( mesh, double_contract( Cmandel, epsMandel, etaMandel ) );
    auto voigtStorageValue = integratedValue( mesh, contract<SymmetricTensorNotation::Voigt>( Cvoigt, epsVoigt, etaVoigt ) );
    auto voigtStorageDdotValue = integratedValue( mesh, ddot<SymmetricTensorNotation::Voigt>( Cvoigt, epsVoigt, etaVoigt ) );
    auto voigtStorageAliasValue = integratedValue( mesh, voigt_contract( Cvoigt, epsVoigt, etaVoigt ) );

    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelContractValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelDdotValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelDoubleContractValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - voigtValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - voigtDdotValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( voigtValue - voigtAliasValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelStorageValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelStorageDdotValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelStorageDoubleContractValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - voigtStorageValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - voigtStorageDdotValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( voigtStorageValue - voigtStorageAliasValue ), 1e-12 );
}

void
checkStorageComponentShellLikeAction()
{
    auto mesh = makeStructuredMesh<3>( 4 );
    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;
    auto Cmandel = isotropic_stiffness<3>( lambda, mu );
    auto Cvoigt = isotropic_stiffness<3, SymmetricTensorNotation::Voigt>( lambda, mu );

    auto gamma33 = 1.0 + Px() - 2.0*Py() + 0.5*Pz();
    auto gamma13 = -2.0 + 3.0*Px() + Py() - Pz();

    auto eps33Tensor = mat<3, 3>( cst( 0.0 ), cst( 0.0 ), cst( 0.0 ),
                                  cst( 0.0 ), cst( 0.0 ), cst( 0.0 ),
                                  cst( 0.0 ), cst( 0.0 ), gamma33 );
    auto eps13Tensor = mat<3, 3>( cst( 0.0 ), cst( 0.0 ), gamma13,
                                  cst( 0.0 ), cst( 0.0 ), cst( 0.0 ),
                                  gamma13, cst( 0.0 ), cst( 0.0 ) );
    auto eps33Mandel = mandel_component<3, 2, 2>( gamma33 );
    auto eps13Mandel = mandel_component<3, 0, 2>( gamma13 );
    auto eps33Voigt = voigt_component<3, 2, 2>( gamma33 );
    auto eps13Voigt = voigt_component<3, 0, 2>( gamma13 );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh,
                                              unmandel( ddot( Cmandel, eps33Mandel ) ) -
                                              tensorStress<3>( lambda, mu, eps33Tensor ) ),
                       1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh,
                                              unmandel( ddot( Cmandel, eps13Mandel ) ) -
                                              tensorStress<3>( lambda, mu, eps13Tensor ) ),
                       1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh,
                                              unvoigt( ddot<SymmetricTensorNotation::Voigt>( Cvoigt, eps33Voigt ) ) -
                                              tensorStress<3>( lambda, mu, eps33Tensor ) ),
                       1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh,
                                              unvoigt( ddot<SymmetricTensorNotation::Voigt>( Cvoigt, eps13Voigt ) ) -
                                              tensorStress<3>( lambda, mu, eps13Tensor ) ),
                       1e-12 );
    BOOST_CHECK_SMALL( std::abs( integratedValue( mesh, ddot( Cmandel, eps33Mandel, eps13Mandel ) ) -
                                 integratedValue( mesh, inner( tensorStress<3>( lambda, mu, eps33Tensor ),
                                                               eps13Tensor ) ) ),
                       1e-12 );
}

void
checkOrthotropic2DConstitutiveBilinearAction()
{
    auto mesh = makeStructuredMesh<2>( 8 );
    auto eps = sampleSymmetricTensor<2>();
    auto eta = sampleSymmetricTensor2<2>();
    auto expected = orthotropicTensorStress2D( eps );
    auto epsMandel = mandel( eps );
    auto etaMandel = mandel( eta );
    auto epsVoigt = voigt( eps );
    auto etaVoigt = voigt( eta );

    auto tensorValue = integratedValue( mesh, inner( expected, eta ) );
    auto mandelStress = contract( orthotropicMandelMatrix2D(), eps );
    auto voigtStress = contract<SymmetricTensorNotation::Voigt>( orthotropicVoigtMatrix2D(), eps );
    auto mandelValue = integratedValue( mesh, ddot( orthotropicMandelMatrix2D(), eps, eta ) );
    auto voigtValue = integratedValue( mesh, ddot<SymmetricTensorNotation::Voigt>( orthotropicVoigtMatrix2D(),
                                                                                   eps,
                                                                                   eta ) );
    auto mandelStorageValue = integratedValue( mesh, ddot( orthotropicMandelMatrix2D(), epsMandel, etaMandel ) );
    auto voigtStorageValue = integratedValue( mesh, ddot<SymmetricTensorNotation::Voigt>( orthotropicVoigtMatrix2D(),
                                                                                          epsVoigt,
                                                                                          etaVoigt ) );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, mandelStress - expected ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, voigtStress - expected ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - voigtValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - mandelStorageValue ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( tensorValue - voigtStorageValue ), 1e-12 );
}

} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( contractions_suite )

BOOST_AUTO_TEST_CASE( constitutive_tensor_action_matches_isotropic_tensor_formula )
{
    checkConstitutiveTensorAction<2>();
    checkConstitutiveTensorAction<3>();
}

BOOST_AUTO_TEST_CASE( constitutive_bilinear_action_matches_tensor_inner_product )
{
    checkConstitutiveBilinearAction<2>();
    checkConstitutiveBilinearAction<3>();
}

BOOST_AUTO_TEST_CASE( constitutive_action_accepts_storage_component_builders )
{
    checkStorageComponentShellLikeAction();
}

BOOST_AUTO_TEST_CASE( constitutive_action_matches_orthotropic_tensor_formula )
{
    checkOrthotropic2DConstitutiveBilinearAction();
}

BOOST_AUTO_TEST_CASE( constitutive_action_constant_evaluate_matches_tensor_formula )
{
    constexpr double lambda = 2.0;
    constexpr double mu = 3.0;

    auto eps = mat<2, 2>( cst( 1.0 ), cst( 4.0 ),
                          cst( 4.0 ), cst( 2.0 ) );
    auto eta = mat<2, 2>( cst( 2.0 ), cst( 5.0 ),
                          cst( 5.0 ), cst( -1.0 ) );

    auto expected = tensorStress<2>( lambda, mu, eps ).evaluate();
    auto mandelStress = contract( isotropic_stiffness<2>( lambda, mu ), eps ).evaluate();
    auto voigtStress = contract<SymmetricTensorNotation::Voigt>( isotropic_stiffness<2, SymmetricTensorNotation::Voigt>( lambda, mu ),
                                                                 eps ).evaluate();

    for ( int i = 0; i < expected.rows(); ++i )
        for ( int j = 0; j < expected.cols(); ++j )
        {
            BOOST_CHECK_CLOSE( mandelStress( i, j ), expected( i, j ), 1e-10 );
            BOOST_CHECK_CLOSE( voigtStress( i, j ), expected( i, j ), 1e-10 );
        }

    auto expectedBilinear = inner( tensorStress<2>( lambda, mu, eps ), eta ).evaluate()( 0, 0 );
    auto mandelBilinear = ddot( isotropic_stiffness<2>( lambda, mu ), eps, eta ).evaluate()( 0, 0 );
    auto voigtBilinear = ddot<SymmetricTensorNotation::Voigt>( isotropic_stiffness<2, SymmetricTensorNotation::Voigt>( lambda, mu ),
                                                               eps,
                                                               eta ).evaluate()( 0, 0 );
    auto voigtStressInner = inner( contract<SymmetricTensorNotation::Voigt>( isotropic_stiffness<2, SymmetricTensorNotation::Voigt>( lambda, mu ),
                                                                             eps ),
                                   eta ).evaluate()( 0, 0 );

    BOOST_CHECK_CLOSE( mandelBilinear, expectedBilinear, 1e-10 );
    BOOST_CHECK_CLOSE( voigtBilinear, expectedBilinear, 1e-10 );
    BOOST_CHECK_CLOSE( voigtBilinear, voigtStressInner, 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END()
