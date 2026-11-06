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

#define BOOST_TEST_MODULE tensor_basis testsuite
#include <feel/feelcore/testsuite.hpp>

#include <cmath>
#include <numbers>

#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

template <typename MeshType, typename ExprType>
double
integratedValue( std::shared_ptr<MeshType> const& mesh, ExprType const& expr )
{
    return integrate( _range=elements( mesh ), _expr=expr ).evaluate()( 0, 0 );
}

template <typename MeshType, typename ExprType>
double
integratedAbsError( std::shared_ptr<MeshType> const& mesh, ExprType const& expr )
{
    return std::abs( integratedValue( mesh, expr ) );
}

template <typename MeshType, typename ExprType>
double
integratedSquaredNorm( std::shared_ptr<MeshType> const& mesh, ExprType const& expr )
{
    return integrate( _range=elements( mesh ), _expr=inner( expr, expr ) ).evaluate()( 0, 0 );
}

auto tensorExpr2D()
{
    return mat<2, 2>( Px() + cst( 1.0 ),
                      Py(),
                      Px()*Py(),
                      cst( 2.0 ) + Px() );
}

} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( tensor_basis_suite )

BOOST_AUTO_TEST_CASE( canonical_tensor_basis_rewrites_matrix_component_access )
{
    auto mesh = unitSquare();
    auto A = tensorExpr2D();

    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, delta<2, 0, 1>() ) - A( 0, 1 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, trans( delta<2, 0, 1>() ) ) - A( 1, 0 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, sym( delta<2, 0, 1>() ) ) -
                                                 cst( 0.5 )*( A( 0, 1 ) + A( 1, 0 ) ) ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( symmetric_and_mandel_tensor_bases_preserve_expected_scaling )
{
    auto mesh = unitSquare();
    auto A = tensorExpr2D();
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, symm_delta<2, 0, 1>() ) -
                                                 ( A( 0, 1 ) + A( 1, 0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, symm_delta<2, 1, 1>() ) - A( 1, 1 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, mandel_delta<2, 0, 1>() ) -
                                                 cst( invSqrt2 )*( A( 0, 1 ) + A( 1, 0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, mandel_delta<2, 1, 1>() ) - A( 1, 1 ) ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( tensor_basis_inner_products_match_their_structured_norms )
{
    auto mesh = unitSquare();

    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( delta<2, 0, 1>(), delta<2, 0, 1>() ) - cst( 1.0 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( delta<2, 0, 1>(), symm_delta<2, 0, 1>() ) - cst( 1.0 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( symm_delta<2, 0, 1>(), symm_delta<2, 0, 1>() ) - cst( 2.0 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( mandel_delta<2, 0, 1>(), mandel_delta<2, 0, 1>() ) - cst( 1.0 ) ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( canonical_tensor_basis_products_extract_expected_rows_columns_and_components )
{
    auto mesh = unitSquare();
    auto A = tensorExpr2D();
    auto v = vec( Px() + 2.0*Py(), cst( 1.0 ) - Py() );

    auto expectedDeltaTimesVector = vec( v( 1, 0 ), cst( 0.0 ) );
    auto expectedMatrixTimesDelta = mat<2, 2>( cst( 0.0 ), A( 0, 0 ),
                                               cst( 0.0 ), A( 1, 0 ) );
    auto expectedDeltaTimesMatrix = mat<2, 2>( A( 1, 0 ), A( 1, 1 ),
                                               cst( 0.0 ), cst( 0.0 ) );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, delta<2, 0, 1>()*v - expectedDeltaTimesVector ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, A*delta<2, 0, 1>() - expectedMatrixTimesDelta ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, delta<2, 0, 1>()*A - expectedDeltaTimesMatrix ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( runtime_tensor_basis_fallback_matches_component_access_and_dense_products )
{
    auto mesh = unitSquare();
    auto A = tensorExpr2D();
    auto v = vec( Px() + 2.0*Py(), cst( 1.0 ) - Py() );

    auto expectedDeltaTimesVector = vec( v( 1, 0 ), cst( 0.0 ) );
    auto expectedMatrixTimesDelta = mat<2, 2>( cst( 0.0 ), A( 0, 0 ),
                                               cst( 0.0 ), A( 1, 0 ) );
    auto expectedDeltaTimesMatrix = mat<2, 2>( A( 1, 0 ), A( 1, 1 ),
                                               cst( 0.0 ), cst( 0.0 ) );
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, delta<2>( 0, 1 ) ) - A( 0, 1 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, symm_delta<2>( 0, 1 ) ) -
                                                 ( A( 0, 1 ) + A( 1, 0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedAbsError( mesh, inner( A, mandel_delta<2>( 0, 1 ) ) -
                                                 cst( invSqrt2 )*( A( 0, 1 ) + A( 1, 0 ) ) ), 1e-12 );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, delta<2>( 0, 1 )*v - expectedDeltaTimesVector ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, A*delta<2>( 0, 1 ) - expectedMatrixTimesDelta ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, delta<2>( 0, 1 )*A - expectedDeltaTimesMatrix ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( symmetric_and_mandel_tensor_basis_products_build_on_structured_delta_products )
{
    auto mesh = unitSquare();
    auto A = tensorExpr2D();
    auto v = vec( Px() + 2.0*Py(), cst( 1.0 ) - Py() );
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    auto expectedSymmTimesVector = vec( v( 1, 0 ), v( 0, 0 ) );
    auto expectedMandelTimesVector = cst( invSqrt2 )*expectedSymmTimesVector;
    auto expectedMatrixTimesSymm = mat<2, 2>( A( 0, 1 ), A( 0, 0 ),
                                              A( 1, 1 ), A( 1, 0 ) );
    auto expectedSymmTimesMatrix = mat<2, 2>( A( 1, 0 ), A( 1, 1 ),
                                              A( 0, 0 ), A( 0, 1 ) );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, symm_delta<2, 0, 1>()*v - expectedSymmTimesVector ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, mandel_delta<2, 0, 1>()*v - expectedMandelTimesVector ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, A*symm_delta<2, 0, 1>() - expectedMatrixTimesSymm ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh, symm_delta<2, 0, 1>()*A - expectedSymmTimesMatrix ), 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
