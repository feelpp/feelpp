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

#define BOOST_TEST_MODULE voigt testsuite
#include <feel/feelcore/testsuite.hpp>

#include <cmath>
#include <numbers>

#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

template <typename MeshType, typename ExprType>
double
integratedSquaredNorm( std::shared_ptr<MeshType> const& mesh, ExprType const& expr )
{
    return integrate( _range=elements( mesh ), _expr=inner( expr, expr ) ).evaluate()( 0, 0 );
}

template <typename MeshType, typename ExprType>
double
integratedValue( std::shared_ptr<MeshType> const& mesh, ExprType const& expr )
{
    return integrate( _range=elements( mesh ), _expr=expr ).evaluate()( 0, 0 );
}

template <int Dim>
auto unitMesh()
{
    if constexpr ( Dim == 2 )
        return unitSquare();
    else
        return unitCube();
}

auto symmetricMatrix2D()
{
    return mat<2, 2>( cst( 1.0 ), cst( 2.0 ),
                      cst( 2.0 ), cst( 3.0 ) );
}

auto symmetricMatrix3D()
{
    return mat<3, 3>( cst( 1.0 ), cst( 2.0 ), cst( 3.0 ),
                      cst( 2.0 ), cst( 4.0 ), cst( 5.0 ),
                      cst( 3.0 ), cst( 5.0 ), cst( 6.0 ) );
}

} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( voigt_suite )

BOOST_AUTO_TEST_CASE( voigt_uses_canonical_symmetric_ordering )
{
    auto mesh2 = unitMesh<2>();
    auto mesh3 = unitMesh<3>();
    constexpr auto sqrt2 = std::numbers::sqrt2_v<double>;

    auto expected2DVoigt = vec( cst( 1.0 ), cst( 2.0 ), cst( 3.0 ) );
    auto expected2DMandel = vec( cst( 1.0 ), cst( sqrt2*2.0 ), cst( 3.0 ) );
    auto expected3DVoigt = vec( cst( 1.0 ), cst( 2.0 ), cst( 3.0 ),
                                cst( 4.0 ), cst( 5.0 ), cst( 6.0 ) );
    auto expected3DMandel = vec( cst( 1.0 ), cst( sqrt2*2.0 ), cst( sqrt2*3.0 ),
                                 cst( 4.0 ), cst( sqrt2*5.0 ), cst( 6.0 ) );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, voigt( symmetricMatrix2D() ) - expected2DVoigt ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, mandel( symmetricMatrix2D() ) - expected2DMandel ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, voigt( symmetricMatrix3D() ) - expected3DVoigt ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, mandel( symmetricMatrix3D() ) - expected3DMandel ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( voigt_and_mandel_roundtrip_to_symmetric_matrices )
{
    auto mesh2 = unitMesh<2>();
    auto mesh3 = unitMesh<3>();
    auto S2 = symmetricMatrix2D();
    auto S3 = symmetricMatrix3D();

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, unvoigt( voigt( S2 ) ) - S2 ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, unmandel( mandel( S2 ) ) - S2 ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, unvoigt( trans( voigt( S2 ) ) ) - S2 ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, unmandel( trans( mandel( S2 ) ) ) - S2 ), 1e-12 );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, unvoigt( voigt( S3 ) ) - S3 ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, unmandel( mandel( S3 ) ) - S3 ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, unvoigt( trans( voigt( S3 ) ) ) - S3 ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, unmandel( trans( mandel( S3 ) ) ) - S3 ), 1e-12 );

    BOOST_CHECK_SMALL( integratedValue( mesh2, inner( mandel( S2 ), mandel( S2 ) ) ) - integratedValue( mesh2, inner( S2, S2 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedValue( mesh3, inner( mandel( S3 ), mandel( S3 ) ) ) - integratedValue( mesh3, inner( S3, S3 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedValue( mesh2, voigt_inner( voigt( S2 ), voigt( S2 ) ) ) - integratedValue( mesh2, inner( S2, S2 ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedValue( mesh3, voigt_inner( voigt( S3 ), voigt( S3 ) ) ) - integratedValue( mesh3, inner( S3, S3 ) ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( unvoigt_and_unmandel_match_manual_reconstruction )
{
    auto mesh2 = unitMesh<2>();
    auto mesh3 = unitMesh<3>();
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    auto v3 = vec( cst( 1.0 ), cst( 2.0 ), cst( 3.0 ) );
    auto m3 = vec( cst( 1.0 ), cst( 2.0*std::numbers::sqrt2_v<double> ), cst( 3.0 ) );
    auto expected2D = symmetricMatrix2D();

    auto v6 = vec( cst( 1.0 ), cst( 2.0 ), cst( 3.0 ),
                   cst( 4.0 ), cst( 5.0 ), cst( 6.0 ) );
    auto m6 = vec( cst( 1.0 ), cst( 2.0*std::numbers::sqrt2_v<double> ), cst( 3.0*std::numbers::sqrt2_v<double> ),
                   cst( 4.0 ), cst( 5.0*std::numbers::sqrt2_v<double> ), cst( 6.0 ) );
    auto expected3D = symmetricMatrix3D();

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, unvoigt( v3 ) - expected2D ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, unmandel( m3 ) - expected2D ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, unvoigt( v6 ) - expected3D ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, unmandel( m6 ) - expected3D ), 1e-12 );

    auto expected2DMandelManual = mat<2, 2>( cst( 1.0 ), cst( 2.0*invSqrt2 ),
                                             cst( 2.0*invSqrt2 ), cst( 3.0 ) );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, unmandel( vec( cst( 1.0 ), cst( 2.0 ), cst( 3.0 ) ) ) - expected2DMandelManual ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( storage_vector_builders_match_tensor_conversions )
{
    auto mesh2 = unitMesh<2>();
    auto mesh3 = unitMesh<3>();

    auto expected2DVoigt = voigt( symmetricMatrix2D() );
    auto expected2DMandel = mandel( symmetricMatrix2D() );
    auto expected3DVoigt = voigt( symmetricMatrix3D() );
    auto expected3DMandel = mandel( symmetricMatrix3D() );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, voigt_vec<2>( cst( 1.0 ), cst( 3.0 ), cst( 2.0 ) ) - expected2DVoigt ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, mandel_vec<2>( cst( 1.0 ), cst( 3.0 ), cst( 2.0 ) ) - expected2DMandel ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, voigt_vec<3>( cst( 1.0 ), cst( 4.0 ), cst( 6.0 ),
                                                                    cst( 2.0 ), cst( 3.0 ), cst( 5.0 ) ) - expected3DVoigt ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, mandel_vec<3>( cst( 1.0 ), cst( 4.0 ), cst( 6.0 ),
                                                                     cst( 2.0 ), cst( 3.0 ), cst( 5.0 ) ) - expected3DMandel ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( storage_basis_vectors_match_expected_canonical_order )
{
    auto mesh2 = unitMesh<2>();
    auto mesh3 = unitMesh<3>();
    constexpr auto sqrt2 = std::numbers::sqrt2_v<double>;

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, voigt_basis<2, 0, 1>() -
                                                     vec( cst( 0.0 ), cst( 1.0 ), cst( 0.0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, mandel_basis<2, 0, 1>() -
                                                      vec( cst( 0.0 ), cst( sqrt2 ), cst( 0.0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, voigt_basis<3, 0, 2>() -
                                                     vec( cst( 0.0 ), cst( 0.0 ), cst( 1.0 ),
                                                          cst( 0.0 ), cst( 0.0 ), cst( 0.0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, mandel_basis<3, 0, 2>() -
                                                      vec( cst( 0.0 ), cst( 0.0 ), cst( sqrt2 ),
                                                           cst( 0.0 ), cst( 0.0 ), cst( 0.0 ) ) ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( storage_component_builders_embed_scalar_components_in_canonical_order )
{
    auto mesh2 = unitMesh<2>();
    auto mesh3 = unitMesh<3>();
    constexpr auto sqrt2 = std::numbers::sqrt2_v<double>;

    auto scalar2 = Px() + 2.0*Py();
    auto scalar3 = Px() - 3.0*Py() + 2.0*Pz();

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, voigt_component<2, 0, 1>( scalar2 ) -
                                                     vec( cst( 0.0 ), scalar2, cst( 0.0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh2, mandel_component<2, 0, 1>( scalar2 ) -
                                                      vec( cst( 0.0 ), cst( sqrt2 )*scalar2, cst( 0.0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, voigt_component<3, 0, 2>( scalar3 ) -
                                                     vec( cst( 0.0 ), cst( 0.0 ), scalar3,
                                                          cst( 0.0 ), cst( 0.0 ), cst( 0.0 ) ) ), 1e-12 );
    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3, mandel_component<3, 0, 2>( scalar3 ) -
                                                      vec( cst( 0.0 ), cst( 0.0 ), cst( sqrt2 )*scalar3,
                                                           cst( 0.0 ), cst( 0.0 ), cst( 0.0 ) ) ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( scale_symm_storage_preserves_canonical_storage_shape )
{
    auto mesh3 = unitMesh<3>();
    auto scalar = Px() - 2.0*Py() + 0.5*Pz();
    auto scale = 3.0 - Px() + Py();
    auto eps = mandel_component<3, 0, 2>( scalar );

    BOOST_CHECK_SMALL( integratedSquaredNorm( mesh3,
                                              scale_symm_storage( scale, eps ) -
                                              mandel_component<3, 0, 2>( scale*scalar ) ),
                       1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
