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

#define BOOST_TEST_MODULE test_shellgeometric
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelvf/vf.hpp>

#include <Eigen/Core>

#include <cmath>
#include <iomanip>
#include <sstream>

using namespace Feel;

namespace
{
using mesh_type = Mesh<Hypercube<3>>;
using mesh_ptrtype = std::shared_ptr<mesh_type>;
using shell_test_vector_type = Eigen::Matrix<double, 3, 1>;
using shell_test_matrix_type = Eigen::Matrix<double, 3, 3>;

constexpr double g_tol = 1.e-10;

struct ShellPatchData
{
    mesh_ptrtype mesh;
    node_type center;
    node_type top;
    node_type bottom;
};

node_type
makeNode( shell_test_vector_type const& x )
{
    node_type n( 3 );
    n( 0 ) = x( 0 );
    n( 1 ) = x( 1 );
    n( 2 ) = x( 2 );
    return n;
}

auto
vectorExpr( shell_test_vector_type const& v )
{
    return vec( cst( v( 0 ) ), cst( v( 1 ) ), cst( v( 2 ) ) );
}

auto
matrixExpr( shell_test_matrix_type const& m )
{
    return mat<3, 3>( cst( m( 0, 0 ) ), cst( m( 0, 1 ) ), cst( m( 0, 2 ) ),
                      cst( m( 1, 0 ) ), cst( m( 1, 1 ) ), cst( m( 1, 2 ) ),
                      cst( m( 2, 0 ) ), cst( m( 2, 1 ) ), cst( m( 2, 2 ) ) );
}

ShellPatchData
createShellPatch( std::string const& caseName,
                  shell_test_vector_type const& origin,
                  shell_test_vector_type const& tangent1,
                  shell_test_vector_type const& tangent2,
                  shell_test_vector_type const& thicknessVector )
{
    auto const p1 = origin - 0.5 * thicknessVector;
    auto const p2 = p1 + tangent1;
    auto const p3 = p2 + tangent2;
    auto const p4 = p1 + tangent2;

    std::ostringstream geoDesc;
    geoDesc << std::setprecision( 16 );
    geoDesc << "Mesh.RecombineAll = 1;\n"
            << "Point(1) = {" << p1( 0 ) << ", " << p1( 1 ) << ", " << p1( 2 ) << ", 1};\n"
            << "Point(2) = {" << p2( 0 ) << ", " << p2( 1 ) << ", " << p2( 2 ) << ", 1};\n"
            << "Point(3) = {" << p3( 0 ) << ", " << p3( 1 ) << ", " << p3( 2 ) << ", 1};\n"
            << "Point(4) = {" << p4( 0 ) << ", " << p4( 1 ) << ", " << p4( 2 ) << ", 1};\n"
            << "Line(1) = {1, 2};\n"
            << "Line(2) = {2, 3};\n"
            << "Line(3) = {3, 4};\n"
            << "Line(4) = {4, 1};\n"
            << "Line Loop(1) = {1, 2, 3, 4};\n"
            << "Plane Surface(1) = {1};\n"
            << "Transfinite Line {1, 3} = 2;\n"
            << "Transfinite Line {2, 4} = 2;\n"
            << "Transfinite Surface {1} = {1, 2, 3, 4};\n"
            << "Recombine Surface {1};\n"
            << "out[] = Extrude {" << thicknessVector( 0 ) << ", "
            << thicknessVector( 1 ) << ", " << thicknessVector( 2 ) << "} {\n"
            << "  Surface{1};\n"
            << "  Layers{1};\n"
            << "  Recombine;\n"
            << "};\n"
            << "Physical Surface(\"Bottom\") = {1};\n"
            << "Physical Surface(\"Pressure\") = {out[0]};\n"
            << "Physical Surface(\"Clamp\") = {out[2], out[3], out[4], out[5]};\n"
            << "Physical Volume(\"Shell\") = {out[1]};\n";

    Feel::Environment::changeRepository( _directory=boost::format( "testsuite/feelvf/%1%/%2%/" )
                                         % Feel::Environment::about().appName()
                                         % caseName );

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=geo( _filename=caseName + ".geo",
                                           _desc=geoDesc.str(),
                                           _dim=3,
                                           _order=1,
                                           _h=1.0 ),
                                _force_rebuild=true );

    auto const center = origin + 0.5 * tangent1 + 0.5 * tangent2;
    return ShellPatchData{ mesh,
                           makeNode( center ),
                           makeNode( center + 0.5 * thicknessVector ),
                           makeNode( center - 0.5 * thicknessVector ) };
}

template <typename ExprT>
double
scalarL2Error( mesh_ptrtype const& mesh, ExprT const& expr )
{
    return normL2( _range=elements( mesh ), _expr=expr );
}

template <typename ExprT>
double
fieldL2Magnitude( mesh_ptrtype const& mesh, ExprT const& expr )
{
    return std::sqrt( integrate( _range=elements( mesh ), _expr=inner( expr ) ).evaluate()( 0, 0 ) );
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( shellgeometric_suite )

BOOST_AUTO_TEST_CASE( flat_shell_patch )
{
    shell_test_vector_type const origin( 0., 0., 0. );
    shell_test_vector_type const tangent1( 2., 0., 0. );
    shell_test_vector_type const tangent2( 0., 1., 0. );
    shell_test_vector_type const thicknessVector( 0., 0., 0.2 );
    shell_test_matrix_type const expectedFrame = shell_test_matrix_type::Identity();
    shell_test_matrix_type const expectedCovariant =
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 0.5, 0.0,
            0.0, 0.0, 0.1 ).finished();
    shell_test_matrix_type const expectedContravariant =
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 2.0, 0.0,
            0.0, 0.0, 10.0 ).finished();
    shell_test_matrix_type const expectedMetric =
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 0.25, 0.0,
            0.0, 0.0, 0.01 ).finished();

    auto patch = createShellPatch( "flat_patch", origin, tangent1, tangent2, thicknessVector );

    auto Zh = Pch<1>( patch.mesh );
    auto Hh = Pdh<0>( patch.mesh );
    auto Nh = Pdhv<0>( patch.mesh );
    auto Rh = Pdhm<0>( patch.mesh );

    auto zetaField = vf::project( _space=Zh, _range=elements( patch.mesh ), _expr=zeta() );
    auto areaField = vf::project( _space=Hh, _range=elements( patch.mesh ), _expr=shellArea0() );
    auto thicknessField = vf::project( _space=Hh, _range=elements( patch.mesh ), _expr=shellThickness() );
    auto normalField = vf::project( _space=Nh, _range=elements( patch.mesh ), _expr=shellNormal() );
    auto frameField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellFrame() );
    auto covariantField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellCovariantBasis0() );
    auto contravariantField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellContravariantBasis0() );
    auto metricField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellMetric0() );
    auto jacobianField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellJacobian0() );
    auto invJ0Mat = mat<3,3>( shellInvJ0_00(), shellInvJ0_01(), shellInvJ0_02(),
                              shellInvJ0_10(), shellInvJ0_11(), shellInvJ0_12(),
                              shellInvJ0_20(), shellInvJ0_21(), shellInvJ0_22() );
    auto invJ0Field = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=invJ0Mat );

    auto ctx = Zh->context();
    ctx.add( patch.center );
    ctx.add( patch.top );
    ctx.add( patch.bottom );
    auto zetaValues = zetaField.evaluate( ctx );

    BOOST_CHECK_SMALL( zetaValues( 0 ), g_tol );
    BOOST_CHECK_SMALL( zetaValues( 1 ) - 1.0, g_tol );
    BOOST_CHECK_SMALL( zetaValues( 2 ) + 1.0, g_tol );

    auto const expectedNormal = shell_test_vector_type( 0., 0., 1. );
    auto areaError = scalarL2Error( patch.mesh, idv( areaField ) - cst( 2.0 ) );
    auto thicknessError = scalarL2Error( patch.mesh, idv( thicknessField ) - cst( 0.2 ) );
    auto normalError = fieldL2Magnitude( patch.mesh, idv( normalField ) - vectorExpr( expectedNormal ) );
    auto frameError = fieldL2Magnitude( patch.mesh, idv( frameField ) - matrixExpr( expectedFrame ) );
    auto covariantError = fieldL2Magnitude( patch.mesh, idv( covariantField ) - matrixExpr( expectedCovariant ) );
    auto contravariantError = fieldL2Magnitude( patch.mesh, idv( contravariantField ) - matrixExpr( expectedContravariant ) );
    auto metricError = fieldL2Magnitude( patch.mesh, idv( metricField ) - matrixExpr( expectedMetric ) );
    auto jacobianError = fieldL2Magnitude( patch.mesh, idv( jacobianField ) - matrixExpr( expectedCovariant ) );
    auto invJ0Error = fieldL2Magnitude( patch.mesh, idv( invJ0Field ) - matrixExpr( expectedContravariant ) );
    auto frameOrthoError = fieldL2Magnitude( patch.mesh, trans( idv( frameField ) ) * idv( frameField ) - eye<3, 3>() );
    auto dualityError = fieldL2Magnitude( patch.mesh,
                                          trans( idv( contravariantField ) ) * idv( covariantField ) - eye<3, 3>() );

    BOOST_CHECK_SMALL( areaError, g_tol );
    BOOST_CHECK_SMALL( thicknessError, g_tol );
    BOOST_CHECK_SMALL( normalError, g_tol );
    BOOST_CHECK_SMALL( frameError, g_tol );
    BOOST_CHECK_SMALL( covariantError, g_tol );
    BOOST_CHECK_SMALL( contravariantError, g_tol );
    BOOST_CHECK_SMALL( metricError, g_tol );
    BOOST_CHECK_SMALL( jacobianError, g_tol );
    BOOST_CHECK_SMALL( invJ0Error, g_tol );
    BOOST_CHECK_SMALL( frameOrthoError, g_tol );
    BOOST_CHECK_SMALL( dualityError, g_tol );
}

BOOST_AUTO_TEST_CASE( rotated_shell_patch )
{
    double const angle = std::acos( -1.0 ) / 6.0;
    shell_test_matrix_type rotation = shell_test_matrix_type::Identity();
    rotation( 0, 0 ) = std::cos( angle );
    rotation( 0, 2 ) = std::sin( angle );
    rotation( 2, 0 ) = -std::sin( angle );
    rotation( 2, 2 ) = std::cos( angle );

    shell_test_vector_type const tangent1 = 2.0 * rotation.col( 0 );
    shell_test_vector_type const tangent2 = 1.0 * rotation.col( 1 );
    shell_test_vector_type const thicknessVector = 0.2 * rotation.col( 2 );
    shell_test_matrix_type expectedFrame = shell_test_matrix_type::Identity();
    expectedFrame.col( 0 ) = rotation.col( 0 );
    expectedFrame.col( 1 ) = rotation.col( 1 );
    expectedFrame.col( 2 ) = rotation.col( 2 );
    shell_test_matrix_type const expectedJacobian =
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 0.5, 0.0,
            0.0, 0.0, 0.1 ).finished();
    shell_test_matrix_type const expectedInvJ0 =
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 2.0, 0.0,
            0.0, 0.0, 10.0 ).finished();
    shell_test_matrix_type const expectedCovariant = rotation * expectedJacobian;
    shell_test_matrix_type const expectedContravariant = rotation *
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 2.0, 0.0,
            0.0, 0.0, 10.0 ).finished();
    shell_test_matrix_type const expectedMetric =
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 0.25, 0.0,
            0.0, 0.0, 0.01 ).finished();

    auto patch = createShellPatch( "rotated_patch", shell_test_vector_type( 0., 0., 0. ), tangent1, tangent2, thicknessVector );

    auto Hh = Pdh<0>( patch.mesh );
    auto Nh = Pdhv<0>( patch.mesh );
    auto Rh = Pdhm<0>( patch.mesh );

    auto areaField = vf::project( _space=Hh, _range=elements( patch.mesh ), _expr=shellArea0() );
    auto thicknessField = vf::project( _space=Hh, _range=elements( patch.mesh ), _expr=shellThickness() );
    auto normalField = vf::project( _space=Nh, _range=elements( patch.mesh ), _expr=shellNormal() );
    auto frameField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellFrame() );
    auto covariantField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellCovariantBasis0() );
    auto contravariantField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellContravariantBasis0() );
    auto metricField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellMetric0() );
    auto jacobianField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellJacobian0() );
    auto invJ0Mat = mat<3,3>( shellInvJ0_00(), shellInvJ0_01(), shellInvJ0_02(),
                              shellInvJ0_10(), shellInvJ0_11(), shellInvJ0_12(),
                              shellInvJ0_20(), shellInvJ0_21(), shellInvJ0_22() );
    auto invJ0Field = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=invJ0Mat );

    shell_test_vector_type expectedNormal = rotation.col( 2 );
    auto areaError = scalarL2Error( patch.mesh, idv( areaField ) - cst( 2.0 ) );
    auto thicknessError = scalarL2Error( patch.mesh, idv( thicknessField ) - cst( 0.2 ) );
    auto normalError = fieldL2Magnitude( patch.mesh, idv( normalField ) - vectorExpr( expectedNormal ) );
    auto frameError = fieldL2Magnitude( patch.mesh, idv( frameField ) - matrixExpr( expectedFrame ) );
    auto covariantError = fieldL2Magnitude( patch.mesh, idv( covariantField ) - matrixExpr( expectedCovariant ) );
    auto contravariantError = fieldL2Magnitude( patch.mesh, idv( contravariantField ) - matrixExpr( expectedContravariant ) );
    auto metricError = fieldL2Magnitude( patch.mesh, idv( metricField ) - matrixExpr( expectedMetric ) );
    auto jacobianError = fieldL2Magnitude( patch.mesh, idv( jacobianField ) - matrixExpr( expectedJacobian ) );
    auto invJ0Error = fieldL2Magnitude( patch.mesh, idv( invJ0Field ) - matrixExpr( expectedInvJ0 ) );
    auto frameOrthoError = fieldL2Magnitude( patch.mesh, trans( idv( frameField ) ) * idv( frameField ) - eye<3, 3>() );
    auto dualityError = fieldL2Magnitude( patch.mesh,
                                          trans( idv( contravariantField ) ) * idv( covariantField ) - eye<3, 3>() );

    BOOST_CHECK_SMALL( areaError, g_tol );
    BOOST_CHECK_SMALL( thicknessError, g_tol );
    BOOST_CHECK_SMALL( normalError, g_tol );
    BOOST_CHECK_SMALL( frameError, g_tol );
    BOOST_CHECK_SMALL( covariantError, g_tol );
    BOOST_CHECK_SMALL( contravariantError, g_tol );
    BOOST_CHECK_SMALL( metricError, g_tol );
    BOOST_CHECK_SMALL( jacobianError, g_tol );
    BOOST_CHECK_SMALL( invJ0Error, g_tol );
    BOOST_CHECK_SMALL( frameOrthoError, g_tol );
    BOOST_CHECK_SMALL( dualityError, g_tol );
}

BOOST_AUTO_TEST_SUITE_END()
