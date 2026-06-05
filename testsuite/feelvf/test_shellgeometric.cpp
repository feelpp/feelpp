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
#include <feel/feeldiscr/product.hpp>
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

using Hh_ptr_t =  Pdh_ptrtype<mesh_type,0>;
// using product_spaceHh_type = ProductSpaces<Hh_ptr_t,Hh_ptr_t,Hh_ptr_t,Hh_ptr_t,Hh_ptr_t,Hh_ptr_t,Hh_ptr_t,Hh_ptr_t>;
// using product_spaceHh_ptrtype = std::shared_ptr<product_spaceHh_type>;
using product_spaceHh_type = dyn_product_space_t<Hh_ptr_t>;
using product_spaceHh_ptrtype = dyn_product_space_ptr_t<Hh_ptr_t>;

using shell_test_vector_type = Eigen::Matrix<double, 3, 1>;
using shell_test_hallquist_type = Eigen::Matrix<double, 8, 1>;
using shell_test_matrix_type = Eigen::Matrix<double, 3, 3>;
using shell_test_matrix2_type = Eigen::Matrix<double, 2, 2>;

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
hallquistExpr( shell_test_hallquist_type const& v )
{
    return vec( cst( v( 0 ) ), cst( v( 1 ) ), cst( v( 2 ) ), cst( v( 3 ) ), cst( v( 4 ) ), cst( v( 5 ) ), cst( v( 6 ) ), cst( v( 7 ) ) );
}

auto
matrixExpr( shell_test_matrix_type const& m )
{
    return mat<3, 3>( cst( m( 0, 0 ) ), cst( m( 0, 1 ) ), cst( m( 0, 2 ) ),
                      cst( m( 1, 0 ) ), cst( m( 1, 1 ) ), cst( m( 1, 2 ) ),
                      cst( m( 2, 0 ) ), cst( m( 2, 1 ) ), cst( m( 2, 2 ) ) );
}

auto
matrix2Expr( shell_test_matrix2_type const& m )
{
    return mat<2, 2>( cst( m( 0, 0 ) ), cst( m( 0, 1 ) ),
                      cst( m( 1, 0 ) ), cst( m( 1, 1 ) ) );
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
    shell_test_hallquist_type const expectedBx =
        ( shell_test_hallquist_type() <<
            -0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125 ).finished();
    shell_test_hallquist_type const expectedBy =
        ( shell_test_hallquist_type() <<
            -0.25, -0.25, 0.25, 0.25, -0.25, -0.25, 0.25, 0.25 ).finished();
    shell_test_hallquist_type const expectedBz =
        ( shell_test_hallquist_type() <<
            -1.25, -1.25, -1.25, -1.25, 1.25, 1.25, 1.25, 1.25 ).finished();
    
    // soit def séparemment puis vec( ..._.T ? ) pour construire expectedVgamma ou def complètement séparemment ou permet d'appeller Vgamma(0) etc pour comparer un par un à ces vecteurs
    shell_test_hallquist_type const expectedVgamma1 =
        ( shell_test_hallquist_type() <<
            0.125, 0.125, -0.125, -0.125, -0.125, -0.125, 0.125, 0.125 ).finished();
    shell_test_hallquist_type const expectedVgamma2 =
        ( shell_test_hallquist_type() <<
            0.125, -0.125, -0.125, 0.125, -0.125, 0.125, 0.125, -0.125 ).finished();
    shell_test_hallquist_type const expectedVgamma3 =
        ( shell_test_hallquist_type() <<
            0.125, -0.125, 0.125, -0.125, 0.125, -0.125, 0.125, -0.125 ).finished();
    shell_test_hallquist_type const expectedVgamma4 =
        ( shell_test_hallquist_type() <<
            -0.125, 0.125, -0.125, 0.125, 0.125, -0.125, 0.125, -0.125 ).finished();

    shell_test_matrix2_type const expectedJai =
        ( shell_test_matrix2_type() <<
            1.0, 0.0,
            0.0, 0.5 ).finished();
    shell_test_matrix_type const expectedInvJai =
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 2.0, 0.0,
            0.0, 0.0, 10.0 ).finished();

    auto patch = createShellPatch( "flat_patch", origin, tangent1, tangent2, thicknessVector );

    auto Zh = Pch<1>( patch.mesh );
    auto Hh = Pdh<0>( patch.mesh );
    auto Nh = Pdhv<0>( patch.mesh );
    auto Rh = Pdhm<0>( patch.mesh );

    // product_spaceHh_ptrtype psh = std::shared_ptr< product( Hh, Hh, Hh, Hh, Hh, Hh, Hh, Hh ) >;    // euh moue...  -> mm espace pour vGamma vec(vGamma(1) .. 4) ? ou def un espace à part de dim (4,8) ?
    // product_spaceHh_ptrtype psh = std::make_shared<product_spaceHh_type>( Hh, Hh, Hh, Hh, Hh, Hh, Hh, Hh );
    // product_spaceHh_type psh = dynProductPtr( 8, Hh );
    auto psh = dynProductPtr( 8, Hh );
    // product_spaceHh_ptrtype psh = std::make_shared<product_spaceHh_type>( 8, Hh );

    // auto psm = product( ? ); pour obtenir un espace matriciel discret 2*2



    auto zetaField = vf::project( _space=Zh, _range=elements( patch.mesh ), _expr=zeta() );
    auto areaField = vf::project( _space=Hh, _range=elements( patch.mesh ), _expr=shellArea0() );
    auto thicknessField = vf::project( _space=Hh, _range=elements( patch.mesh ), _expr=shellThickness() );
    auto normalField = vf::project( _space=Nh, _range=elements( patch.mesh ), _expr=shellNormal() );
    auto frameField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellFrame() );
    auto covariantField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellCovariantBasis0() );
    auto contravariantField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellContravariantBasis0() );
    auto metricField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellMetric0() );
    auto jacobianField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellJacobian0() );
    auto invJaField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellInvJa() );
    auto invJbField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellInvJb() );
    auto invJcField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellInvJc() );
    auto invJdField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellInvJd() );

    // auto JaField = vf::project( _space=psm, _range=elements( patch.mesh ), _expr=shellJa() );
    // auto JbField = vf::project( _space=psm, _range=elements( patch.mesh ), _expr=shellJb() );
    // auto JcField = vf::project( _space=psm, _range=elements( patch.mesh ), _expr=shellJc() );
    // auto JdField = vf::project( _space=psm, _range=elements( patch.mesh ), _expr=shellJd() );



    // auto bxField = vf::project( _space=psh, _range=elements( patch.mesh ), _expr=shellBx() );
    // for... bxField(i).on(_range=elements( patch.mesh ), _expr=shellBx(i)); // faut alors permettre d'appeller composante par composante dans test_shellgeometric.cpp
    // ex.add(bx(i))   pour exporter   -> composante par comp

    // auto byField = vf::project( _space=psh, _range=elements( patch.mesh ), _expr=shellBy() );
    // auto bzField = vf::project( _space=psh, _range=elements( patch.mesh ), _expr=shellBz() );

    auto ctx = Zh->context();
    ctx.add( patch.center );
    ctx.add( patch.top );
    ctx.add( patch.bottom );
    auto zetaValues = zetaField.evaluate( ctx );        // commun pour tous les éléments - pareil pour hallquist

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
    auto frameOrthoError = fieldL2Magnitude( patch.mesh, trans( idv( frameField ) ) * idv( frameField ) - eye<3, 3>() );
    auto dualityError = fieldL2Magnitude( patch.mesh, trans( idv( contravariantField ) ) * idv( covariantField ) - eye<3, 3>() );
    auto invJaError = fieldL2Magnitude( patch.mesh, idv( invJaField ) - matrixExpr( expectedInvJai ) );
    auto invJbError = fieldL2Magnitude( patch.mesh, idv( invJbField ) - matrixExpr( expectedInvJai ) );
    auto invJcError = fieldL2Magnitude( patch.mesh, idv( invJcField ) - matrixExpr( expectedInvJai ) );
    auto invJdError = fieldL2Magnitude( patch.mesh, idv( invJdField ) - matrixExpr( expectedInvJai ) );

    // auto JaError = fieldL2Magnitude( patch.mesh, idv( JaField ) - matrixExpr( expectedJai ) );
    // auto JbError = fieldL2Magnitude( patch.mesh, idv( JbField ) - matrixExpr( expectedJai ) );
    // auto JcError = fieldL2Magnitude( patch.mesh, idv( JcField ) - matrixExpr( expectedJai ) );
    // auto JdError = fieldL2Magnitude( patch.mesh, idv( JdField ) - matrixExpr( expectedJai ) );

    // auto bxError = fieldL2Magnitude( patch.mesh, shellBx() - hallquistExpr( expectedBx ) );
    // auto bxError = scalarL2Error( patch.mesh, idv( shellBx() ) - hallquistExpr( expectedBx ) );  
    // auto bxError = fieldL2Magnitude( patch.mesh, idv( bxField ) - hallquistExpr( expectedBx ) );
    // auto byError = fieldL2Magnitude( patch.mesh, idv( byField ) - hallquistExpr( expectedBy ) );
    // auto bzError = fieldL2Magnitude( patch.mesh, idv( bzField ) - hallquistExpr( expectedBz ) );


    BOOST_CHECK_SMALL( areaError, g_tol );
    BOOST_CHECK_SMALL( thicknessError, g_tol );
    BOOST_CHECK_SMALL( normalError, g_tol );
    BOOST_CHECK_SMALL( frameError, g_tol );
    BOOST_CHECK_SMALL( covariantError, g_tol );
    BOOST_CHECK_SMALL( contravariantError, g_tol );
    BOOST_CHECK_SMALL( metricError, g_tol );
    BOOST_CHECK_SMALL( jacobianError, g_tol );
    BOOST_CHECK_SMALL( frameOrthoError, g_tol );
    BOOST_CHECK_SMALL( dualityError, g_tol );
    BOOST_CHECK_SMALL( invJaError, g_tol );
    BOOST_CHECK_SMALL( invJbError, g_tol );
    BOOST_CHECK_SMALL( invJcError, g_tol );
    BOOST_CHECK_SMALL( invJdError, g_tol );

    // BOOST_CHECK_SMALL( JaError, g_tol );
    // BOOST_CHECK_SMALL( JbError, g_tol );
    // BOOST_CHECK_SMALL( JcError, g_tol );
    // BOOST_CHECK_SMALL( JdError, g_tol );

    // BOOST_CHECK_SMALL( bxError, g_tol );
    // BOOST_CHECK_SMALL( byError, g_tol );
    // BOOST_CHECK_SMALL( bzError, g_tol );
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
    shell_test_matrix_type const expectedInvJai =
        ( shell_test_matrix_type() <<
            1.0, 0.0, 0.0,
            0.0, 2.0, 0.0,
            0.0, 0.0, 10.0 ).finished();

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
    auto invJaField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellInvJa() );
    auto invJbField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellInvJb() );
    auto invJcField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellInvJc() );
    auto invJdField = vf::project( _space=Rh, _range=elements( patch.mesh ), _expr=shellInvJd() );

    shell_test_vector_type expectedNormal = rotation.col( 2 );
    auto areaError = scalarL2Error( patch.mesh, idv( areaField ) - cst( 2.0 ) );
    auto thicknessError = scalarL2Error( patch.mesh, idv( thicknessField ) - cst( 0.2 ) );
    auto normalError = fieldL2Magnitude( patch.mesh, idv( normalField ) - vectorExpr( expectedNormal ) );
    auto frameError = fieldL2Magnitude( patch.mesh, idv( frameField ) - matrixExpr( expectedFrame ) );
    auto covariantError = fieldL2Magnitude( patch.mesh, idv( covariantField ) - matrixExpr( expectedCovariant ) );
    auto contravariantError = fieldL2Magnitude( patch.mesh, idv( contravariantField ) - matrixExpr( expectedContravariant ) );
    auto metricError = fieldL2Magnitude( patch.mesh, idv( metricField ) - matrixExpr( expectedMetric ) );
    auto jacobianError = fieldL2Magnitude( patch.mesh, idv( jacobianField ) - matrixExpr( expectedJacobian ) );
    auto frameOrthoError = fieldL2Magnitude( patch.mesh, trans( idv( frameField ) ) * idv( frameField ) - eye<3, 3>() );
    auto dualityError = fieldL2Magnitude( patch.mesh, trans( idv( contravariantField ) ) * idv( covariantField ) - eye<3, 3>() );
    auto invJaError = fieldL2Magnitude( patch.mesh, idv( invJaField ) - matrixExpr( expectedInvJai ) );
    auto invJbError = fieldL2Magnitude( patch.mesh, idv( invJbField ) - matrixExpr( expectedInvJai ) );
    auto invJcError = fieldL2Magnitude( patch.mesh, idv( invJcField ) - matrixExpr( expectedInvJai ) );
    auto invJdError = fieldL2Magnitude( patch.mesh, idv( invJdField ) - matrixExpr( expectedInvJai ) );



    BOOST_CHECK_SMALL( areaError, g_tol );
    BOOST_CHECK_SMALL( thicknessError, g_tol );
    BOOST_CHECK_SMALL( normalError, g_tol );
    BOOST_CHECK_SMALL( frameError, g_tol );
    BOOST_CHECK_SMALL( covariantError, g_tol );
    BOOST_CHECK_SMALL( contravariantError, g_tol );
    BOOST_CHECK_SMALL( metricError, g_tol );
    BOOST_CHECK_SMALL( jacobianError, g_tol );
    BOOST_CHECK_SMALL( frameOrthoError, g_tol );
    BOOST_CHECK_SMALL( dualityError, g_tol );
    BOOST_CHECK_SMALL( invJaError, g_tol );
    BOOST_CHECK_SMALL( invJbError, g_tol );
    BOOST_CHECK_SMALL( invJcError, g_tol );
    BOOST_CHECK_SMALL( invJdError, g_tol );
}

BOOST_AUTO_TEST_SUITE_END()
