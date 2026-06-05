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
 * \file test_sb9_common.hpp
 * \brief Shared mesh, field, and energy helpers for SB9 Feel++ VF tests.
 */
#ifndef FEELPP_TEST_SB9_COMMON_HPP
#define FEELPP_TEST_SB9_COMMON_HPP 1

#include <cmath>
#include <iomanip>
#include <sstream>

#include <Eigen/Core>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelvf/vf.hpp>

/**
 * \namespace Feel::Tests::SB9
 * \brief Reusable fixtures and helper expressions for SB9 shell tests.
 */
namespace Feel::Tests::SB9
{
/// Q1 hexahedral mesh type used by the local SB9 shell patches.
using mesh_type = Mesh<Hypercube<3>>;
/// Shared pointer to the SB9 test mesh type.
using mesh_ptrtype = std::shared_ptr<mesh_type>;
/// Three-component vector used to define patch geometry.
using test_vector_type = Eigen::Matrix<double, 3, 1>;
/// Three-by-three matrix used for patch rotations.
using matrix_type = Eigen::Matrix<double, 3, 3>;

/// Default tolerance for zero-energy rigid mode checks.
inline constexpr double g_tol = 1.e-9;

/**
 * \brief Create a single-element shell patch mesh from geometric vectors.
 *
 * The generated mesh starts from the quadrilateral mid-surface defined by
 * \p origin, \p tangent1, and \p tangent2, then extrudes it by
 * \p thicknessVector to produce one recombined hexahedral shell cell.
 *
 * \param caseName Name used for the generated `.geo` file and test directory.
 * \param origin Center point of the lower-left mid-surface corner.
 * \param tangent1 First in-plane edge vector.
 * \param tangent2 Second in-plane edge vector.
 * \param thicknessVector Through-thickness extrusion vector.
 * \return Generated one-cell shell mesh.
 */
inline mesh_ptrtype
createShellPatch( std::string const& caseName,
                  test_vector_type const& origin,
                  test_vector_type const& tangent1,
                  test_vector_type const& tangent2,
                  test_vector_type const& thicknessVector )
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
            << "Physical Volume(\"Shell\") = {out[1]};\n";

    Environment::changeRepository( _directory=boost::format( "testsuite/feelvf/%1%/%2%/" )
                                   % Environment::about().appName()
                                   % caseName );

    return createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename=caseName + ".geo",
                                      _desc=geoDesc.str(),
                                      _dim=3,
                                      _order=1,
                                      _h=1.0 ),
                           _force_rebuild=true );
}

/**
 * \brief Compute the quadratic energy associated with a closed bilinear form.
 *
 * The helper closes the form and evaluates `testElement^T A trialElement`
 * through the assembled matrix.
 *
 * \tparam FormType Feel++ bilinear form type.
 * \tparam TestElement Test-space element type.
 * \tparam TrialElement Trial-space element type.
 * \param form Bilinear form to close and evaluate.
 * \param testElement Test element used on the left side.
 * \param trialElement Trial element used on the right side.
 * \return Matrix energy for the two supplied elements.
 */
template <typename FormType, typename TestElement, typename TrialElement>
double
formEnergy( FormType& form,
            TestElement const& testElement,
            TrialElement const& trialElement )
{
    form.close();
    return form.matrixPtr()->energy( testElement, trialElement );
}

/**
 * \brief Return a constant rigid translation field used by SB9 tests.
 * \return Feel++ vector expression with three constant components.
 */
inline auto
rigidTranslationField()
{
    return vf::vec( 1.3, -0.7, 0.4 );
}

/**
 * \brief Return a small rigid rotation field around the origin.
 *
 * The expression has the form `omega x x` and should produce zero strain
 * energy for shell kinematic operators that preserve rigid body modes.
 *
 * \return Feel++ vector expression for a rigid rotation.
 */
inline auto
rigidRotationField()
{
    constexpr double wx = 0.31;
    constexpr double wy = -0.27;
    constexpr double wz = 0.19;

    return vf::vec( wy*vf::Pz() - wz*vf::Py(),
                    wz*vf::Px() - wx*vf::Pz(),
                    wx*vf::Py() - wy*vf::Px() );
}

/**
 * \brief Create a flat rectangular shell patch.
 *
 * \param caseName Name used for generated mesh artifacts.
 * \return One-cell shell mesh aligned with the coordinate axes.
 */
inline mesh_ptrtype
createFlatShellPatch( std::string const& caseName )
{
    return createShellPatch( caseName,
                             test_vector_type( 0.0, 0.0, 0.0 ),
                             test_vector_type( 2.0, 0.0, 0.0 ),
                             test_vector_type( 0.0, 1.0, 0.0 ),
                             test_vector_type( 0.0, 0.0, 0.2 ) );
}

/**
 * \brief Create a shell patch rotated in the x-z plane.
 *
 * \param caseName Name used for generated mesh artifacts.
 * \return One-cell shell mesh with rotated in-plane and thickness directions.
 */
inline mesh_ptrtype
createRotatedShellPatch( std::string const& caseName )
{
    double const angle = std::acos( -1.0 ) / 6.0;
    matrix_type rotation = matrix_type::Identity();
    rotation( 0, 0 ) = std::cos( angle );
    rotation( 0, 2 ) = std::sin( angle );
    rotation( 2, 0 ) = -std::sin( angle );
    rotation( 2, 2 ) = std::cos( angle );

    return createShellPatch( caseName,
                             test_vector_type( 0.0, 0.0, 0.0 ),
                             2.0 * rotation.col( 0 ),
                             rotation.col( 1 ),
                             0.2 * rotation.col( 2 ) );
}

/**
 * \brief Create an axis-aligned unit shell patch.
 *
 * The patch has unit in-plane lengths and unit thickness, which makes simple
 * extension energy checks easy to compare against closed-form values.
 *
 * \param caseName Name used for generated mesh artifacts.
 * \return One-cell axis-aligned unit shell mesh.
 */
inline mesh_ptrtype
createAxisAlignedUnitPatch( std::string const& caseName )
{
    return createShellPatch( caseName,
                             test_vector_type( 0.0, 0.0, 0.0 ),
                             test_vector_type( 1.0, 0.0, 0.0 ),
                             test_vector_type( 0.0, 1.0, 0.0 ),
                             test_vector_type( 0.0, 0.0, 1.0 ) );
}
} // namespace Feel::Tests::SB9

#endif
