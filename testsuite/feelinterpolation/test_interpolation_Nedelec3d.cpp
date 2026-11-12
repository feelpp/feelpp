/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
    SPDX-FileContributor: Cécile Daversin <daversin@math.unistra.fr>

    SPDX-FileCopyrightText: 2007-2012 Université Joseph Fourier
    SPDX-FileCopyrightText: 2011-2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_interpolation_Nedelec3d.cpp
 * @brief Lowest-order Nedelec interpolation on affine tetrahedra.
 */
#define BOOST_TEST_MODULE test_interpolation_Nedelec3d
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelfilters/loadmesh.hpp>

#include "test_nedelec_interpolation.hpp"

using namespace Feel;

namespace
{
/** @brief Return command-line options for the 3D interpolation cases. */
po::options_description makeOptions()
{
    po::options_description options(
        "Nedelec tetrahedron interpolation options" );
    options.add_options()(
        "meshes-3d",
        po::value<std::vector<std::string>>(),
        "one-cell tetrahedron mesh paths" );
    return options;
}

/** @brief Return application metadata for the 3D regression. */
AboutData makeAbout()
{
    AboutData about(
        "test_interpolation_Nedelec3d",
        "test_interpolation_Nedelec3d", "0.2",
        "NED1 edge-moment interpolation on tetrahedra",
        AboutData::License_GPL,
        "Copyright (c) 2026 University of Strasbourg" );
    about.addAuthor(
        "Cécile Daversin", "developer",
        "daversin@math.unistra.fr", "" );
    about.addAuthor(
        "Christophe Prud'homme", "developer",
        "christophe.prudhomme@feelpp.org", "" );
    return about;
}
} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_CASE( tetrahedron_ned1_edge_moments )
{
    using mesh_type = Mesh<Simplex<3, 1>>;

    for ( auto const& meshPath :
          vsoption( _name = "meshes-3d" ) )
    {
        BOOST_TEST_CONTEXT( "mesh=" << meshPath )
        {
            auto mesh = loadMesh(
                _mesh = new mesh_type,
                _filename = meshPath + ".msh" );
            auto edgeMesh = createSubmesh(
                _mesh = mesh,
                _range = boundaryedges( mesh ) );
            auto space = Ned1h<0>( mesh );
            auto const expression = vec(
                cst( 1.0 ), cst( 1.0 ), cst( 1.0 ) );
            std::vector<std::string> const edgeMarkers{
                "zAxis", "yAxis", "yzAxis",
                "xyAxis", "xzAxis", "xAxis" };

            test::checkLowestOrderNedelecInterpolation(
                mesh, space, edgeMarkers, expression,
                [&]( std::string const& marker )
                {
                    return markedelements(
                        edgeMesh, marker );
                },
                []( auto const& range )
                {
                    return test::canonicalEdgeOrientation(
                        range );
                } );
        }
    }
}
