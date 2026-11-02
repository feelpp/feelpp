/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 SPDX-FileCopyrightText: 2019-2026 University of Strasbourg
 SPDX-License-Identifier: LGPL-3.0-or-later
 */

// give a name to the testsuite
#define BOOST_TEST_MODULE importer mesh testsuite

#include <feel/feelcore/testsuite.hpp>
#include <boost/test/data/test_case.hpp>
#include <array>
#include <type_traits>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <fstream>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_importer_mesh )

namespace bdata = boost::unit_test::data;
using namespace Feel;

template <typename Callable>
void runForDim( int dim, Callable&& c )
{
    switch ( dim )
    {
    case 1:
        c( std::integral_constant<int,1>{} );
        break;
    case 2:
        c( std::integral_constant<int,2>{} );
        break;
    case 3:
        c( std::integral_constant<int,3>{} );
        break;
    default:
        BOOST_FAIL( "Unsupported dimension " << dim );
    }
}

BOOST_DATA_TEST_CASE( test_scale, bdata::make( std::array<int,3>{ { 1,2,3 } } ), dim )
{
    BOOST_TEST_CONTEXT( "dim=" << dim )
    {
        runForDim( dim, []( auto d )
        {
    using mesh_type = Mesh<Simplex<decltype( d )::value>>;
    auto mesh = loadMesh( _mesh=new mesh_type);
    double scalingUsed = 10;
    auto meshWithScale = loadMesh( _mesh=new mesh_type,
                                   _scale=scalingUsed );

    double volume = mesh->measure();
    double volumeWithScale = meshWithScale->measure();
    BOOST_CHECK_CLOSE( volume*std::pow(scalingUsed, mesh_type::nDim), volumeWithScale, 1e-10 );

    // save mesh with hdf5 format
    double scalingUsed2 = 1./3.;
    std::string filenameHDF5 = (fs::current_path()/(boost::format("test_scale_mesh_%1%d.json")%mesh_type::nDim).str()).string();
    mesh->saveHDF5( filenameHDF5, scalingUsed2 );
    // load mesh with hdf5 format
    double scalingUsed3 = 5;
    auto meshHDF5 = loadMesh( _mesh=new mesh_type,_filename=filenameHDF5,_scale=scalingUsed3);
    double volumeHDF5 = meshHDF5->measure();
    BOOST_CHECK_CLOSE( volume*std::pow(scalingUsed2*scalingUsed3, mesh_type::nDim), volumeHDF5, 1e-10 );
        } );
    }
}

BOOST_AUTO_TEST_CASE( runtime_order_legacy_importers )
{
    using mesh_type = Mesh<Simplex<3, Dynamic>>;
    auto const prefix = ( fs::current_path() /
                          ( "runtime_import_" + std::to_string( Environment::rank() ) ) ).string();
    {
        std::ofstream( prefix + ".nodes" ) << "1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n";
        std::ofstream( prefix + ".cells" ) << "1 1 2 3 4\n";
        std::ofstream( prefix + ".faces" ) << "1 1 1 2 3\n";
        std::ofstream( prefix + ".arm" )
            << "<mesh><coordinates file=\"" << prefix << ".nodes\"/>"
            << "<element_set name=\"volume\" topology=\"tet4\" file=\"" << prefix << ".cells\"/>"
            << "<surface_set name=\"wall\" topology=\"tria3\" file=\"" << prefix << ".faces\"/></mesh>\n";
        std::ofstream( prefix + ".dat" )
            << ".NOE\nI 1 X 0 Y 0 Z 0\nI 2 X 1 Y 0 Z 0\n"
            << "I 3 X 0 Y 1 Z 0\nI 4 X 0 Y 0 Z 1\n"
            << ".MAI\nI 1 ATT 1 N 1 2 3 4\n.FIN\n";
    }
    for ( auto const& extension : { ".arm", ".dat" } )
    {
        auto mesh = loadMesh(
            _mesh = new mesh_type( RuntimeOrder{1}, "", Environment::worldCommSeqPtr() ),
            _worldcomm = Environment::worldCommSeqPtr(), _filename = prefix + extension );
        BOOST_CHECK_EQUAL( mesh->order(), 1 );
        BOOST_CHECK_EQUAL( mesh->numGlobalElements(), 1 );
        BOOST_CHECK_EQUAL( mesh->beginElement()->second.nPoints(), 4 );
        BOOST_CHECK_CLOSE( mesh->measure(), 1. / 6., 1e-10 );
    }
}

#if defined(FEELPP_HAS_HDF5)
BOOST_AUTO_TEST_CASE( runtime_order_hdf5_loading_and_repartition )
{
    using mesh_type = Mesh<Simplex<2, Dynamic>>;
    GeoTool::Node center( 0, 0 ), radius( 1, 0 );
    GeoTool::Circle circle( 0.3, "runtime_circle", center, radius );
    circle.setMarker( _type="line", _name="boundary", _markerAll=true );
    circle.setMarker( _type="surface", _name="domain", _markerAll=true );
    auto mesh = circle.createMesh( _mesh=new mesh_type( RuntimeOrder{2} ),
                                   _name="runtime_circle", _force_rebuild=true, _straighten=false );
    auto traceMesh = createSubmesh( _mesh=mesh, _range=boundaryfaces( mesh ) );
    for ( int order : { 0, 3 } )
    {
        auto traceSpace = Pdh<Dynamic>( traceMesh, RuntimeOrder{order} );
        auto value = traceSpace->element();
        auto exact = expr( order == 0 ? "2" : "1+x+y:x:y" );
        value.on( _range=boundaryfaces( mesh ), _expr=exact );
        // The closed trace has no boundary vertices: GEOMAP_OPT would use
        // straight segments for its interior cells instead of the curved map.
        BOOST_CHECK_SMALL( normL2( _range=elements( traceMesh ),
                                  _expr=idv( value ) - exact,
                                  _geomap=GeomapStrategyType::GEOMAP_HO ), 1e-10 );
        if ( order > 0 )
        {
            auto constants = Pdh<Dynamic>( traceMesh, RuntimeOrder{0} );
            auto constant = constants->element();
            // Identical C++ element types can still have incompatible runtime
            // coefficient layouts; neither arithmetic direction may index them.
            BOOST_CHECK_THROW( value += constant, std::invalid_argument );
            BOOST_CHECK_THROW( constant -= value, std::invalid_argument );
        }
    }
    auto const filename = ( fs::current_path() / "runtime_circle.json" ).string();
    mesh->saveHDF5( filename );
    for ( bool repartition : { false, true } )
    {
        auto loaded = loadMesh( _mesh=new mesh_type( RuntimeOrder{1} ),
                               _filename=filename, _straighten=false,
                               _rebuild_partitions=repartition,
                               _rebuild_partitions_filename="runtime_circle_repartitioned.json" );
        BOOST_CHECK_EQUAL( loaded->order(), 2 );
        BOOST_CHECK_EQUAL( loaded->numGlobalElements(), mesh->numGlobalElements() );
        BOOST_CHECK_CLOSE( loaded->measure(), mesh->measure(), 1e-9 );
        for ( auto const& [id, element] : loaded->elements() )
            BOOST_CHECK_EQUAL( element.nPoints(), 6 );
    }
}
#endif

BOOST_AUTO_TEST_SUITE_END()
