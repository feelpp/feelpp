/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/
#define BOOST_TEST_MODULE test_masterworker_laplacian
#include <feel/feelalg/backend.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

namespace test_masterworker_laplacian
{
using namespace Feel;
using namespace Feel::vf;

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline po::options_description
makeOptions()
{
    po::options_description desc_options( "test_twolaplaciansdistributed options" );
    desc_options.add_options()( "hsize", po::value<double>()->default_value( 0.02 ), "mesh size" )( "1d-hsize", po::value<double>()->default_value( 0.02 ), "mesh size 1d" );
    return desc_options.add( Feel::feel_options() ).add( backend_options( "backend1" ) );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline AboutData
makeAbout()
{
    AboutData about( "Test_masterworker_laplacian",
                     "Test_masterworker_laplacian",
                     "0.1",
                     "test master worker laplacian",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2020 Feel++ Consortium" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@feelpp.org", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

void run()
{
    typedef Mesh<Simplex<2, 1, 2>> mesh_type;
    double meshSize = doption( _name = "hsize" );

    const int nTotalProc = Environment::worldComm().size();

    auto worldCommDomain1 = Environment::worldCommPtr();
    auto worldCommDomain2 = Environment::worldCommPtr();

    if ( nTotalProc > 1 )
    {
        std::vector<int> MapWorld( nTotalProc, 1 );
        MapWorld[Environment::masterRank()] = 0;
        auto worldCommColored = WorldComm( MapWorld );
        //auto worldCommColored = WorldComm( Environment::isMasterRank()?0:1 );
        //worldCommColored.showMe();
        worldCommDomain1 = worldCommColored.subWorldComm( 0 );
        worldCommDomain2 = worldCommColored.subWorldComm( 1 );
    }
#if 0    
    worldCommDomain1->showMe();
    worldCommDomain2->showMe();
    Environment::worldCommPtr()->showMe();
#endif    
    //---------------------------------------------------------------------------------------//

    if ( worldCommDomain2->isActive() )
    {
        std::cerr << "rank: { local: " << worldCommDomain2->localRank() << "},{ global: " << worldCommDomain2->globalRank() << "} } " << std::endl;
        GeoTool::Node x21( -1, 0.3 );
        GeoTool::Node x22( -1 + 0.5, 0.3 );
        GeoTool::Circle Omega2( meshSize, "Omega", x21, x22 );
        Omega2.setMarker( _type = "line", _name = "Paroi", _markerAll = true );
        Omega2.setMarker( _type = "surface", _name = "Omega", _markerAll = true );
        auto mesh2 = Omega2.createMesh( _mesh = new mesh_type,
                                        _name = "omega2_" + mesh_type::shape_type::name(),
                                        _partitions = worldCommDomain2->localSize(),
                                        _worldcomm = worldCommDomain2 );
        //---------------------------------------------------------------------------------------//
        // functionspaces
        auto Xh2 = Pch<1>( mesh2 );
        auto u2 = Xh2->element();
        //---------------------------------------------------------------------------------------//
        // backends
        auto backend2 = backend( _rebuild = true, _worldcomm = Xh2->worldCommPtr() );
        //---------------------------------------------------------------------------------------//
        // init matrix and vectors
        auto A2 = backend2->newMatrix( _test = Xh2, _trial = Xh2 );
        auto F2 = backend2->newVector( Xh2 );

        //---------------------------------------------------------------------------------------//
        // assembly
        form2( _test = Xh2, _trial = Xh2, _matrix = A2 ) +=
            integrate( _range = elements( mesh2 ), _expr = gradt( u2 ) * trans( grad( u2 ) ) );
        form1( _test = Xh2, _vector = F2 ) +=
            integrate( _range = elements( mesh2 ), _expr = id( u2 ) );
        form2( _test = Xh2, _trial = Xh2, _matrix = A2 ) +=
            on( _range = boundaryfaces( mesh2 ),
                _element = u2, _rhs = F2,
                _expr = cst( 0. ) );

        //---------------------------------------------------------------------------------------//
        // solve
        backend2->solve( _matrix = A2, _solution = u2, _rhs = F2 /*,_prec=prec2*/ );
        double l2NormU2 = u2.l2Norm();
        BOOST_CHECK( l2NormU2 > 1e-5 );

        //---------------------------------------------------------------------------------------//
        // exports
        auto myexporter2 = exporter( _mesh = mesh2, _name = "MyExportDomain2" );
        myexporter2->step( 0 )->add( "u2", u2 );
        myexporter2->save();
    }
    Environment::worldCommPtr()->barrier();
}

} // namespace test_masterworker
 
FEELPP_ENVIRONMENT_WITH_OPTIONS( test_masterworker_laplacian::makeAbout(),
                                 test_masterworker_laplacian::makeOptions() )

BOOST_AUTO_TEST_SUITE( test_masterworker_laplacian )

BOOST_AUTO_TEST_CASE( test_masterworker_laplacian_1 )
{
    test_masterworker_laplacian::run();
}
BOOST_AUTO_TEST_SUITE_END()
