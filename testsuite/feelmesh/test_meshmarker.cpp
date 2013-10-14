
#define BOOST_TEST_MODULE test_meshmarker
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>



using namespace Feel;
using namespace Feel::vf;

namespace test_meshmarker
{

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_meshmarker options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ;
    return desc_options.add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "Test_meshmarker" ,
                     "Test_meshmarker" ,
                     "0.1",
                     "test meshmarker",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2012 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}

void
test_meshmarker( Application_ptrtype test_app )
{
    typedef Mesh<Simplex<3,1,3> > mesh_type;

    auto meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1(-1,-1,-1);
    GeoTool::Node x2( 1,-1,-1);
    GeoTool::Node x3( 1, 1,-1);
    GeoTool::Node x4(-1, 1,-1);
    GeoTool::Node x5(-1,-1, 1);
    GeoTool::Node x6( 1,-1, 1);
    GeoTool::Node x7( 1, 1, 1);
    GeoTool::Node x8(-1, 1, 1);
    GeoTool::Hexahedron H(meshSize,"UnHexa",x1,x2,x3,x4,x5,x6,x7,x8);
    H.setMarker(_type="surface",_name="GammaDirichlet",_marker2=true,_marker3=true,_marker4=true,_marker5=true,_marker6=true);
    H.setMarker(_type="surface",_name="GammaNeumann",_marker1=true);
    H.setMarker(_type="volume",_name="OmegaFluid",_markerAll=true);

    auto mesh = H.createMesh(_mesh= new mesh_type,
                             _name="un_cube" );


    mesh->updateMarker3WithRange(markedfaces(mesh,"GammaNeumann"),1);
    mesh->updateMarkersFromFaces();

    auto submesh = createSubmesh(mesh,marked3elements(mesh,1));

    auto myexporter = Exporter<mesh_type>::New( test_app->vm(), "test_meshmarker_export" );
    //myexporter->step(0)->setMesh( mesh );
    myexporter->step(0)->setMesh( submesh );
    //myexporter->step(0)->add( "test2dP1mesh_uP1", u );
    myexporter->save();
}


} // namespace test_meshmarker

/**
 * main code
 */
FEELPP_ENVIRONMENT_WITH_OPTIONS( test_meshmarker::makeAbout(), test_meshmarker::makeOptions() )

BOOST_AUTO_TEST_SUITE( test_meshmarker )

BOOST_AUTO_TEST_CASE( test_meshmarker1 )
{
    using namespace Feel::vf;
    using namespace test_meshmarker;

    Application_ptrtype test_app( new Application_type( boost::unit_test::framework::master_test_suite().argc,
                                  boost::unit_test::framework::master_test_suite().argv,
                                  test_meshmarker::makeAbout(),
                                  test_meshmarker::makeOptions()
                                                      ) );

    test_app->changeRepository( boost::format( "/testsuite/feelmesh/%1%/" )
                                % test_app->about().appName() );

    test_meshmarker::test_meshmarker( test_app);
}

BOOST_AUTO_TEST_SUITE_END()

