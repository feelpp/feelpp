#define BOOST_TEST_MODULE test_normal3d
#include <testsuite/testsuite.hpp>


#include <feel/options.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <iostream>

using namespace Feel;

namespace test_normal3d
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
    po::options_description desc_options( "test_normal3d options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ( "geomap", po::value<int>()->default_value( ( int )GeomapStrategyType::GEOMAP_OPT ), "geomap (0=opt, 1=p1, 2=ho)" )
    ( "straighten", po::value<int>()->default_value( 1 ), "straighten mesh" )
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
    AboutData about( "Test_normal3d" ,
                     "Test_normal3d" ,
                     "0.1",
                     "test normal 3d",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}

template <uint32_type orderGeo = 1>
void
runtest( Application_ptrtype test_app )
{
    BOOST_MESSAGE( "================================================================================\n"
                   << "Order: " << orderGeo << "\n" );

    typedef Simplex<3,orderGeo,3> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    typedef bases<Lagrange<2+orderGeo,Vectorial,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //typedef typename space_type::element_type element_type;

    //-----------------------------------------------------------//

    double meshSize = test_app->vm()["hsize"].as<double>();
    bool exportResults = test_app->vm()["exporter.export"].as<bool>();
    int straighten = test_app->vm()["straighten"].as<int>();
    GeomapStrategyType geomap = ( GeomapStrategyType )test_app->vm()["geomap"].as<int>();
    GeoTool::Node Centre( 0,0,0 );
    GeoTool::Node Rayon( 1 );
    GeoTool::Node Dir( 1,0,0 );
    GeoTool::Node Lg( 3,0,0 );
    GeoTool::Cylindre C( meshSize,"UnCylindre",Centre,Dir,Rayon,Lg );
    C.setMarker( _type="surface",_name="Inlet",_marker1=true );
    C.setMarker( _type="surface",_name="Outlet",_marker2=true );
    C.setMarker( _type="surface",_name="Wall",_marker3=true );
    C.setMarker( _type="volume",_name="OmegaFluid",_markerAll=true );

    std::ostringstream __ostrData;
    __ostrData<<convex_type::name()<<"h"<<meshSize;

    auto mesh_ = C.createMesh(_mesh=new mesh_type,_name= "domain"+__ostrData.str(), _straighten=false );
    auto mesh=mesh_;

    if ( straighten )
        mesh = straightenMesh( mesh_, Environment::worldComm(), false, true );

    //-----------------------------------------------------------//

    auto Xh = space_type::New( mesh );
    auto u = Xh->element();

    BOOST_MESSAGE( "testing Gauss formula on ( 1, 1, 1 )\n" );
    u = project( Xh,elements( mesh ),vec( cst( 1.0 ),cst( 1.0 ),cst( 1.0 ) ) );
    auto value1 = integrate( _range=elements( mesh ),_expr=divv( u )/*, _quad=_Q<8>(),_quad1=_Q<8>(),*/,_geomap=geomap ).evaluate()( 0,0 );
    auto value2 = integrate( _range=boundaryfaces( mesh ),_expr=trans( idv( u ) )*N()/*,_quad=_Q<8>(),_quad1=_Q<8>()*/,_geomap=geomap ).evaluate()( 0,0 );
    BOOST_MESSAGE( "\n value (div) =" << value1 << "\n value (n) =" << value2 <<"\n" );
    BOOST_CHECK_SMALL( value1, 1e-12 );
    BOOST_CHECK_SMALL( value2, 1e-12 );

    BOOST_MESSAGE( "testing Gauss formula on ( cos(M_PI*Px()/5.),cos(M_PI*Py()/5.),cos(M_PI*Py()/5.))\n" );
    u = project( Xh,elements( mesh ),vec( cos( M_PI*Px()/5. ),cos( M_PI*Py()/5. ),cos( M_PI*Py()/5. ) ) );
    value1 = integrate( _range=elements( mesh ),_expr=divv( u ), _quad=_Q<8>(),_quad1=_Q<8>(),_geomap=geomap ).evaluate()( 0,0 );
    value2 = integrate( _range=boundaryfaces( mesh ),_expr=trans( idv( u ) )*N(),_quad=_Q<8>(),_quad1=_Q<8>(),_geomap=geomap ).evaluate()( 0,0 );
    BOOST_MESSAGE( "\n value (div) =" << value1 << "\n value (n) =" << value2 <<"\n" );
    BOOST_CHECK_CLOSE( value1, value2, 1e-8 );
    //BOOST_CHECK_SMALL( value1-value2,1e-8);


    if ( exportResults )
    {
        auto nnn = Xh->element();
        nnn = project( Xh,markedfaces( mesh,"Wall" ),N() );
        auto UNexporter = Exporter<mesh_type>::New( "gmsh"/*test_app->vm()*/, "ExportOOOO"+__ostrData.str() );
        UNexporter->step( 0 )->setMesh( mesh );
        UNexporter->step( 0 )->add( "u", u );
        UNexporter->step( 0 )->add( "n", nnn );
        UNexporter->save();
    }
}

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_normal3d::makeAbout(), test_normal3d::makeOptions() )

BOOST_AUTO_TEST_SUITE( normal3d )

BOOST_AUTO_TEST_CASE( normal3d )
{

    using namespace test_normal3d;

    auto test_app = Application_ptrtype( new Application_type );

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/" )
                                % test_app->about().appName()
                              );


    runtest<1>( test_app );
    runtest<2>( test_app );
    //runtest<3>( test_app );
    //runtest<4>( test_app );


}

BOOST_AUTO_TEST_SUITE_END()



