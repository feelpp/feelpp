#define BOOST_TEST_MODULE test_geotool
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelfilters/straightenmesh_impl.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>

using namespace Feel;
using namespace Feel::vf;



namespace test_geotool
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
    po::options_description desc_options( "test_geotool options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.03 ), "mesh size" )
    ( "hsize3d", po::value<double>()->default_value( 0.1 ), "mesh size" )
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
    AboutData about( "test_geotool" ,
                     "test_geotool" ,
                     "0.1",
                     "test_geotool",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@feelpp.org", "" );
    return about;
}

//-----------------------------------------------------//

void runRectangle()
{
    typedef Mesh<Simplex<2,1,2> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runRectangle" );

    double meshSize = option(_name="hsize").as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary1",_marker3=true);
    R.setMarker(_type="line",_name="Boundary2",_marker1=true,_marker2=true,_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);

    auto mesh = R.createMesh(_mesh=new mesh_type,
                             _name="domainRectangle" );

    auto area = integrate( _range=elements( mesh ),
                           _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( area, 4., 1e-9 );

    auto lg1 = integrate( _range=markedfaces( mesh,"Boundary1" ),
                          _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( lg1, 4., 1e-9 );

    auto lg2 = integrate( _range=markedfaces( mesh,"Boundary2" ),
                          _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( lg2, 6., 1e-9 );

    auto Xh = Pch<3>(mesh);
}

//-----------------------------------------------------//

void runRectangleWithPerfo()
{
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
    typedef Mesh<Simplex<2,2,2> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runRectangleWithPerfo" );

    double meshSize = option(_name="hsize").as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary1",_markerAll=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);

    GeoTool::Node xc(0.2,0.2);
    GeoTool::Node xr(0.3,0.2);
    GeoTool::Circle Circ(meshSize/3.,"OmegaCircle",xc,xr);
    Circ.setMarker(_type="line",_name="Boundary2",_markerAll=true);
    Circ.setMarker(_type="surface",_name="Omega",_markerAll=true);

    GeoTool::Node xc2(0.7,0.7);
    GeoTool::Node xr2(0.9,0.7);
    GeoTool::Circle Circ2(meshSize/3.,"OmegaCircle2",xc2,xr2);
    Circ2.setMarker(_type="line",_name="Boundary2",_markerAll=true);
    Circ2.setMarker(_type="surface",_name="Omega",_markerAll=true);

    auto mesh = (R-(Circ+Circ2)).createMesh(_mesh=new mesh_type,
                                            _name="domainRectangleWithPerfo" );


    auto area = integrate( _range=elements( mesh ),
                           _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( area, 1.-M_PI*0.1*0.1-M_PI*0.2*0.2, 1e-5 );

    auto lg1 = integrate( _range=markedfaces( mesh,"Boundary1" ),
                          _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( lg1, 4., 1e-9 );

    auto lg2 = integrate( _range=markedfaces( mesh,"Boundary2" ),
                          _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( lg2, 2*M_PI*(0.1+0.2), 1e-5 );

    auto Xh = Pch<3>(mesh);
#endif
}

//-----------------------------------------------------//

void runRectangleWithConformalInternalShape()
{
    typedef Mesh<Simplex<2,1,2> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runRectangleWithConformalInternalShape" );

    double meshSize = option(_name="hsize").as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary1",_markerAll=true);
    R.setMarker(_type="surface",_name="Omega1",_markerAll=true);

    GeoTool::Node xc(0.5,0.5);
    GeoTool::Node xr(0.4,0.5);
    GeoTool::Circle Circ(meshSize/3.,"OmegaCircle",xc,xr);
    Circ.setMarker(_type="line",_name="Boundary2",_markerAll=true);
    Circ.setMarker(_type="surface",_name="Omega2",_markerAll=true);

    GeoTool::Node xr2(0.3,0.5);
    GeoTool::Circle Circ2(meshSize/3.,"OmegaCircle2",xc,xr2);
    Circ2.setMarker(_type="line",_name="Boundary3",_markerAll=true);
    Circ2.setMarker(_type="surface",_name="Omega3",_markerAll=true);

    GeoTool::Node xc3(0.2,0.2);
    GeoTool::Node xr3(0.3,0.2);
    GeoTool::Circle Circ3(meshSize/3.,"OmegaCircle3",xc3,xr3);
    Circ3.setMarker(_type="line",_name="Boundary4",_markerAll=true);
    Circ3.setMarker(_type="surface",_name="Omega4",_markerAll=true);


    auto mesh = ((R-(Circ2+Circ3))+(Circ2-Circ)+Circ+Circ3).createMesh(_mesh=new mesh_type,
                                                                       _name="domainRectangleWithConformalInternalShape" );

    auto area = integrate( _range=elements( mesh ),
                           _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( area, 1, 1e-9 );

    auto lg1 = integrate( _range=markedfaces( mesh,"Boundary1" ),
                          _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( lg1, 4., 1e-9 );

    auto lg2 = integrate( _range=markedfaces( mesh,"Boundary2" ),
                          _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( lg2, 2*M_PI*0.1, 1e-1 );

    auto lg3 = integrate( _range=markedfaces( mesh,"Boundary3" ),
                          _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( lg3, 2*M_PI*0.2, 1e-1 );

    auto Xh = Pch<3>(mesh);
}

//-----------------------------------------------------//

void runSpecial1a()
{
    typedef Mesh<Simplex<2,1,2> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runSpecial1a" );

    double meshSize = option(_name="hsize").as<double>()*4;

    GeoTool::Node x1( -0.3,-2.3 );
    GeoTool::Node x2( 11.3,2.3 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary1",_marker3=true);
    R.setMarker(_type="line",_name="Boundary2",_marker1=true,_marker2=true,_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);

    GeoTool::Special_1a S( meshSize,"OMEGA2",x1 );
    S.setMarker(_type="line",_name="Boundary1",_markerAll=true);
    S.setMarker(_type="surface",_name="Omega",_markerAll=true);

    auto mesh = (R-S).createMesh(_mesh=new mesh_type,
                                 _name="domainSpecial1a" );
    auto Xh = Pch<3>(mesh);
}

//-----------------------------------------------------//

void runHexahedron()
{

    typedef Mesh<Simplex<3,1,3> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runHexahedron" );

    double meshSize = option(_name="hsize3d").as<double>();

    GeoTool::Node x1(-1,-1,-1);
    GeoTool::Node x2( 1,-1,-1);
    GeoTool::Node x3( 1, 1,-1);
    GeoTool::Node x4(-1, 1,-1);
    GeoTool::Node x5(-1,-1, 1);
    GeoTool::Node x6( 1,-1, 1);
    GeoTool::Node x7( 1, 1, 1);
    GeoTool::Node x8(-1, 1, 1);
    GeoTool::Hexahedron H( meshSize,"Hex",x1,x2,x3,x4,x5,x6,x7,x8);
    H.setMarker(_type="surface",_name="Crush",_marker5=true);
    H.setMarker(_type="surface",_name="Free",_marker1=true);
    H.setMarker(_type="surface",_name="Fixed",
                _marker2=true,
                _marker3=true,
                _marker4=true,
                _marker6=true);
    H.setMarker(_type="volume",_name="OmegaSolid",_markerAll=true);

    auto mesh = H.createMesh(_mesh=new mesh_type,
                             _name="domainHexahedron" );

    auto area = integrate( _range=elements( mesh ),
                           _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( area, 8., 1e-9 );

    typedef Mesh<Simplex<2,1,3> > mesh_surf_type;
    auto meshSurf = H.createMesh(_mesh=new mesh_surf_type,
                                 _name="domainHexahedronSurf" );
    auto Xh = Pch<3>(meshSurf);
}

//-----------------------------------------------------//

void runCylinder()
{
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
    typedef Mesh<Simplex<3,2,3> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runCylinder" );

    double meshSize = option(_name="hsize3d").as<double>();

    GeoTool::Node Centre(0,0,0);
    GeoTool::Node Rayon( 0.5);
    GeoTool::Node Dir(1,0,0);
    GeoTool::Node Lg(5,0,0);
    GeoTool::Cylindre C( meshSize,"Cyl",Centre,Dir,Rayon,Lg);
    C.setMarker(_type="surface",_name="Inlet",_marker1=true);
    C.setMarker(_type="surface",_name="Outlet",_marker2=true);
    C.setMarker(_type="surface",_name="Wall",_marker3=true);
    C.setMarker(_type="volume",_name="Omega",_markerAll=true);

    auto mesh = C.createMesh(_mesh=new mesh_type,
                             _name="domainCylinder" );

    auto area = integrate( _range=elements( mesh ),
                           _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( area, M_PI*0.5*0.5*5, 1e-3 );

    typedef Mesh<Simplex<2,2,3> > mesh_surf_type;
    auto meshSurf = C.createMesh(_mesh=new mesh_surf_type,
                             _name="domainCylinderSurf" );
    auto Xh = Pch<3>(meshSurf);
#endif
}

//-----------------------------------------------------//

void runCubeWithPerfo()
{
    typedef Mesh<Simplex<3,1,3> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runCubeWithPerfo" );

    double meshSize = option(_name="hsize3d").as<double>();

    GeoTool::Node x1(0,0,0);
    GeoTool::Node x2(1,1,1);
    GeoTool::Cube C( meshSize,"Cube1",x1,x2);
    C.setMarker(_type="surface",_name="Boundary1",_marker1=true);
    C.setMarker(_type="surface",_name="Boundary2",_marker2=true);
    C.setMarker(_type="surface",_name="Boundary3",_marker3=true,_marker4=true,_marker5=true,_marker6=true);
    C.setMarker(_type="volume",_name="Omega",_markerAll=true);

    GeoTool::Node x3(0.3,0.3,0.3);
    GeoTool::Node x4(0.7,0.7,0.7);
    GeoTool::Cube C2( meshSize,"Cube2",x3,x4);
    C2.setMarker(_type="surface",_name="Boundary4",_markerAll=true);
    C2.setMarker(_type="volume",_name="Omega",_markerAll=true);

    auto mesh = (C-C2).createMesh(_mesh=new mesh_type,
                                  _name="domainCubeWithPerfo" );

    auto area = integrate( _range=elements( mesh ),
                           _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( area, 1-0.4*0.4*0.4, 1e-8 );

    typedef Mesh<Simplex<2,1,3> > mesh_surf_type;
    auto meshSurf = (C/*-C2*/).createMesh(_mesh=new mesh_surf_type,
                                      _name="domainCubeWithPerfoSurf" );
    auto Xh = Pch<3>(meshSurf);
}

//-----------------------------------------------------//

void runCubeWithConformalInternalShape()
{
    typedef Mesh<Simplex<3,1,3> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runCubeWithConformalInternalShape" );

    double meshSize = option(_name="hsize3d").as<double>();

    GeoTool::Node x1(0,0,0);
    GeoTool::Node x2(1,1,1);
    GeoTool::Cube C( meshSize,"Cube1",x1,x2);
    C.setMarker(_type="surface",_name="Boundary1",_marker1=true);
    C.setMarker(_type="surface",_name="Boundary2",_marker2=true);
    C.setMarker(_type="surface",_name="Boundary3",_marker3=true,_marker4=true,_marker5=true,_marker6=true);
    C.setMarker(_type="volume",_name="Omega",_markerAll=true);

    GeoTool::Node x3(0.3,0.3,0.3);
    GeoTool::Node x4(0.7,0.7,0.7);
    GeoTool::Cube C2( meshSize,"Cube2",x3,x4);
    C2.setMarker(_type="surface",_name="Boundary4",_markerAll=true);
    C2.setMarker(_type="volume",_name="Omega",_markerAll=true);

    auto mesh = ((C-C2)+C2).createMesh(_mesh=new mesh_type,
                                       _name="domainCubeWithConformalInternalShape",
                                       _optimize3d_netgen=false // else bug with last gmsh-devel(svn.rev:13736)
                                       );

    auto area = integrate( _range=elements( mesh ),
                           _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( area, 1., 1e-8 );

    typedef Mesh<Simplex<2,1,3> > mesh_surf_type;
    auto meshSurf = ((C-C2)+C2).createMesh(_mesh=new mesh_surf_type,
                                       _name="domainCubeWithConformalInternalShapeSurf",
                                       _optimize3d_netgen=false // else bug with last gmsh-devel(svn.rev:13736)
                                       );
    auto Xh = Pch<3>(meshSurf);
}

//-----------------------------------------------------//
//-----------------------------------------------------//
//-----------------------------------------------------//

void runGeneric3d(typename GeoTool::GeoGMSHTool geo, std::string name)
{
    typedef Mesh<Simplex<3,1,3> > mesh_type;
    auto mesh = geo.createMesh(_mesh=new mesh_type,
                               _name=name );

    typedef Mesh<Simplex<2,1,3> > mesh_surf_type;
    auto meshSurf = geo.createMesh(_mesh=new mesh_surf_type,
                              _name=name+"Surf" );
    auto Xh = Pch<3>(meshSurf);
}


#if 0
void runSphereHollow()
{
    typedef Mesh<Simplex<3,1,3> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runSphereHollow" );

    double meshSize = option(_name="hsize3d").as<double>()/3.;

   GeoTool::Sphere S1(meshSize,"S1",
                       GeoTool::Node(0.5,0,0),
                       GeoTool::Node(0.1,0,0) );
    S1.setMarker(_type="surface",_name="WallFBM",_markerAll=true);
    S1.setMarker(_type="volume",_name="Sphere",_markerAll=true);
    GeoTool::Sphere S2(meshSize,"S2",
                       GeoTool::Node(0.5,0,0),
                       GeoTool::Node(0.2,0,0) );
    S2.setMarker(_type="surface",_name="FictitiousBoundary",_markerAll=true);
    S2.setMarker(_type="volume",_name="Sphere",_markerAll=true);

    runGeneric3d(S2-S1,"domainSphere");
}
#endif

#if 0
void runTetrahedron()
{
    typedef Mesh<Simplex<3,1,3> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runTetrahedron" );

    double meshSize = option(_name="hsize3d").as<double>()/3.;

    GeoTool::Node x1( -1,-1,-1 );
    GeoTool::Node x2(  1,-1,-1 );
    GeoTool::Node x3(  0, 1,-1 );
    GeoTool::Node x4(  0, 0, 1 );
    GeoTool::Tetrahedron R( meshSize,"OMEGA",x1,x2,x3,x4 );
    R.setMarker(_type="surface",_name="Boundary1",_marker1=true);
    R.setMarker(_type="surface",_name="Boundary2",_marker2=true);
    R.setMarker(_type="surface",_name="Boundary3",_marker3=true);
    R.setMarker(_type="surface",_name="Boundary4",_marker4=true);
    R.setMarker(_type="volume",_name="Omega",_markerAll=true);

    runGeneric3d(R,"domainTetra");
}
#endif

#if 0
void runSpecial3D_1()
{
    typedef Mesh<Simplex<3,1,3> > mesh_type;

    if ( Environment::worldComm().isMasterRank() ) BOOST_MESSAGE( "runSpecial3D_1" );

    double meshSize = option(_name="hsize3d").as<double>()/3.;

    GeoTool::Special3D_1 S( meshSize/3.,"shapeStruct");
    S.setMarker(_type="surface",_name="fsiwall",_marker1=true);
    S.setMarker(_type="surface",_name="wallcylinder",_marker2=true);
    S.setMarker(_type="volume",_name="OmegaStruct",_markerAll=true);

    runGeneric3d(S,"domainSpecial3D_1");
}
#endif


} // namespace test_geotool

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_geotool::makeAbout(), test_geotool::makeOptions() )

BOOST_AUTO_TEST_SUITE( interp_geotool )

BOOST_AUTO_TEST_CASE( interp_geotool )
{
    test_geotool::runRectangle();
    test_geotool::runRectangleWithPerfo();
    test_geotool::runRectangleWithConformalInternalShape();
    test_geotool::runSpecial1a();
    test_geotool::runHexahedron();
    test_geotool::runCylinder();
    test_geotool::runCubeWithPerfo();
    test_geotool::runCubeWithConformalInternalShape();
#if 0
    test_geotool::runSphereHollow();
    test_geotool::runTetrahedron();
    test_geotool::runSpecial3D_1();
#endif

}

BOOST_AUTO_TEST_SUITE_END()



