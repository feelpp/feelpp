
#if 1
#define BOOST_TEST_MODULE test_geotool2
#include <testsuite/testsuite.hpp>
#endif

#if defined(USE_BOOST_TEST)
#define TEST_GEOTOOL2_CHECK_CLOSE(val1,val2,tol) \
    BOOST_CHECK_CLOSE( val1, val2, tol )         \
    /**/
#else
#define TEST_GEOTOOL2_CHECK_CLOSE(val1,val2,tol) \
  CHECK( math::abs( (val1)-(val2) ) < tol ) << "error " << " val1 " << val1 << " val2 " << val2 << " tol " << tol \
  /**/
#endif


#if defined(USE_BOOST_TEST)
#define TEST_GEOTOOL2_MESSAGE(msg)                  \
  if ( Environment::worldComm().isMasterRank() )    \
      BOOST_MESSAGE( msg );                         \
  /**/
#else
#define TEST_GEOTOOL2_MESSAGE(msg)                  \
  if ( Environment::worldComm().isMasterRank() )    \
      std::cout << msg << "\n";                     \
  /**/
#endif

#include <feel/options.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>

using namespace Feel;
using namespace Feel::vf;



namespace test_geotool2
{

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_geotool options" );
    desc_options.add_options()
    ( "hsize2d", po::value<double>()->default_value( 0.07 ), "mesh size" )
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
    AboutData about( "test_geotool2" ,
                     "test_geotool2" ,
                     "0.1",
                     "test_geotool",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@feelpp.org", "" );
    return about;
}

//-----------------------------------------------------//


void runFusion2d( bool keepInterface )
{
    std::string keepInterfaceStr( (keepInterface)? "with-interface" : "without-interface" );

    TEST_GEOTOOL2_MESSAGE("runFusion2d-"+keepInterfaceStr)

    double meshSize = doption(_name="hsize2d");
    typedef Mesh<Simplex<2,1,2> > mesh_type;

    GeoTool::Node x1a( 0,0 );
    GeoTool::Node x2a( 1,1 );
    GeoTool::Rectangle Ra( meshSize,"OMEGA",x1a,x2a );
    Ra.setMarker(_type="line",_name="Boundary1",_marker1=true,_marker3=true,_marker4=true);
    Ra.setMarker(_type="line",_name="InternalInterface",_marker2=true);
    Ra.setMarker(_type="surface",_name="Omega1",_markerAll=true);

    GeoTool::Node x1b( 1,0 );
    GeoTool::Node x2b( 2,1 );
    GeoTool::Rectangle Rb( meshSize/2.,"OMEGA2",x1b,x2b );
    Rb.setMarker(_type="line",_name="Boundary1",_marker1=true,_marker3=true);
    Rb.setMarker(_type="line",_name="InternalInterface",_marker2=true,_marker4=true);
    Rb.setMarker(_type="surface",_name="Omega1",_markerAll=true);

    GeoTool::Node x1e( 2,0 );
    GeoTool::Node x2e( 3,1 );
    GeoTool::Rectangle Rc( meshSize/3.,"OMEGA11",x1e,x2e );
    Rc.setMarker(_type="line",_name="Boundary2",_marker1=true,_marker3=true);
    Rc.setMarker(_type="line",_name="InternalInterface",_marker2=true,_marker4=true);
    Rc.setMarker(_type="surface",_name="Omega2",_markerAll=true);

    GeoTool::Node x1f( 3,0 );
    GeoTool::Node x2f( 4,1 );
    GeoTool::Rectangle Rd( meshSize/2.,"OMEG",x1f,x2f );
    Rd.setMarker(_type="line",_name="Boundary2",_marker1=true,_marker2=true,_marker3=true);
    Rd.setMarker(_type="line",_name="InternalInterface",_marker4=true);
    Rd.setMarker(_type="surface",_name="Omega2",_markerAll=true);

    auto mesh1 = (Ra+Rb+Rc+Rd).
        fusion(Ra,2,Rb,4,keepInterface).
        fusion(Rb,2,Rc,4,true).
        fusion(Rc,2,Rd,4,keepInterface).
        createMesh(_mesh=new mesh_type,
                   _name="domain1"+keepInterfaceStr );

    auto area1elt1 = integrate( _range=markedelements( mesh1, "Omega1" ),
                                _expr=cst(1.) ).evaluate()( 0,0 );
    TEST_GEOTOOL2_CHECK_CLOSE( area1elt1, 2., 1e-9 );
    auto area1elt2 = integrate( _range=markedelements( mesh1, "Omega2" ),
                                _expr=cst(1.) ).evaluate()( 0,0 );
    TEST_GEOTOOL2_CHECK_CLOSE( area1elt2, 2., 1e-9 );
    auto area1face1 = integrate( _range=markedfaces( mesh1, "Boundary1" ),
                                 _expr=cst(1.) ).evaluate()( 0,0 );
    TEST_GEOTOOL2_CHECK_CLOSE( area1face1, 5., 1e-9 );
    auto area1face2 = integrate( _range=markedfaces( mesh1, "Boundary2" ),
                                 _expr=cst(1.) ).evaluate()( 0,0 );
    TEST_GEOTOOL2_CHECK_CLOSE( area1face2, 5., 1e-9 );

    if ( keepInterface )
    {
        auto area1face3 = integrate( _range=markedfaces( mesh1, "InternalInterface" ),
                                     _expr=cst(1.) ).evaluate()( 0,0 );
        TEST_GEOTOOL2_CHECK_CLOSE( area1face3, 3., 1e-9 );
    }

    GeoTool::Node xca(0.5,0.5);
    GeoTool::Node xra(0.7,0.5);
    GeoTool::Circle Ca(meshSize/6.,"OmegaCircle1",xca,xra);
    Ca.setMarker(_type="line",_name="BoundaryCirc1",_markerAll=true);
    Ca.setMarker(_type="surface",_name="Omega",_markerAll=true);

    GeoTool::Node xcb(1+0.5,0.5);
    GeoTool::Node xrb(1+0.7,0.5);
    GeoTool::Circle Cb(meshSize/6.,"OmegaCircle2",xcb,xrb);
    Cb.setMarker(_type="line",_name="BoundaryCirc2",_markerAll=true);
    Cb.setMarker(_type="surface",_name="Omega",_markerAll=true);

    auto mesh2 = ((Ra-Ca)+(Rb-Cb)).
        fusion(Ra,2,Rb,4,keepInterface).
        createMesh(_mesh=new mesh_type,
                   _name="domain2"+keepInterfaceStr );

    auto area2elt = integrate( _range=markedelements( mesh2, "Omega1" ),
                                _expr=cst(1.) ).evaluate()( 0,0 );
    TEST_GEOTOOL2_CHECK_CLOSE( area2elt, 2.0-2*M_PI*0.04, 1e-2 );
    auto area2face1 = integrate( _range=boundaryfaces( mesh2 ),
                                _expr=cst(1.) ).evaluate()( 0,0 );
    TEST_GEOTOOL2_CHECK_CLOSE( area2face1, 6+4*M_PI*0.2, 1e-2 );


    GeoTool::Node x1g( 0.2,0.3 );
    GeoTool::Node x2g( 0.5,0.7 );
    GeoTool::Rectangle Rg( meshSize/3.,"OMEGA11g",x1g,x2g );
    Rg.setMarker(_type="line",_name="BoundaryRect1",_marker1=true,_marker3=true,_marker4=true);
    Rg.setMarker(_type="line",_name="InternalInterfaceRect1",_marker2=true);
    Rg.setMarker(_type="surface",_name="Omega2",_markerAll=true);

    GeoTool::Node x1h( 0.5,0.3 );
    GeoTool::Node x2h( 0.8,0.7 );
    GeoTool::Rectangle Rh( meshSize/2.,"OMEGh",x1h,x2h );
    Rh.setMarker(_type="line",_name="BoundaryRect1",_marker1=true,_marker2=true,_marker3=true);
    Rh.setMarker(_type="line",_name="InternalInterfaceRect1",_marker4=true);
    Rh.setMarker(_type="surface",_name="Omega2",_markerAll=true);

    GeoTool::Node x1i( 1+0.2,0.3 );
    GeoTool::Node x2i( 1+0.5,0.7 );
    GeoTool::Rectangle Ri( meshSize/3.,"OMEGA11i",x1i,x2i );
    Ri.setMarker(_type="line",_name="BoundaryRect2",_marker1=true,_marker3=true,_marker4=true);
    Ri.setMarker(_type="line",_name="InternalInterfaceRect2",_marker2=true);
    Ri.setMarker(_type="surface",_name="Omega2",_markerAll=true);

    GeoTool::Node x1j( 1+0.5,0.3 );
    GeoTool::Node x2j( 1+0.8,0.7 );
    GeoTool::Rectangle Rj( meshSize/2.,"OMEGj",x1j,x2j );
    Rj.setMarker(_type="line",_name="BoundaryRect2",_marker1=true,_marker2=true,_marker3=true);
    Rj.setMarker(_type="line",_name="InternalInterfaceRect2",_marker4=true);
    Rj.setMarker(_type="surface",_name="Omega2",_markerAll=true);

    auto mesh3 =
        //((Ra-(Rg+Rh))+Rb).
        //((Ra-(Rg+Rh))+(Rb-(Ri+Rj))).
        ((Ra-(Rg+Rh))+(Rb-(Ri+Rj))+Rc+Rd).
        fusion(Ra,2,Rb,4,keepInterface).
        fusion(Rg,2,Rh,4,false/*keepInterface*/).
        fusion(Ri,2,Rj,4,false/*keepInterface*/).
        fusion(Rb,2,Rc,4,keepInterface).
        fusion(Rc,2,Rd,4,keepInterface).
        createMesh(_mesh=new mesh_type,
                   _name="domain3"+keepInterfaceStr );

    if ( keepInterface )
    {
        auto area3elt1 = integrate( _range=markedelements( mesh3, "Omega1" ),
                                    _expr=cst(1.) ).evaluate()( 0,0 );
        TEST_GEOTOOL2_CHECK_CLOSE( area3elt1, 2.0-2*0.6*0.4, 1e-9 );
    }
    auto area3elt2 = integrate( _range=elements( mesh3 ),
                                _expr=cst(1.) ).evaluate()( 0,0 );
    TEST_GEOTOOL2_CHECK_CLOSE( area3elt2, 4.0-2*0.6*0.4, 1e-9 );

    auto area3face = integrate( _range=boundaryfaces( mesh3 ),
                                _expr=cst(1.) ).evaluate()( 0,0 );
    TEST_GEOTOOL2_CHECK_CLOSE( area3face, 10.0 + 4*0.6 + 4*0.4, 1e-9 );
}


GeoTool::GeoGMSHTool
createCube(std::string nameObj,double l1,double meshSize,double centerX,double centerY,double centerZ, int type)
{
    GeoTool::Node x1( centerX-l1,centerY-l1,centerZ-l1);
    GeoTool::Node x2( centerX+l1,centerY-l1,centerZ-l1);
    GeoTool::Node x3( centerX+l1,centerY+l1,centerZ-l1);
    GeoTool::Node x4( centerX-l1,centerY+l1,centerZ-l1);
    GeoTool::Node x5( centerX-l1,centerY-l1,centerZ+l1);
    GeoTool::Node x6( centerX+l1,centerY-l1,centerZ+l1);
    GeoTool::Node x7( centerX+l1,centerY+l1,centerZ+l1);
    GeoTool::Node x8( centerX-l1,centerY+l1,centerZ+l1);
    GeoTool::Hexahedron H( meshSize,nameObj,x1,x2,x3,x4,x5,x6,x7,x8);

    if ( type == 1 )
    {
        H.setMarker(_type="surface",_name="Interface",_marker4=true);
        H.setMarker(_type="surface",_name="Inlet",_marker6=true);
        H.setMarker(_type="surface",_name="Wall",
                    _marker1=true,
                    _marker2=true,
                    _marker3=true,
                    _marker5=true);
        H.setMarker(_type="volume",_name="Omega1",_markerAll=true);
    }
    else
    {
        H.setMarker(_type="surface",_name="Interface",_marker6=true);
        H.setMarker(_type="surface",_name="Inlet",_marker4=true);
        H.setMarker(_type="surface",_name="Wall",
                    _marker1=true,
                    _marker2=true,
                    _marker3=true,
                    _marker5=true);
        H.setMarker(_type="volume",_name="Omega2",_markerAll=true);
    }
    return H;
}

void runFusion3d( bool keepInterface )
{
    double meshSize = doption(_name="hsize3d");
    typedef Mesh<Simplex<3,1,3> > mesh_type;
    std::string keepInterfaceStr( (keepInterface)? "with-interface" : "without-interface" );

    auto H1 = createCube("MyHex1",0.5,meshSize,0.0,0.0,0.0,1);
    auto H2 = createCube("MyHex2",0.5,meshSize/2.,1.0,0.0,0.0,2);

    auto mesh =
        (H1+H2).
        fusion(H1,4,H2,6,keepInterface).
        createMesh(_mesh=new mesh_type,
                   _name="domain3d_0_"+keepInterfaceStr );

}

} // namespace test_geotool2

#if defined(USE_BOOST_TEST)

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_geotool2::makeAbout(), test_geotool2::makeOptions() )

BOOST_AUTO_TEST_SUITE( interp_geotool2 )

BOOST_AUTO_TEST_CASE( interp_geotool2 )
{
    test_geotool2::runFusion2d( true );
    test_geotool2::runFusion2d( false );
    test_geotool2::runFusion3d( true );

}

BOOST_AUTO_TEST_SUITE_END()

#else

int main(int argc, char**argv )
{
    using namespace Feel;

	Environment env( _argc=argc, _argv=argv,
                     _desc=test_geotool2::makeOptions(),
                     _about=test_geotool2::makeAbout() );

    test_geotool2::runFusion2d( true );
    test_geotool2::runFusion2d( false );
    test_geotool2::runFusion3d( true );

    return 0;
}
#endif


