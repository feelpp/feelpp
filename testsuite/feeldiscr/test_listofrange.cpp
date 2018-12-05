#define BOOST_TEST_MODULE test_listofrange
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( listofrange )

BOOST_AUTO_TEST_CASE( listofrange1 )
{
    typedef Mesh<Simplex<2,1,2> > mesh_type;
    GeoTool::Rectangle R( doption(_name="gmsh.hsize"),"OMEGA1",
                          GeoTool::Node(0,0),GeoTool::Node(2,1) );
    R.setMarker(_type="line",_name="Boundary1",_marker1=true);
    R.setMarker(_type="line",_name="Boundary2",_marker2=true);
    R.setMarker(_type="line",_name="Boundary3",_marker3=true);
    R.setMarker(_type="line",_name="Boundary4",_marker4=true);
    R.setMarker(_type="surface",_name="Omega1",_markerAll=true);
    GeoTool::Rectangle R2( doption(_name="gmsh.hsize"),"OMEGA2",
                          GeoTool::Node(3,3),GeoTool::Node(4,5) );
    R2.setMarker(_type="line",_name="BoundaryOthers",_markerAll=true);
    R2.setMarker(_type="surface",_name="Omega2",_markerAll=true);
    auto mesh = (R+R2).createMesh(_mesh=new mesh_type,_name= "domain" );

    std::list<std::string> myElementMarkers = { "Omega1","Omega2" };
    std::list<std::string> myFaceMarkers = { "Boundary1","Boundary2","Boundary3","Boundary4","BoundaryOthers" };
    double measSurfaceRef = 4.;
    double measBoundaryRef = 12.;

    // integrate evaluated
    double measSurface = integrate( _range=markedelements( mesh,myElementMarkers ), _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( measSurface-measSurfaceRef,1e-12 );
    double measBoundary = integrate( _range=markedfaces( mesh,myFaceMarkers ), _expr=cst(1.) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( measBoundary-measBoundaryRef,1e-12 );

    // create submesh
    auto submesh1 = createSubmesh( mesh, markedelements( mesh,myElementMarkers ) );
    BOOST_CHECK_SMALL( submesh1->measure()-measSurfaceRef,1e-12 );
    auto submesh2 = createSubmesh( mesh, markedfaces( mesh,myFaceMarkers ) );
    BOOST_CHECK_SMALL( submesh2->measure()-measBoundaryRef,1e-12 );

    // projection
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    u.on(_range=markedelements( mesh,myElementMarkers ),_expr=cst(3.) );
    double intu = integrate( _range=markedelements( mesh,myElementMarkers ), _expr=idv(u) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( intu-3*measSurfaceRef,1e-12 );
    auto v = Vh->element();
    v.on(_range=markedfaces( mesh,myFaceMarkers ),_expr=cst(4.) );
    double intv = integrate( _range=markedfaces( mesh,myFaceMarkers ), _expr=idv(v) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( intv-4*measBoundaryRef,1e-12 );

    // bilinear/linear form
    auto l1 = form1( _test=Vh );
    l1 = integrate(_range=elements( mesh ),
                   _expr=id(u));
    l1+= integrate(_range=boundaryfaces( mesh ),
                   _expr=id(u));
    auto a1 = form2( _trial=Vh, _test=Vh);
    a1 = integrate(_range=elements( mesh ),
                   _expr=gradt(u)*trans(grad(u)) );
    a1+= integrate(_range=boundaryfaces( mesh ),
                   _expr=inner(gradt(u)*N(),id(u)) );
    a1+=on(_range=boundaryfaces( mesh ), _rhs=l1, _element=u, _expr=cst(0.) );
    a1.matrixPtr()->close();

    auto l2 = form1( _test=Vh );
    l2 = integrate(_range=markedelements( mesh,myElementMarkers ),
                  _expr=id(u));
    l2+= integrate(_range=markedfaces( mesh,myFaceMarkers ),
                  _expr=id(u));
    auto a2 = form2( _trial=Vh, _test=Vh);
    a2 = integrate(_range=markedelements( mesh,myElementMarkers ),
                  _expr=gradt(u)*trans(grad(u)) );
    a2+= integrate(_range=markedfaces( mesh,myFaceMarkers ),
                   _expr=inner(gradt(u)*N(),id(u)) );
    a2+=on(_range=markedfaces( mesh,myFaceMarkers ), _rhs=l2, _element=u, _expr=cst(0.) );
    a2.matrixPtr()->close();

    double evall1Norm1 = a1.matrixPtr()->l1Norm();
    double evallinftyNorm1 = a1.matrixPtr()->linftyNorm();
    double evall1Norm2 = a2.matrixPtr()->l1Norm();
    double evallinftyNorm2 = a2.matrixPtr()->linftyNorm();
    BOOST_CHECK_SMALL( evall1Norm1-evall1Norm2,1e-12 );
    BOOST_CHECK_SMALL( evallinftyNorm1-evallinftyNorm2,1e-12 );
    a2.matrixPtr()->addMatrix(-1.0, a1.matrixPtr() );
    evall1Norm2 = a2.matrixPtr()->l1Norm();
    evallinftyNorm2 = a2.matrixPtr()->linftyNorm();
    BOOST_CHECK_SMALL( evall1Norm2, 1e-12 );
    BOOST_CHECK_SMALL( evallinftyNorm2, 1e-12 );

    // interpolation operator
    auto Wh = Pch<3>( mesh );
    auto w = Wh->element();
    auto opI = opInterpolation( _domainSpace=Vh,_imageSpace=Wh,
                                _range=markedfaces( mesh,myFaceMarkers ) );
    opI->apply( v,w );
    double intw = integrate( _range=markedfaces( mesh,myFaceMarkers ), _expr=idv(w) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( intw-4*measBoundaryRef,1e-12 );

    auto VhSub = Pch<1>( submesh1 );
    auto uSub = VhSub->element( cst(5.) );
    auto opI2 = opInterpolation( _domainSpace=VhSub,_imageSpace=Vh,
                                 _range=markedelements( mesh,myElementMarkers ) );
    opI2->apply( uSub,u );
    intu = integrate( _range=markedelements( mesh,myElementMarkers ), _expr=idv(u) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( intu-5*measSurfaceRef,1e-12 );

    auto WhSub = Pch<1>( submesh2 );
    auto wSub = WhSub->element( cst(6.) );
    auto opI3 = opInterpolation( _domainSpace=WhSub,_imageSpace=Vh,
                                 _range=markedfaces( mesh,myFaceMarkers ) );
    opI3->apply( wSub,u );
    intw = integrate( _range=markedfaces( mesh,myFaceMarkers ), _expr=idv(u) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( intw-6*measBoundaryRef,1e-12 );

}


BOOST_AUTO_TEST_SUITE_END()



