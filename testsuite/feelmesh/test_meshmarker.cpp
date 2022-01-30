
#define BOOST_TEST_MODULE test_meshmarker
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/geotool.hpp>
#include <feel/feelvf/vf.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_meshmarker )

BOOST_AUTO_TEST_CASE( test_meshmarker1 )
{
    using namespace Feel;
    typedef Mesh<Simplex<3,1,3> > mesh_type;
    double l = 0.3;
    GeoTool::Node x1(-l,-l,-l);
    GeoTool::Node x2( l,-l,-l);
    GeoTool::Node x3( l, l,-l);
    GeoTool::Node x4(-l, l,-l);
    GeoTool::Node x5(-l,-l, l);
    GeoTool::Node x6( l,-l, l);
    GeoTool::Node x7( l, l, l);
    GeoTool::Node x8(-l, l, l);
    GeoTool::Hexahedron H(doption(_name="gmsh.hsize"),"UnHexa",x1,x2,x3,x4,x5,x6,x7,x8);
    H.setMarker(_type="point",_name="markPoints",_markerAll=true);
    H.setMarker(_type="line",_name="markLines",_markerAll=true);
    H.setMarker(_type="surface",_name="GammaDirichlet",_marker2=true,_marker3=true,_marker4=true,_marker5=true,_marker6=true);
    H.setMarker(_type="surface",_name="GammaNeumann",_marker1=true);
    H.setMarker(_type="volume",_name="OmegaFluid",_markerAll=true);

    auto mesh = H.createMesh(_mesh= new mesh_type,
                             _name="un_cube" );

    double intSurf1 = integrate(_range=markedfaces(mesh, "GammaDirichlet" ),_expr=cst(1.) ).evaluate()(0,0);
    double intSurf2 = integrate(_range=markedfaces(mesh, "GammaNeumann" ),_expr=cst(1.) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intSurf1,5*(2*l*2*l),1e-12 );
    BOOST_CHECK_CLOSE( intSurf2,1*(2*l*2*l),1e-12 );

    auto submeshFaces = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,"GammaNeumann"));
    double intSubmeshFaces = integrate(_range=elements(submeshFaces),_expr=cst(1.) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intSubmeshFaces,1*(2*l*2*l),1e-12 );
    double intMarkedLinesSubmeshFaces = integrate(_range=markedfaces(submeshFaces,"markLines"),_expr=cst(1.) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intMarkedLinesSubmeshFaces,4*2*l,1e-12 );
    size_type nMarkedPointsSubmeshFaces = nelements( markedpoints(submeshFaces,"markPoints"), true );
    if ( Environment::numberOfProcessors() > 1 )
        BOOST_CHECK( nMarkedPointsSubmeshFaces >= 4 );
    else
        BOOST_CHECK( nMarkedPointsSubmeshFaces == 4 );

    auto submeshLines = createSubmesh(_mesh=mesh,_range=markededges(mesh,"markLines"),_update=size_type(0));
    double intSubmeshLines = integrate(_range=elements(submeshLines),_expr=cst(1.) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intSubmeshLines,12*2*l,1e-12 );
    size_type nMarkedPointsSubmeshLines = nelements( markedpoints(submeshLines,"markPoints"), true );
    if ( Environment::numberOfProcessors() > 1 )
        BOOST_CHECK( nMarkedPointsSubmeshLines >= 8 );
    else
        BOOST_CHECK( nMarkedPointsSubmeshLines == 8 );

    size_type nMarkedPoints = nelements( markedpoints(mesh,"markPoints"), true );
    if ( Environment::numberOfProcessors() > 1 )
        BOOST_CHECK( nMarkedPoints >= 8 );
    else
        BOOST_CHECK( nMarkedPoints == 8 );

    mesh->updateMarker3WithRange(markedfaces(mesh,"GammaNeumann"),1);
    mesh->updateMarkersFromFaces();
    auto submesh = createSubmesh(_mesh=mesh,_range=marked3elements(mesh,1));
    double intSubmesh = integrate(_range=markedfaces(submesh,"GammaNeumann"),_expr=cst(1.) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intSubmesh,1*(2*l*2*l),1e-12 );
    size_type nMarkedPointsSubmesh = nelements( markedpoints(submesh,"markPoints"), true );
    if ( Environment::numberOfProcessors() > 1 )
        BOOST_CHECK( nMarkedPointsSubmesh >= 4 );
    else
        BOOST_CHECK( nMarkedPointsSubmesh == 4 );

}

BOOST_AUTO_TEST_SUITE_END()

