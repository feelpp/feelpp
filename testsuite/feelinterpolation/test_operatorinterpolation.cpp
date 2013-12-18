
#define BOOST_TEST_MODULE test_operatorinterpolation
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>


using namespace Feel;

namespace test_operatorinterpolation
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
    po::options_description desc_options( "test_operatorinterpolation options" );
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
    AboutData about( "Test_OperatorInterpolation" ,
                     "Test_OperatorInterpolation" ,
                     "0.1",
                     "test composant of a field vectorial",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}


template <uint32_type OrderGeo>
void
test2dTo1d()
{

    typedef Backend<double> backend_type;

    // typedef Mesh<Simplex<1,OrderGeo,2> > mesh_1d_type;
    // typedef Mesh<Simplex<2,OrderGeo,2> > mesh_2d_type;

    typedef Mesh<Hypercube<1,OrderGeo,2> > mesh_1d_type;
    typedef Mesh<Hypercube<2,OrderGeo,2> > mesh_2d_type;

    //typedef Mesh<Simplex<1,OrderGeo,2> > mesh_1d_type;
    //typedef Mesh<Simplex<2,OrderGeo,2> > mesh_2d_type;


    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_2d_type;
    typedef FunctionSpace<mesh_2d_type, basis_2d_type> space_2d_type;
    typedef typename space_2d_type::element_type element_2d_type;

    typedef bases<Lagrange<5,Vectorial,Continuous,PointSetFekete> > basis_1d_type;
    typedef FunctionSpace<mesh_1d_type, basis_1d_type> space_1d_type;
    typedef typename space_1d_type::element_type element_1d_type;

    //-----------------------------------------------------------//

    auto meshSize = option(_name="hsize").template as<double>();
    BOOST_TEST_MESSAGE( "meshSize=" << meshSize );

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 2,1 );
    GeoTool::Rectangle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh2d = C.createMesh( _mesh=new mesh_2d_type,
                                _name="test2dTo1d_domain"+mesh_2d_type::shape_type::name() );

#if 0
    GeoTool::Node x3( 0,0 );
    GeoTool::Node x4( 2,1 );
    GeoTool::Line L( meshSize, "Line",x3,x4 );
    L.setMarker( _type="point",_name="Sortie",_markerAll=true );
    L.setMarker( _type="line",_name="Omega1d",_markerAll=true );
    auto mesh1d = L.createMesh( _mesh=new mesh_1d_type,
                                _name="test2dTo1d_domain1d"+mesh_1d_type::shape_type::name(),_straighten=0 );
#else
    // create the 1d mesh from a submesh of 2d mesh because gmsh don't partion 1d mesh
    GeoTool::Node x3( 0,0 );
    GeoTool::Node x4( 2,1 );
    GeoTool::Node x5( 0,1 );
    GeoTool::Triangle T( meshSize,"OMEGA",x3,x4,x5);
    T.setMarker( _type="line",_name="Boundary1",_marker1=true );
    T.setMarker( _type="line",_name="Boundary2",_marker2=true,_marker3=true );
    T.setMarker( _type="surface",_name="Omega2dFor1d",_markerAll=true );
    auto mesh2dFor1d = T.createMesh( _mesh=new mesh_2d_type,
                                     _name="test2dTo1d_domain2dFor1d"+mesh_2d_type::shape_type::name(), _straighten=0 );
    auto mesh1d = createSubmesh(mesh2dFor1d,markedfaces(mesh2dFor1d,"Boundary1"));
#endif

    //-----------------------------------------------------------//

    auto Xh1d = space_1d_type::New( _mesh=mesh1d );
    auto Xh2d = space_2d_type::New( _mesh=mesh2d );

    auto exprProj = vec( cos( M_PI*Px() ),sin( M_PI*Py() ) );
    //auto exprProj = vec( 2*Px()*Py(), 0.25*(1.0-Py())*Px() );

    auto u1dinterp = Xh1d->element();
    auto u2dproj = vf::project( _space=Xh2d,
                                _range=elements( mesh2d ),
                                _expr=exprProj );
    auto u1dproj = vf::project( _space=Xh1d,
                                _range=elements( mesh1d ),
                                _expr=exprProj );

    //-----------------------------------------------------------//

    auto opI=opInterpolation( _domainSpace=Xh2d,
                              _imageSpace=Xh1d );
    opI->apply( u2dproj,u1dinterp );

    auto s = integrate( _range=elements( mesh1d ),
                        _expr=inner( idv( u1dinterp )-idv( u1dproj ), idv( u1dinterp )-idv( u1dproj ) ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s,1e-8 );
    BOOST_TEST_MESSAGE( "s=" << s );

} // test2dTo1d

//---------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------//

template <uint32_type OrderGeo>
void
test2dTo2d()
{
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;

    //typedef bases<Lagrange<3,Scalar,Continuous,PointSetFekete> > basis_1_type;
    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_1_type;
    typedef FunctionSpace<mesh_type, basis_1_type> space_1_type;
    typedef typename space_1_type::element_type element_1_type;

    //typedef bases<Lagrange<4,Scalar,Continuous,PointSetFekete> > basis_2_type;
    typedef bases<Lagrange<4,Vectorial,Continuous,PointSetFekete> > basis_2_type;
    typedef FunctionSpace<mesh_type, basis_2_type> space_2_type;
    typedef typename space_2_type::element_type element_2_type;

    //-------------------------------------------------------
    //case 1 : same mesh
    //-------------------------------------------------------
    WorldComm myWorldComm = Environment::worldComm();
    auto meshSize = option(_name="hsize").template as<double>();
    LOG(INFO) << "meshSize=" << meshSize << "\n";
    BOOST_TEST_MESSAGE( "meshSize=" << meshSize );
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.6,0 );
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name="test2dTo2d_domain"+mesh_type::shape_type::name(),
                              _partitions=myWorldComm.localSize() );

    auto Xh1 = space_1_type::New( _mesh=mesh );
    auto Xh2 = space_2_type::New( _mesh=mesh );
    auto u1 = Xh1->element();
    auto u2 = Xh2->element();
    auto u2a = Xh2->element();

    auto exprProj = vec( cos( M_PI*Px() ),sin( M_PI*Py() ) );
    u1 = vf::project( _space=Xh1,
                      _range=elements( mesh ),
                      _expr=exprProj );

    auto mybackend = backend(_rebuild=true);

    auto opI=opInterpolation( _domainSpace=Xh1,
                              _imageSpace=Xh2,
                              _backend=mybackend );
    opI->apply( u1,u2 );

    auto s1 = integrate( _range=elements( mesh ),
                         _expr=trans( idv( u1 ) )*idv( u1 ) ).evaluate()( 0,0 );
    auto s2 = integrate( _range=elements( mesh ),
                         _expr=trans( idv( u2 ) )*idv( u2 ) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( s1, s2,1e-8 );
    BOOST_TEST_MESSAGE( "s1=" << s1 << " s2=" << s2 << " s1-s2 = " << s1-s2 );

    auto opIa=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2,
                               _range=boundaryfaces( mesh ),
                               _backend=mybackend );
    opIa->apply( u1,u2a );

    auto s2a = integrate( _range=boundaryfaces( mesh ),
                          _expr=trans( idv( u1 )-idv( u2a ) )*( idv( u1 )-idv( u2a ) ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s2a,1e-8 );
    BOOST_TEST_MESSAGE( "s2a=" << s2a );

    //-------------------------------------------------------//
    //case 2 : with interpolation tool
    //-------------------------------------------------------//
    GeoTool::Node x3(0.4,0);
    GeoTool::Circle C2( meshSize/2.,"OMEGA",x1,x3);
    C2.setMarker(_type="line",_name="Boundary",_markerAll=true);
    C2.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh2 = C2.createMesh(_mesh=new mesh_type,
                               _name="test2dTo2d_domain2"+mesh_type::shape_type::name(),
                               _partitions=myWorldComm.localSize() );
    auto Xh2bis = space_2_type::New(_mesh=mesh2);
    auto u2bis = Xh2bis->element();
    auto u2bisbis = Xh2bis->element();
    auto u2bisproj = vf::project( _space=Xh2bis,
                                  _range=elements( mesh2 ),
                                  _expr=exprProj );
    // opInterp on elements
    auto opI2=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=elements( mesh2 ),
                               _backend=mybackend );
    opI2->apply( u1,u2bis );
    auto s3 = integrate( _range=elements( mesh2 ),
                         _expr=inner( idv( u2bis )-idv( u2bisproj ), idv( u2bis )-idv( u2bisproj ) )
                         /*_geomap=GeomapStrategyType::GEOMAP_HO*/ ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s3,1e-8 );
    BOOST_TEST_MESSAGE( "s3=" << s3 );

    // opInterp on faces
    auto opI3=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=boundaryfaces( mesh2 ),
                               _backend=mybackend );
    opI3->apply( u1,u2bisbis );
    auto s4 = integrate( _range=boundaryfaces( mesh2 ),
                         _expr=inner( idv( u2bisbis )-idv( u2bisproj ), idv( u2bisbis )-idv( u2bisproj ) ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s4,1e-8 );
    BOOST_TEST_MESSAGE( "s4=" << s4 );

} // test2dTo2d


#if defined(FEELPP_HAS_VTK)
template <uint32_type OrderGeo>
void
test2dOpLagrangeP1()
{
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;
    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > basis_P1_type;
    typedef FunctionSpace<mesh_type, basis_P1_type> space_P1_type;

    WorldComm myWorldComm;
    auto meshSize = option(_name="hsize").template as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.6,0 );
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name="test2dOpLagrangeP1_domain"+mesh_type::shape_type::name(),
                              _partitions=myWorldComm.localSize() );

    auto Xh = space_type::New( _mesh=mesh );
    auto exprProj = vec( cos( M_PI*Px() ),sin( M_PI*Py() ) );
    auto u = vf::project( _space=Xh,
                          _range=elements( mesh ),
                          _expr=exprProj );

    auto mybackend = backend(_rebuild=true);

    auto opLagP1 = lagrangeP1(_space=Xh);
    auto meshLagP1 = opLagP1->mesh();

    auto XhLagP1 = space_P1_type::New( _mesh=meshLagP1 );
    auto uLagP1interp = XhLagP1->element();

    auto uLagP1proj = vf::project( _space=XhLagP1,
                          _range=elements( meshLagP1 ),
                          _expr=exprProj );

    auto opI=opInterpolation( _domainSpace=Xh,
                              _imageSpace=XhLagP1,
                              _range=elements( meshLagP1 ) );
    opI->apply( u,uLagP1interp );

    auto s1 = integrate(_range=elements(meshLagP1),
                        _expr=inner( idv(uLagP1interp)-idv(uLagP1proj) , idv(uLagP1interp)-idv(uLagP1proj) ) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s1,1e-6);

#if 0
    auto myexporter = exporter( _mesh=meshLagP1, _name="test2dOpLagrangeP1_MyExport" );
    myexporter->add( "test2dOpLagrangeP1_uLagP1", uLagP1 );
    myexporter->save();
#endif
}


#endif // FEELPP_HAS_VTK


} // namespace test_operatorinterpolation

/**
 * main code
 */
FEELPP_ENVIRONMENT_WITH_OPTIONS( test_operatorinterpolation::makeAbout(),
                                 test_operatorinterpolation::makeOptions() )

BOOST_AUTO_TEST_SUITE( interp_operatorinterpolation )

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_2d_geo1 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_2d_geo1" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/2d2dgeo1" )
                                   % Environment::about().appName() );

    test_operatorinterpolation::test2dTo2d<1>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_2d_geo1 done" );
}

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_2d_geo2 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_2d_geo2" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/2d2dgeo2" )
                                   % Environment::about().appName() );

    test_operatorinterpolation::test2dTo2d<2>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_2d_geo2 done" );
}

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_1d_geo1 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_1d_geo1" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/2d1dgeo1" )
                                   % Environment::about().appName() );

    //test_operatorinterpolation::test2dTo2d<3>(test_app);
    test_operatorinterpolation::test2dTo1d<1>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_1d_geo1 done" );
}

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_1d_geo2 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_1d_geo2" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/2d1dgeo2" )
                                   % Environment::about().appName() );

    test_operatorinterpolation::test2dTo1d<2>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_2d_1d_geo2 done" );
}

#if defined(FEELPP_HAS_VTK)
BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_oplagp1_geo1 )
{
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_oplagp1_geo1" );
    using namespace test_operatorinterpolation;

    Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/oplagp1geo1" )
                                   % Environment::about().appName() );


    test_operatorinterpolation::test2dOpLagrangeP1<1>();
    BOOST_TEST_MESSAGE( "interp_operatorinterpolation_oplagp1_geo1 done" );
}
#endif // FEELPP_HAS_VTK


BOOST_AUTO_TEST_SUITE_END()

