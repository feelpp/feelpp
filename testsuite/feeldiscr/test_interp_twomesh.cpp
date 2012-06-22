#define BOOST_TEST_MODULE interp_twomesh tests
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <feel/options.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
//#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <iostream>

using namespace Feel;
using namespace Feel::vf;



namespace test_interp_twomesh
{

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

enum type_champs
{
    ScalarTest = 0,
    VectorialTest
};


Application_ptrtype test_app;
/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_interp_twomesh options" );
    desc_options.add_options()
    ( "hsize1", po::value<double>()->default_value( 0.1 ), "mesh size1" ) //0.05
    ( "hsize2", po::value<double>()->default_value( 0.7 ), "mesh size2" ) //0.03
    ( "gmsh", po::value<float>()->default_value( 2.1 ), " version of gmsh(2.0 or 2.1)" )
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
    AboutData about( "Test_Interp_TwoMesh" ,
                     "Test_Interp_TwoMesh" ,
                     "0.1",


                     "Test1 verify geomap such as phi*ph^{-1}=Id ; Test2: on gamma, ||..v(u1)-..v(u2)|| < epsilon",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;

}

/*_________________________________________________*
 * run_test_geomap
 *_________________________________________________*/

template<uint32_type Dim,uint32_type OrderChamp,uint32_type OrderGeo>
void
run_test_geomap( Application_ptrtype & test_app )
{

    Log() << "[testGeoMap] start test geomap inverse\n";

    //-----------------------------------------------------------------------------------//

    typedef Mesh<Simplex<Dim,OrderGeo,Dim> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    //-----------------------------------------------------------------------------------//

    double meshSize1 = test_app->vm()["hsize1"].template as<double>();
    double meshSize2 = test_app->vm()["hsize2"].template as<double>();

    //-----------------------------------------------------------------------------------//

    //Geometry
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R1( meshSize1,"LEFT",x1,x2 );

    GeoTool::Node x_centre( 0.5,0.5 );
    GeoTool::Node x_bord( 0.8,0.5 );
    GeoTool::Circle C( meshSize2,"uncercle",x_centre,x_bord );

    //-----------------------------------------------------------------------------------//

    //Mesh
    //mesh_ptrtype mesh = createMesh<Dim,OrderGeo>( "mesh_geomap_" + mesh_type::shape_type::name(), (R1-C).geoStr() );
    mesh_ptrtype mesh = ( R1-C ).createMesh<mesh_type>( "mesh_geomap_" + mesh_type::shape_type::name() );

    //-----------------------------------------------------------------------------------//

    //geomap
    typedef typename mesh_type::gm_type gm_type;
    typedef typename boost::shared_ptr<gm_type> gm_ptrtype;
    gm_ptrtype __gm = mesh->gm();

    //geomap context
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename gm_type::template Context<vm::POINT, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;

    //-----------------------------------------------------------------------------------//

    //pts de Gauss
    typename _Q<OrderChamp>::template apply<mesh_type::nDim,typename mesh_type::value_type, Simplex >::type im;

    //-----------------------------------------------------------------------------------//

    // iterators sur les elements du maillage
    typename mesh_type::element_iterator el_it;
    typename mesh_type::element_iterator el_en;
    boost::tie( boost::tuples::ignore, el_it, el_en ) = Feel::elements( *mesh );

    typedef typename mesh_type::gm_type::precompute_type geopc_type;
    typedef typename mesh_type::gm_type::precompute_ptrtype geopc_ptrtype;
    geopc_ptrtype __geopc( new geopc_type( __gm, im.points() ) );

    for ( ; el_it!=el_en; ++el_it )
    {
        //init and apply the geometric transformation phi
        gmc_ptrtype __c( new gmc_type( __gm,
                                       *el_it,//mesh->element( *el_it ),
                                       __geopc ) );

        //init inverse geometric transformation phi^{-1}
        typename mesh_type::Inverse::gic_type gic( __gm, *el_it );

        for ( uint32_type i=0; i< __c->xReal().size2(); ++i )
        {

            // get phi for one point
            typename mesh_type::node_type n = ublas::column( __c->xReal(), i );

            //compute phi^{-1}
            gic.setXReal( n );

            //verify that : phi°phi^{-1}=Id
            for ( uint32_type j=0; j<n.size(); ++j )
            {
                double err = std::abs( ublas::column( im.points(),i )( j ) - gic.xRef()( j ) );
                //if ( err> 1e-9) std::cout<< "\nProb : "<< err;
                BOOST_CHECK( err<1e-9 );
            }

        }

    }

    //-----------------------------------------------------------------------------------//

    Log() << "[testGeoMap] finish test geomap inverse\n";

} //end run

/*_________________________________________________*
 * run_interp_boundary
 *_________________________________________________*/

template<uint32_type Dim,uint32_type OrderChamp,uint32_type OrderGeo>
void
test_interp_boundary( boost::tuple<
                      boost::shared_ptr< Mesh<Simplex<Dim,OrderGeo,Dim> > >,
                      boost::shared_ptr< Mesh<Simplex<Dim,OrderGeo,Dim> > > > __mesh,
                      boost::mpl::int_<ScalarTest>
                    )
{


    typedef Mesh<Simplex<Dim,OrderGeo,Dim> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    typedef bases<Lagrange<OrderChamp,Scalar,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef bases<Lagrange<OrderChamp,Vectorial,Continuous,PointSetFekete> > basis_v_type;
    typedef FunctionSpace<mesh_type, basis_v_type> space_v_type;
    typedef boost::shared_ptr<space_v_type> space_v_ptrtype;
    typedef typename space_v_type::element_type element_v_type;
    //BOOST_MPL_ASSERT_MSG( ( boost::is_same<mpl::int_<space_v_type::rank>, mpl::int_<1> > ), INVALID_RANK_SHOULD_BE_1_VECTORIAL, (mpl::int_<space_v_type::rank>, space_v_type, basis_v_type) );
    //-----------------------------------------------------------------------------------//

    mesh_ptrtype mesh1 = boost::get<0>( __mesh );
    mesh_ptrtype mesh2 = boost::get<1>( __mesh );

    space_ptrtype Xh1 = space_type::New( mesh1 );
    space_ptrtype Xh2 = space_type::New( mesh2 );
    space_v_ptrtype Yh1 = space_v_type::New( mesh1 );
    space_v_ptrtype Yh2 = space_v_type::New( mesh2 );

    element_type u1( Xh1, "u1" );
    element_type u2( Xh2, "u2" );
    element_v_type v1( Yh1, "v1" );
    element_v_type v2( Yh2, "v2" );

    //-----------------------------------------------------------------------------------//

    //AUTO (e ,exp(Px()*Py())*sin(2*M_PI*Px()));
    //AUTO (f , sin(2*M_PI*Px())*sin(2*M_PI*Py()));
    //AUTO (g , sin(0.5*M_PI*Px())*(cos(0.5*M_PI*Py())));
    //AUTO (g , Px()*(Px()-1)*Py()*(Py()-1) );
    AUTO ( h , ( Px()+1 )*( Px()-1 )*( Py()+1 )*( Py()-1 ) );
    //AUTO (h , cst(1.)*Px());

    u1 = vf::project( Xh1, elements( mesh1 ), h );
    u2 = vf::project( Xh2, elements( mesh2 ), h );
    v1 = vf::project( Yh1, elements( mesh1 ), vec( h,h ) );
    v2 = vf::project( Yh2, elements( mesh2 ), vec( h,h ) );

    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//

    double  __errId = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                            ( idv( u1 )-idv( u2 ) )*( idv( u1 )-idv( u2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] idv(u1) error : " << __errId << "\n";
    BOOST_CHECK_SMALL( __errId,1e-7 );
    //-----------------------------------------------------------------------------------//
    double  __errGrad = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                   ( gradv( u1 )-gradv( u2 ) )*trans( gradv( u1 )-gradv( u2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] gradv(u1) error : " << __errGrad <<"\n";
    BOOST_CHECK_SMALL( __errGrad,1e-5 );
    //-----------------------------------------------------------------------------------//
    double  __errDiv = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                  ( divv( v1 )-divv( v2 ) )*( divv( v1 )-divv( v2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] divv(u1) error : " << __errDiv << "\n";
    BOOST_CHECK_SMALL( __errDiv,1e-5 );
    //-----------------------------------------------------------------------------------//
    /*
    std::cout << "\ndiv 1 :" << std::sqrt(integrate( markedfaces(mesh1, mesh1->markerName("Interface")),
                                            (divv(v1))*(divv(v1)) ).evaluate()(0,0));
    std::cout << "\ndiv 2" << std::sqrt(integrate( markedfaces(mesh1, mesh1->markerName("Interface")),
                                                   print((divv(v2))*(divv(v2)),"div2") ).evaluate()(0,0));
    */
    //double  __errDiv = std::sqrt(integrate( markedfaces(mesh1, mesh1->markerName("Interface")),
    //                                      print((divv(v1)-divv(v2))*(divv(v1)-divv(v2)),"div=") ).evaluate()(0,0));
    //std::cout <<"\n" << std::sqrt(integrate( elements(mesh1),
    //                                         print((divv(v1))   ,"div=") ).evaluate()(0,0)) << "\n";
    //-----------------------------------------------------------------------------------//
    double  __errDivx = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                   ( dxv( v1 )+dyv( v1 )-dxv( v2 )-dyv( v2 ) )*( dxv( v1 )+dyv( v1 )-dxv( v2 )-dyv( v2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] dxv+dyv(v1) error : " << __errDivx << "\n";
    BOOST_CHECK_SMALL( __errDivx,1e-5 );
    //-----------------------------------------------------------------------------------//
    double  __errDx = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                            ( dxv( u1 )-dxv( u2 ) )*( dxv( u1 )-dxv( u2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] dxv(u1) error : " << __errDx << "\n";
    BOOST_CHECK_SMALL( __errDx,1e-5 );
    //-----------------------------------------------------------------------------------//
    double  __errDy = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                            ( dyv( u1 )-dyv( u2 ) )*( dyv( u1 )-dyv( u2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] dyv(u1) error : " << __errDy << "\n";
    BOOST_CHECK_SMALL( __errDy,1e-5 );
    //-----------------------------------------------------------------------------------//
#if 0
    double  __errDz = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                            ( dzv( u1 )-dzv( u2 ) )*( dzv( u1 )-dzv( u2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] dzv(u1) error : " << __errDz << "\n";
    BOOST_CHECK_SMALL( __errDz,1e-7 );
#endif
    //-----------------------------------------------------------------------------------//
    double  __errCurl = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                   trans( curlv( v1 )-curlv( v2 ) )*( curlv( v1 )-curlv( v2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] curlv(u1) error : " << __errCurl << "\n";
    BOOST_CHECK_SMALL( __errCurl,1e-7 );
    //-----------------------------------------------------------------------------------//
    double  __errCurlx = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                    ( curlxv( v1 )-curlxv( v2 ) )*( curlxv( v1 )-curlxv( v2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] curlxv(u1) error : " << __errCurlx << "\n";
    BOOST_CHECK_SMALL( __errCurlx,1e-5 );
    //-----------------------------------------------------------------------------------//
    double  __errCurly = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                    ( curlyv( v1 )-curlyv( v2 ) )*( curlyv( v1 )-curlyv( v2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] curlyv(u1) error : " << __errCurly << "\n";
    BOOST_CHECK_SMALL( __errCurly,1e-5 );
    //-----------------------------------------------------------------------------------//

#if 0
    double  __errCurlz = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                    ( curlzv( v1 )-curlzv( v2 ) )*( curlzv( v1 )-curlzv( v2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] curlzv(u1) error : " << __errCurlz << "\n";
    //-----------------------------------------------------------------------------------//

#endif
#if 0
    double  __errHess;

    if ( OrderGeo==1 )
        __errHess = std::sqrt( integrate( markedfaces( mesh1, mesh1->markerName( "Interface" ) ),
                                          trace( ( hessv( u1 )-hessv( u2 ) )*trans( hessv( u1 )-hessv( u2 ) ) ) ).evaluate()( 0,0 ) );

    Log() << "[testBoundary] hessv(u1) error : " << __errHess << "\n";
    BOOST_CHECK_SMALL( __errHess,1e-7 );
#endif
    //-----------------------------------------------------------------------------------//
    double  __errId2 = std::sqrt( integrate( markedfaces( mesh2, mesh2->markerName( "Interface" ) ),
                                  ( idv( u1 )-idv( u2 ) )*( idv( u1 )-idv( u2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] idv(u2) error : " << __errId2 << "\n";
    BOOST_CHECK_SMALL( __errId2,1e-7 );
    //-----------------------------------------------------------------------------------//
    double  __errGrad2 = std::sqrt( integrate( markedfaces( mesh2, mesh2->markerName( "Interface" ) ),
                                    ( gradv( u1 )-gradv( u2 ) )*trans( gradv( u1 )-gradv( u2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] gradv(u2) error : " << __errGrad2 << "\n";
    BOOST_CHECK_SMALL( __errGrad2,1e-5 );
    //-----------------------------------------------------------------------------------//
    double  __errDiv2 = std::sqrt( integrate( markedfaces( mesh2, mesh2->markerName( "Interface" ) ),
                                   ( divv( v1 )-divv( v2 ) )*( divv( v1 )-divv( v2 ) ) ).evaluate()( 0,0 ) );
    Log() << "[testBoundary] divv(u2) error : " << __errDiv2 << "\n";
    BOOST_CHECK_SMALL( __errDiv2,1e-5 );
    //-----------------------------------------------------------------------------------//

}

/*_________________________________________________*
 * run_test_interp
 *_________________________________________________*/

template<uint32_type Dim,uint32_type OrderChamp,uint32_type OrderGeo>
void
run_test_interp( Application_ptrtype & test_app )
{
    typedef Mesh<Simplex<Dim,OrderGeo,Dim> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    //-----------------------------------------------------------------------------------//

    double meshSize1 = test_app->vm()["hsize1"].template as<double>();
    double meshSize2 = test_app->vm()["hsize2"].template as<double>();

    //-----------------------------------------------------------------------------------//

    //Geometry
    GeoTool::Node x1( -1,-1 );
    GeoTool::Node x2( 0,1 );
    GeoTool::Rectangle R1( meshSize1,"LEFT",x1,x2 );
    R1.setMarker( _type="line",_name="Interface",_marker2=true );
    R1.setMarker( _type="line",_name="Bord",_marker1=true,_marker3=true,_marker4=true );
    R1.setMarker( _type="surface",_name="Omega1",_marker1=true );

    GeoTool::Node x3( 0,-1 );
    GeoTool::Node x4( 1,1 );
    GeoTool::Rectangle R2( meshSize2,"RIGHT",x3,x4 );
    R2.setMarker( _type="line",_name="Interface",_marker4=true );
    R2.setMarker( _type="line",_name="Bord",_marker1=true,_marker2=true,_marker3=true );
    R2.setMarker( _type="surface",_name="Omega2",_marker1=true );

    //-----------------------------------------------------------------------------------//

    GeoTool::Node x21( 0,0 );
    GeoTool::Node x22( 1,1 );
    GeoTool::Rectangle R( meshSize1,"RIGHT",x21,x22 );
    R.setMarker( _type="line",_name="Bord",_marker1=true,_marker2=true,_marker3=true,_marker4=true );
    R.setMarker( _type="surface",_name="Omega3",_marker1=true );

    GeoTool::Node x_centre( 0.5,0.5 );
    GeoTool::Node x_bord( 0.7,0.5 );
    GeoTool::Circle C( meshSize2,"DISQUE",x_centre,x_bord );
    C.setMarker( _type="line",_name="Interface",_marker1=true );
    //C.setMarker(_type="surface",_name="Omega2",_marker1=true);

    GeoTool::Node x2_centre( 0.5,0.5 );
    GeoTool::Node x2_bord( 0.7,0.5 );
    GeoTool::Circle C2( meshSize2,"DISQUE",x2_centre,x2_bord );
    C2.setMarker( _type="line",_name="Interface",_marker1=true );
    C2.setMarker( _type="surface",_name="Omega4",_marker1=true );

    //-----------------------------------------------------------------------------------//

    //Mesh
    mesh_ptrtype mesh1 = R1.createMesh<mesh_type>( "mesh1_interp" + mesh_type::shape_type::name() );
    mesh_ptrtype mesh2 = R2.createMesh<mesh_type>( "mesh2_interp" + mesh_type::shape_type::name() );
    mesh_ptrtype mesh3 = ( R-C ).createMesh<mesh_type>( "mesh3_interp" + mesh_type::shape_type::name() );
    mesh_ptrtype mesh4 = C2.createMesh<mesh_type>( "mesh4_interp" + mesh_type::shape_type::name() );

    std::list<boost::tuple<mesh_ptrtype,mesh_ptrtype> > __listMesh;
    __listMesh.push_back( boost::make_tuple( mesh1,mesh2 ) );
    __listMesh.push_back( boost::make_tuple( mesh3,mesh4 ) );

    //-----------------------------------------------------------------------------------//

    typename std::list<boost::tuple<mesh_ptrtype,mesh_ptrtype> >::iterator itMesh = __listMesh.begin();
    typename std::list<boost::tuple<mesh_ptrtype,mesh_ptrtype> >::iterator itMesh_end = __listMesh.end();

    int i=1;

    for ( ; itMesh!=itMesh_end ; ++itMesh, ++i )
    {
        Log() << "[run] start test interpolation mesh pair" << i << "\n";
        test_interp_boundary<Dim,OrderChamp,OrderGeo>( *itMesh, boost::mpl::int_<ScalarTest>() );
        Log() << "[run] finish test interpolation mesh pair" << i << "\n";
    }

}

/*_________________________________________________*
 * runTestExport
 *_________________________________________________*/

template<uint32_type Dim,uint32_type OrderChamp,uint32_type OrderGeo>
void
run_test_export( Application_ptrtype & test_app )
{

    Log() << "[testExport] starts\n";

    //-----------------------------------------------------------------------------------//

    typedef Mesh<Simplex<Dim,OrderGeo,Dim> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    typedef bases<Lagrange<OrderChamp,Scalar,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef Exporter<mesh_type> exporter_type;
    typedef boost::shared_ptr<exporter_type> exporter_ptrtype;

    //-----------------------------------------------------------------------------------//

    double meshSize1 = test_app->vm()["hsize1"].template as<double>();
    double meshSize2 = test_app->vm()["hsize2"].template as<double>();

    //-----------------------------------------------------------------------------------//

    //Geometry
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( meshSize1,"rectangle",x1,x2 );

    GeoTool::Node x_centre( 0.5,0.5 );
    GeoTool::Node x_bord( 0.7,0.5 );
    GeoTool::Circle C( meshSize2,"cercle",x_centre,x_bord );

    //-----------------------------------------------------------------------------------//

    //Mesh
    mesh_ptrtype mesh = ( R-C ).createMesh<mesh_type>( "mesh_test_export" );

    space_ptrtype Xh = space_type::New( mesh );

    element_type u( Xh, "u" );

    //AUTO (e ,exp(Px()*Py())*sin(2*M_PI*Px()));
    //AUTO (f , sin(2*M_PI*Px())*sin(2*M_PI*Py()));
    //AUTO (g , sin(0.5*M_PI*Px())*(cos(0.5*M_PI*Py())));
    //AUTO (g , Px()*(Px()-1)*Py()*(Py()-1) );
    AUTO ( h , ( Px()+1 )*( Px()-1 )*( Py()+1 )*( Py()-1 ) );
    //AUTO (h , cst(1.)*Px());

    u = vf::project( Xh, elements( mesh ), h );

    //-----------------------------------------------------------------------------------//

    exporter_ptrtype __exporter( exporter_type::New( "gmsh", test_app->about().appName()
                                 + "_"
                                 + mesh_type::shape_type::name() ) );

    __exporter->step( 0 )->setMesh( mesh );
    __exporter->step( 0 )->add( "u_scal", u );
    __exporter->save();
    Log() << "[testExport] finish\n";

} // run_test_export


}//end namespace test_interp_twomesh



/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * MAIN
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/



BOOST_AUTO_TEST_SUITE( interp_twomesh_testsuite )
Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                       boost::unit_test::framework::master_test_suite().argv );
BOOST_AUTO_TEST_CASE( interp_twomesh_geomap )
{

    using namespace test_interp_twomesh;

    test_app = Application_ptrtype( new Application_type( boost::unit_test::framework::master_test_suite().argc,
                                    boost::unit_test::framework::master_test_suite().argv,
                                    makeAbout(),
                                    makeOptions()
                                                        ) );

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/geomap/%1%/" )
                                % test_app->about().appName()
                              );

#if 1
    BOOST_MESSAGE(   "\n[main] ----------TEST_GEOMAP_START----------\n" );
    BOOST_MESSAGE(   "[main] ----------------<2,6,1>---------------\n" );
    run_test_geomap<2,6,1>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,6,2>---------------\n" );
    run_test_geomap<2,6,2>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,6,3>---------------\n" );
    run_test_geomap<2,6,3>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,6,4>---------------\n" );
    run_test_geomap<2,6,4>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,6,5>---------------\n" );
    run_test_geomap<2,6,5>( test_app );
    BOOST_MESSAGE(   "[main] ----------TEST_GEOMAP_FINISH----------\n\n" );
#endif
}

// problem with this test case
#if 1
BOOST_AUTO_TEST_CASE( interp_twomesh_interp )
{

    using namespace test_interp_twomesh;

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/interp/%1%/" )
                                % test_app->about().appName()
                              );

    BOOST_MESSAGE(   "[main] ----------TEST_INTERP_START-----------\n" );
    BOOST_MESSAGE(   "[main] ----------------<2,7,1>---------------\n" );
    run_test_interp<2,7,1>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,8,2>---------------\n" );
    run_test_interp<2,8,2>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,9,3>---------------\n" );
    run_test_interp<2,9,3>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,10,4>---------------\n" );
    run_test_interp<2,10,4>( test_app );
    //BOOST_MESSAGE(   "[main] ----------------<2,11,5>---------------\n");
    //run_test_interp<2,11,5>(test_app);
    BOOST_MESSAGE(   "[main] ----------TEST_INTERP_FINISH----------\n" );
}
#endif
#if 0
BOOST_AUTO_TEST_CASE( interp_twomesh_export )
{

    using namespace test_interp_twomesh;

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/export/%1%/" )
                                % test_app->about().appName()
                              );

    BOOST_MESSAGE(   "[main] ----------TEST_EXPORT_START-----------\n" );
    BOOST_MESSAGE(   "[main] ----------------<2,7,1>---------------\n" );
    run_test_export<2,7,1>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,8,2>---------------\n" );
    run_test_export<2,8,2>( test_app );
    BOOST_MESSAGE(   "[main] ----------------<2,9,3>---------------\n" );
    run_test_export<2,9,3>( test_app );
    BOOST_MESSAGE(   "[main] ----------TEST_EXPORT_FINISH----------\n" );


}
#endif
BOOST_AUTO_TEST_SUITE_END()
