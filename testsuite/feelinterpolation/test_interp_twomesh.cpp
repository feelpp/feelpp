
#if 1
#define BOOST_TEST_MODULE interp_twomesh tests
#include <testsuite/testsuite.hpp>
#endif

#include <feel/options.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <iostream>

using namespace Feel;
using namespace Feel::vf;



namespace test_interp_twomesh
{
  /*
typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

enum type_champs
{
    ScalarTest = 0,
    VectorialTest
};
  */

  //Application_ptrtype test_app;
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
    ( "hsize2", po::value<double>()->default_value( 0.3 ), "mesh size2" ) //0.03
    ( "userelation", po::value<bool>()->default_value( false ), "mesh size2" ) //0.03
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
    AboutData about( "test_interp_twomesh" ,
                     "test_interp_twomesh" ,
                     "0.1",
                     "Test1 verify geomap such as phi*ph^{-1}=Id ; Test2: on gamma, ||..v(u1)-..v(u2)|| < epsilon",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@feelpp.org", "" );
    return about;

}

/*_________________________________________________*
 * run_test_geomap
 *_________________________________________________*/

template<uint32_type Dim,uint32_type OrderChamp,uint32_type OrderGeo>
void
run_test_geomap()
{

    LOG(INFO) << "[testGeoMap] start test geomap inverse\n";

    //-----------------------------------------------------------------------------------//

    typedef Mesh<Simplex<Dim,OrderGeo,Dim> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    //-----------------------------------------------------------------------------------//

    //Geometry
    double meshSize1 = doption(_name="hsize1");
    double meshSize2 = doption(_name="hsize2");
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R1( meshSize1,"LEFT",x1,x2 );
    GeoTool::Node x_centre( 0.5,0.5 );
    GeoTool::Node x_bord( 0.8,0.5 );
    GeoTool::Circle C( meshSize2,"uncercle",x_centre,x_bord );
    //Mesh
    mesh_ptrtype mesh = ( R1-C ).createMesh(_mesh=new mesh_type,
                                            _name= "mesh_geomap_" + mesh_type::shape_type::name() );

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

    //init inverse geometric transformation phi^{-1}
    typename mesh_type::Inverse::gic_type gic( __gm, *el_it );

    typedef typename mesh_type::gm_type::precompute_type geopc_type;
    typedef typename mesh_type::gm_type::precompute_ptrtype geopc_ptrtype;
    geopc_ptrtype __geopc( new geopc_type( __gm, im.points() ) );

    for ( ; el_it!=el_en; ++el_it )
    {
        //init and apply the geometric transformation phi
        gmc_ptrtype __c( new gmc_type( __gm,
                                       *el_it,//mesh->element( *el_it ),
                                       __geopc ) );
#if 0
        //init inverse geometric transformation phi^{-1}
        typename mesh_type::Inverse::gic_type gic( __gm, *el_it );
#else
        gic.update( *el_it );
#endif
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
#if USE_BOOST_TEST
                BOOST_CHECK( err<1e-9 );
#endif
            }

        }

    }

    //-----------------------------------------------------------------------------------//

    LOG(INFO) << "[testGeoMap] finish test geomap inverse\n";

} //end run

/*_________________________________________________*
 * run_interp_boundary
 *_________________________________________________*/

template<uint32_type Dim,uint32_type OrderChamp,uint32_type OrderGeo>
void
test_interp_boundary( boost::tuple<
                      boost::shared_ptr< Mesh<Simplex<Dim,OrderGeo,Dim> > >,
                      boost::shared_ptr< Mesh<Simplex<Dim,OrderGeo,Dim> > >,
                      bool> __mesh/*,
                                    boost::mpl::int_<ScalarTest>*/
                    )
{
    auto mesh1 = boost::get<0>( __mesh );
    auto mesh2 = boost::get<1>( __mesh );
    auto Xh1 = Pch<OrderChamp>( mesh1 );
    auto Xh2 = Pch<OrderChamp>( mesh2 );
    auto Yh1 = Pchv<OrderChamp>( mesh1 );
    auto Yh2 = Pchv<OrderChamp>( mesh2 );
    //-----------------------------------------------------------------------------------//

    //AUTO (e ,exp(Px()*Py())*sin(2*M_PI*Px()));
    //AUTO (f , sin(2*M_PI*Px())*sin(2*M_PI*Py()));
    //AUTO (g , sin(0.5*M_PI*Px())*(cos(0.5*M_PI*Py())));
    //AUTO (g , Px()*(Px()-1)*Py()*(Py()-1) );
    //AUTO ( h , ( Px()+1 )*( Px()-1 )*( Py()+1 )*( Py()-1 ) );

    //auto h = ( Px()+1 )*( Px()-1 )*( Py()+1 )*( Py()-1 );
    //auto h = Px()*Py();
    auto theexpr = Px()*Px()*Py()*Py();

    auto u1 = vf::project( _space=Xh1, _range=elements( mesh1 ), _expr=theexpr );
    auto u2 = vf::project( _space=Xh2, _range=elements( mesh2 ), _expr=theexpr );
    auto v1 = vf::project( _space=Yh1, _range=elements( mesh1 ), _expr=theexpr*one() );
    auto v2 = vf::project( _space=Yh2, _range=elements( mesh2 ), _expr=theexpr*one() );

    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    bool useConformalIntegration=boost::get<2>( __mesh );//true;//false;
    double  __errId = normL2( _range=markedfaces( mesh1, "Interface" ),
                              _expr=idv( u1 )-idv( u2, useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] idv(u1) error : " << __errId );
    BOOST_CHECK_SMALL( __errId,1e-12 );
#else
    LOG(INFO) << "[testBoundary] idv(u1) error : " << __errId << "\n";
    std::cout << "[testBoundary] idv(u1) error : " << __errId << "\n";
#endif
    //-----------------------------------------------------------------------------------//
    double __errGrad = normL2( _range=markedfaces( mesh1, "Interface" ),
                               _expr= gradv( u1 )-gradv( u2 , useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] gradv(u1) error : " << __errGrad );
    BOOST_CHECK_SMALL( __errGrad,1e-12 );
#else
    LOG(INFO) << "[testBoundary] gradv(u1) error : " << __errGrad <<"\n";
    std::cout << "[testBoundary] gradv(u1) error : " << __errGrad <<"\n";
#endif
    //-----------------------------------------------------------------------------------//

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
    //-----------------------------------------------------------------------------------//
    double  __errDx = normL2( _range=markedfaces( mesh1, "Interface" ),
                              _expr=dxv( u1 )-dxv( u2, useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] dxv(u1) error : " << __errDx );
    BOOST_CHECK_SMALL( __errDx,1e-12 );
#else
    LOG(INFO) << "[testBoundary] dxv(u1) error : " << __errDx << "\n";
    std::cout << "[testBoundary] dxv(u1) error : " << __errDx << "\n";
#endif
    //-----------------------------------------------------------------------------------//
    double  __errDy = normL2( _range=markedfaces( mesh1, "Interface" ),
                              _expr=dyv( u1 )-dyv( u2, useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] dyv(u1) error : " << __errDy );
    BOOST_CHECK_SMALL( __errDy,1e-12 );
#else
    LOG(INFO) << "[testBoundary] dyv(u1) error : " << __errDy << "\n";
    std::cout << "[testBoundary] dyv(u1) error : " << __errDy << "\n";
#endif
    //-----------------------------------------------------------------------------------//
    if ( OrderGeo==1 )
      {
        double __errHess = normL2( _range=markedfaces( mesh1, "Interface" ),
                            _expr= hessv( u1 )-hessv( u2, useConformalIntegration ) );
#if USE_BOOST_TEST
        BOOST_MESSAGE( "[testBoundary] hessv(u1) error : " << __errHess );
        BOOST_CHECK_SMALL( __errHess,1e-9 );
#else
        std::cout << "[testBoundary] hessv(u1) error : " << __errHess << "\n";
#endif
      }
    //-----------------------------------------------------------------------------------//

    double  __errIdVec = normL2( _range=markedfaces( mesh1, "Interface" ),
                                 _expr=idv( v1 )-idv( v2, useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] idv vectorial error : " << __errIdVec );
    BOOST_CHECK_SMALL( __errIdVec,1e-12 );
#else
    std::cout << "[testBoundary] idv vectorial error : " << __errIdVec <<"\n";
#endif
    //-----------------------------------------------------------------------------------//

    double __errGradVec = normL2( _range=markedfaces( mesh1, "Interface" ),
                                  _expr= gradv( v1 )-gradv( v2 , useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] gradv vectorial error : " << __errGradVec );
    BOOST_CHECK_SMALL( __errGradVec,1e-12 );
#else
    std::cout << "[testBoundary] gradv vectorial error : " << __errGradVec <<"\n";
#endif
    //-----------------------------------------------------------------------------------//
    double  __errDiv = normL2( _range=markedfaces( mesh1, "Interface" ),
                               _expr=divv( v1 )-divv( v2 , useConformalIntegration) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] divv(u1) error : " << __errDiv );
    BOOST_CHECK_SMALL( __errDiv,1e-12 );
#else
    LOG(INFO) << "[testBoundary] divv(u1) error : " << __errDiv << "\n";
    std::cout << "[testBoundary] divv(u1) error : " << __errDiv << "\n";
#endif
    //-----------------------------------------------------------------------------------//
    double  __errDivx = normL2( _range=markedfaces( mesh1, "Interface" ),
                                _expr=dxv( v1 )+dyv( v1 )-dxv( v2, useConformalIntegration )-dyv( v2, useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] dxv+dyv(v1) error : " << __errDivx );
    BOOST_CHECK_SMALL( __errDivx,1e-12 );
#else
    LOG(INFO) << "[testBoundary] dxv+dyv(v1) error : " << __errDivx << "\n";
    std::cout << "[testBoundary] dxv+dyv(v1) error : " << __errDivx << "\n";
#endif

    //-----------------------------------------------------------------------------------//
    double  __errCurl = normL2( _range=markedfaces( mesh1, "Interface" ),
                                _expr=curlv( v1 )-curlv( v2, useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] curlv(v1) error : " << __errCurl );
    BOOST_CHECK_SMALL( __errCurl,1e-12 );
#else
    LOG(INFO) << "[testBoundary] curlv(u1) error : " << __errCurl << "\n";
    std::cout << "[testBoundary] curlv(u1) error : " << __errCurl << "\n";
#endif
    //-----------------------------------------------------------------------------------//
    double  __errCurlx = normL2( _range=markedfaces( mesh1, "Interface" ),
                                 _expr=curlxv( v1 )-curlxv( v2, useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] curlxv(v1) error : " << __errCurlx );
    BOOST_CHECK_SMALL( __errCurlx,1e-12 );
#else
    LOG(INFO) << "[testBoundary] curlxv(u1) error : " << __errCurlx << "\n";
    std::cout << "[testBoundary] curlxv(u1) error : " << __errCurlx << "\n";
#endif
    //-----------------------------------------------------------------------------------//
    double  __errCurly = normL2( _range=markedfaces( mesh1, "Interface" ),
                                 _expr=curlyv( v1 )-curlyv( v2, useConformalIntegration ) );
#if USE_BOOST_TEST
    BOOST_MESSAGE( "[testBoundary] curlyv(v1) error : " << __errCurly );
    BOOST_CHECK_SMALL( __errCurly,1e-12 );
#else
    LOG(INFO) << "[testBoundary] curlyv(u1) error : " << __errCurly << "\n";
    std::cout << "[testBoundary] curlyv(u1) error : " << __errCurly << "\n";
#endif

}

/*_________________________________________________*
 * run_test_interp
 *_________________________________________________*/

template<uint32_type Dim,uint32_type OrderChamp,uint32_type OrderGeo>
void
run_test_interp()
{
    typedef Mesh<Simplex<Dim,OrderGeo,Dim> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    //-----------------------------------------------------------------------------------//

    double meshSize1 = doption(_name="hsize1");
    double meshSize2 = doption(_name="hsize2");
    bool userelation = boption(_name="userelation");
    //-----------------------------------------------------------------------------------//

    //Geometry
    GeoTool::Node x1( -1,-1 );
    GeoTool::Node x2( 0,1 );
    GeoTool::Rectangle R1( meshSize1,"Omega1",x1,x2 );
    R1.setMarker( _type="line",_name="Interface",_marker2=true );
    R1.setMarker( _type="line",_name="Boundary",_marker1=true,_marker3=true,_marker4=true );
    R1.setMarker( _type="surface",_name="Omega1",_markerAll=true );

    GeoTool::Node x3( 0,-1 );
    GeoTool::Node x4( 1,1 );
    GeoTool::Rectangle R2( meshSize2,"Omega2",x3,x4 );
    R2.setMarker( _type="line",_name="Interface",_marker4=true );
    R2.setMarker( _type="line",_name="Boundary",_marker1=true,_marker2=true,_marker3=true );
    R2.setMarker( _type="surface",_name="Omega2",_markerAll=true );

    //-----------------------------------------------------------------------------------//

    GeoTool::Node x21( 0,0 );
    GeoTool::Node x22( 1,1 );
    GeoTool::Rectangle R( meshSize1,"Omega3",x21,x22 );
    R.setMarker( _type="line",_name="Boundary",_markerAll=true );
    R.setMarker( _type="surface",_name="Omega3",_markerAll=true );

    GeoTool::Node x_centre( 0.5,0.5 );
    GeoTool::Node x_bord( 0.7,0.5 );
    GeoTool::Circle C( meshSize2,"Omega4",x_centre,x_bord );
    C.setMarker( _type="line",_name="Interface",_markerAll=true );
    C.setMarker(_type="surface",_name="Omega4",_markerAll=true);
    /*
    GeoTool::Node x2_centre( 0.5,0.5 );
    GeoTool::Node x2_bord( 0.7,0.5 );
    GeoTool::Circle C2( meshSize2,"DISQUE",x2_centre,x2_bord );
    C2.setMarker( _type="line",_name="Interface",_marker1=true );
    C2.setMarker( _type="surface",_name="Omega4",_marker1=true );*/

    //-----------------------------------------------------------------------------------//

    //Mesh
    auto mesh1 = R1.createMesh(_mesh=new mesh_type,
                               _name= "mesh1_interp" + mesh_type::shape_type::name() );
    auto mesh2 = R2.createMesh(_mesh=new mesh_type,
                               _name="mesh2_interp" + mesh_type::shape_type::name() );
#if 0
    auto mesh3 = /*( R-C )*/C.createMesh(_mesh=new mesh_type,
                                    _name="mesh3_interp" + mesh_type::shape_type::name() );
#if 1
    auto mesh4 = /*( R-C )*/C.createMesh(_mesh=new mesh_type,
                                         _name="mesh4_interp" + mesh_type::shape_type::name() );
#else

    size_type ctx = userelation? EXTRACTION_KEEP_MESH_RELATION:0;
    auto mesh4 = createSubmesh(mesh3,elements(mesh3), (size_type)ctx);
#endif
#else
    auto meshBase = ((R-C)+C).createMesh(_mesh=new mesh_type,
                                         _name="meshBase_interp" + mesh_type::shape_type::name() );
    size_type ctx = userelation? EXTRACTION_KEEP_MESH_RELATION:0;
    auto mesh3 = createSubmesh(meshBase,markedelements(meshBase,"Omega3"), (size_type)ctx);
    auto mesh4 = createSubmesh(meshBase,markedelements(meshBase,"Omega4"), (size_type)ctx);
#endif
    std::list<boost::tuple<mesh_ptrtype,mesh_ptrtype,bool> > __listMesh;
    __listMesh.push_back( boost::make_tuple( mesh1,mesh2, meshSize1==meshSize2 ) );
    __listMesh.push_back( boost::make_tuple( mesh3,mesh4, true ) );

    //-----------------------------------------------------------------------------------//

    auto itMesh = __listMesh.begin();
    auto const itMesh_end = __listMesh.end();
    for ( int i=1 ; itMesh!=itMesh_end ; ++itMesh, ++i )
    {
        LOG(INFO) << "[run] start test interpolation mesh pair" << i << "\n";
        test_interp_boundary<Dim,OrderChamp,OrderGeo>( *itMesh/*, boost::mpl::int_<ScalarTest>()*/ );
        LOG(INFO) << "[run] finish test interpolation mesh pair" << i << "\n";
    }

}

/*_________________________________________________*
 * runTestExport
 *_________________________________________________*/



}//end namespace test_interp_twomesh



/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * MAIN
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_interp_twomesh::makeAbout(), test_interp_twomesh::makeOptions() )

BOOST_AUTO_TEST_SUITE( interp_twomesh_testsuite )

#if 0
BOOST_AUTO_TEST_CASE( interp_twomesh_geomap )
{
    using namespace test_interp_twomesh;

    BOOST_MESSAGE(   "\n[main] ----------TEST_GEOMAP_START----------\n" );
    BOOST_MESSAGE(   "[main] ----------------<2,6,1>---------------\n" );
    run_test_geomap<2,6,1>();
    BOOST_MESSAGE(   "[main] ----------------<2,6,2>---------------\n" );
    run_test_geomap<2,6,2>();
    BOOST_MESSAGE(   "[main] ----------------<2,6,3>---------------\n" );
    run_test_geomap<2,6,3>();
    BOOST_MESSAGE(   "[main] ----------------<2,6,4>---------------\n" );
    run_test_geomap<2,6,4>();
    //BOOST_MESSAGE(   "[main] ----------------<2,6,5>---------------\n" );
    //run_test_geomap<2,6,5>();
    BOOST_MESSAGE(   "[main] ----------TEST_GEOMAP_FINISH----------\n\n" );
}
#endif

// problem with this test case
#if 1
BOOST_AUTO_TEST_CASE( interp_twomesh_interp )
{

    using namespace test_interp_twomesh;

    BOOST_MESSAGE(   "[main] ----------TEST_INTERP_START-----------\n" );
    BOOST_MESSAGE(   "[main] ----------------<2,4,1>---------------\n" );
    run_test_interp<2,4,1>();
    BOOST_MESSAGE(   "[main] ----------------<2,4,2>---------------\n" );
    run_test_interp<2,4,2>();
    //BOOST_MESSAGE(   "[main] ----------------<2,9,3>---------------\n" );
    //run_test_interp<2,9,3>();

    //BOOST_MESSAGE(   "[main] ----------------<2,10,4>---------------\n" );
    //run_test_interp<2,10,4>();
    //BOOST_MESSAGE(   "[main] ----------------<2,11,5>---------------\n");
    //run_test_interp<2,11,5>();
    BOOST_MESSAGE(   "[main] ----------TEST_INTERP_FINISH----------\n" );
}
#endif
BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=test_interp_twomesh::makeOptions(),
                     _about=test_interp_twomesh::makeAbout() );

    Environment::changeRepository( _directory=boost::format( "/testsuite/feeldiscr/interp/%1%/" )
                                   % env.about().appName()
                                   );

    std::cout << "[main] ----------TEST_INTERP_START-----------\n";
    std::cout << "[main] ----------------<2,4,1>---------------\n";
    test_interp_twomesh::run_test_interp<2,4,1>();
    std::cout << "[main] ----------------<2,4,2>---------------\n";
    test_interp_twomesh::run_test_interp<2,4,2>();
    std::cout << "[main] ----------TEST_INTERP_FINISH----------\n";

    return 0;
}
#endif
