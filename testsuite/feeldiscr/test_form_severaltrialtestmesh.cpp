/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#define BOOST_TEST_MODULE test_form_severaltrialtestmesh
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;


#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

//#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/matrixblock.hpp>

#include <feel/feeldiscr/createsubmesh.hpp>

//#include <benchmarks/fbm/integratelocal2.hpp>




namespace test_form_severaltrialtestmesh
{

using namespace Feel;
using namespace Feel::vf;

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_form_severaltrialtestmesh options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.02 ), "mesh size" )
    ( "1d-hsize", po::value<double>()->default_value( 0.02 ), "mesh size 1d" )
    ;
    return desc_options.add( Feel::feel_options() ).add( backend_options( "laplacian" ) );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "Test_Form_SeveralTrialTestMesh" ,
                     "Test_Form_SeveralTrialTestMesh" ,
                     "0.1",
                     "test bilinear and linear form for several trial/test mesh ",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}

/*_________________________________________________*
 * Run
 *_________________________________________________*/

template <uint16_type OrderPoly>
void run( Application_ptrtype & theApp )
{
    //using namespace Feel;

    theApp->changeRepository( boost::format( "/fsitest/%1%/P%2%/" )
                              % theApp->about().appName()
                              % OrderPoly );


    /* change parameters below */
    const int nDim = 2;
    const int nOrderPoly = OrderPoly;
    const int nOrderGeo = 1;

    //--------------------------------------------------------------------------------------------------//

    typedef Mesh< Simplex<nDim, 1,nDim> > mesh_type;
    typedef Mesh< Simplex<1, 1, nDim> > mesh_1d_type;
    //typedef Mesh< Simplex<1, 2, nDim> > mesh_1d_type2;

    double meshSize = theApp->vm()["hsize"].as<double>();

    double meshSize1d = theApp->vm()["1d-hsize"].as<double>();

    // mesh
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle Omega( meshSize,"Omega",x1,x2 );
    Omega.setMarker( _type="line",_name="Paroi",_markerAll=true );
    Omega.setMarker( _type="surface",_name="Omega",_markerAll=true );

    auto mesh = Omega.createMesh<mesh_type>( "omega_"+ mesh_type::shape_type::name() );

#if 1
    auto meshParoi = createSubmesh( mesh,markedfaces( mesh,"Paroi" ) );
#else
    GeoTool::Node c1( 0,0 );
    GeoTool::Node c2( 1,0 );
    GeoTool::Node c3( 1,1 );
    GeoTool::Node c4( 0,1 );
    GeoTool::Line L1( meshSize1d,"line1",c1,c2 );
    GeoTool::Line L2( meshSize1d,"line2",c2,c3 );
    GeoTool::Line L3( meshSize1d,"line3",c3,c4 );
    GeoTool::Line L4( meshSize1d,"line4",c4,c1 );
    /*L1.setMarker(_type="point",_name="corners",_markerAll=true);
    L1.setMarker(_type="line",_name="OmegaParoi",_markerAll=true);
    L2.setMarker(_type="point",_name="corners",_markerAll=true);
    L2.setMarker(_type="line",_name="OmegaParoi",_markerAll=true);
    L3.setMarker(_type="point",_name="corners",_markerAll=true);
    L3.setMarker(_type="line",_name="OmegaParoi",_markerAll=true);
    L4.setMarker(_type="point",_name="corners",_markerAll=true);
    L4.setMarker(_type="line",_name="OmegaParoi",_markerAll=true);*/

    auto meshParoi2 = ( L1+L2+L3+L4 ).createMesh<mesh_1d_type>( "trace_"+ mesh_type::shape_type::name() );
#endif

    //--------------------------------------------------------------------------------------------------//

    typedef Lagrange<nOrderPoly, Scalar,Continuous,PointSetFekete> basis_u_type;
    typedef Lagrange<nOrderPoly, Scalar,Continuous> basis_l_type;

    typedef FunctionSpace<mesh_type, bases<basis_u_type> > space_u_type;
    typedef FunctionSpace<mesh_1d_type, bases<basis_l_type> > space_l_type;

    auto Xh_u = space_u_type::New( mesh );
    auto Xh_l = space_l_type::New( meshParoi );

    auto U_u = Xh_u->element();
    auto U_v = Xh_u->element();

    auto U_lambda = Xh_l->element();
    auto U_nu = Xh_l->element();

    //--------------------------------------------------------------------------------------------------//

    typedef Backend<double> backend_type;
    auto backend = backend_type::build( theApp->vm(),"laplacian" );

    auto A_uu = backend->newMatrix( _trial=Xh_u, _test=Xh_u, _diag_is_nonzero=false );
    auto A_ul = backend->newMatrix( _trial=Xh_l, _test=Xh_u, _buildGraphWithTranspose=true );

    auto A_lu = backend->newMatrix( _trial=Xh_u, _test=Xh_l, _diag_is_nonzero=false );
    auto A_ll = backend->newMatrix( _trial=Xh_l, _test=Xh_l, _diag_is_nonzero=false );

    auto F_u = backend->newVector( Xh_u );
    auto F_l = backend->newVector( Xh_l );

    //--------------------------------------------------------------------------------------------------//

    auto pi = M_PI;
    auto g = sin( pi*Px() )*cos( pi*Py() ); //*cos(pi*Pz());
    auto f = 2*pi*pi*g;


    auto solanal = vf::project( Xh_u,elements( mesh ), g );

    // volume force
    //auto f = cst(1.);

    form2( Xh_u, Xh_u, A_uu ) =
        integrate( elements( mesh ), gradt( U_v )*trans( grad( U_v ) ) );
    form2( Xh_u, Xh_u, A_uu ) +=
        integrate( markedfaces( mesh,"Paroi" ), - gradt( U_v )*N()*id( U_v ) );
    A_uu->close();


    form2( Xh_u, Xh_l, A_ul ) =
        integrate( elements(meshParoi), id(U_u)*idt(U_lambda) );
        //integrate( markedfaces( mesh,"Paroi" ), id( U_u )*idt( U_lambda ) );
    A_ul->close();

    form2( Xh_l, Xh_u, A_lu ) =
        integrate( elements( meshParoi ), idt( U_u )*id( U_nu ) );
    //integrate( markedfaces(mesh,"Paroi"), idt(U_u)*id(U_nu) );
    A_lu->close();

    form1( Xh_l,  F_l ) =
        integrate( elements( meshParoi ), g*id( U_nu ) ,_Q<15>() );
    //integrate( markedfaces(mesh,"Paroi"), g*id(U_nu) );
    F_l->close();


    //form2( Xh_l, Xh_l, A_ll, _init=true );
    //A_ll->close();

    //--------------------------------------------//

    form1( Xh_u, F_u, _init=true ) =
        integrate( elements( mesh ), trans( f )*id( U_v ),_Q<15>() );
    F_u->close();

    //--------------------------------------------------------------------------------------------------//

    std::cout << "\n start create MatBlock" << std::endl;
    boost::timer time;

    auto myb = Blocks<2,2>()<< A_uu << A_ul
               << A_lu << A_ll;
    auto AbB = backend->newBlockMatrix( myb );
    AbB->close();

    //form2( Xh_u, Xh_u, AbB ) +=
    //    on( markedfaces(mesh, "Paroi") , U_u, F_u, one()-one() );

    std::cout << "\n time elapsed for create MatBlock " << time.elapsed() << std::endl;

    //--------------------------------------------------------------------------------------------------//

    auto FbB = backend->newVector( F_u->size()+F_l->size(),F_u->size()+F_l->size() );
    auto UbB = backend->newVector( F_u->size()+F_l->size(),F_u->size()+F_l->size() );

    for ( size_type i = 0 ; i < F_u->size(); ++ i )
        FbB->set( i, ( *F_u )( i ) );

    for ( size_type i = 0 ; i < F_l->size(); ++ i )
        FbB->set( F_u->size()+i, ( *F_l )( i ) );

    //--------------------------------------------------------------------------------------------------//

    std::cout << "\n solve system start "<< std::endl;
    backend->solve( _matrix=AbB,
                    _solution=UbB,
                    _rhs=FbB,
                    _pcfactormatsolverpackage="umfpack" );
    std::cout << "\n solve system finish "<< std::endl;

    //--------------------------------------------------------------------------------------------------//

    for ( size_type i = 0 ; i < U_u.size(); ++ i )
        U_u.set( i, ( *UbB )( i ) );

    for ( size_type i = 0 ; i < U_lambda.size(); ++ i )
        U_lambda.set( i, ( *UbB )( U_u.size()+i ) );

    //--------------------------------------------------------------------------------------------------//

    double L2error2 =integrate( elements( mesh ),
                                ( idv( U_u )-g )*( idv( U_u )-g ) ).evaluate()( 0,0 );
    double L2error =   math::sqrt( L2error2 );

    BOOST_CHECK_SMALL( L2error,1e-5 );


    auto exporter = Exporter<mesh_type>::New( theApp->vm(), "Export" );

    exporter->step( 0 )->setMesh( mesh );
    exporter->step( 0 )->add( "u", U_u );
    exporter->step( 0 )->add( "solanal", solanal );
    exporter->save();





} // run

} //namespace test_form_severaltrialtestmesh


/*_________________________________________________*
 * Main
 *_________________________________________________*/

BOOST_AUTO_TEST_SUITE( form_severaltrialtestmesh )

typedef Feel::Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                       boost::unit_test::framework::master_test_suite().argv );

BOOST_AUTO_TEST_CASE( form_severaltrialtestmesh1 )
{
    auto theApp = Application_ptrtype( new Application_type( boost::unit_test::framework::master_test_suite().argc,
                                       boost::unit_test::framework::master_test_suite().argv,
                                       test_form_severaltrialtestmesh::makeAbout(),
                                       test_form_severaltrialtestmesh::makeOptions()
                                                           ) );

    if ( theApp->vm().count( "help" ) )
    {
        std::cout << theApp->optionsDescription() << "\n";
        exit( 0 );
    }

    test_form_severaltrialtestmesh::run<2>( theApp );

}
BOOST_AUTO_TEST_SUITE_END()


#if 0
int
main( int argc, char** argv )
{
    auto theApp = Application_ptrtype( new Application_type( argc,argv,
                                       test_form_severaltrialtestmesh::makeAbout(),
                                       test_form_severaltrialtestmesh::makeOptions() ) );


    if ( theApp->vm().count( "help" ) )
    {
        std::cout << theApp->optionsDescription() << "\n";
        exit( 0 );
    }

    test_form_severaltrialtestmesh::run<2>( theApp );

}
#endif
