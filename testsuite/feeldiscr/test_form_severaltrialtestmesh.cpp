/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#define USE_BOOST_TEST 1


#define BOOST_TEST_MODULE test_form_severaltrialtestmesh
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

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
 * Create mesh
 *_________________________________________________*/
boost::shared_ptr<Mesh< Simplex<2,1,2> > >
createTheMesh( Application_ptrtype & theApp, mpl::int_<2> /**/ )
{
    typedef Mesh< Simplex<2,1,2> > mesh_type;

    double meshSize = theApp->vm()["hsize"].as<double>();

#if 0
    // mesh
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle Omega( meshSize,"Omega",x1,x2 );
    Omega.setMarker( _type="line",_name="Paroi",_markerAll=true );
    Omega.setMarker( _type="surface",_name="Omega",_markerAll=true );
#else
    // mesh
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.5,0.5 );
    GeoTool::Circle Omega( meshSize,"Omega",x1,x2 );
    Omega.setMarker( _type="line",_name="Paroi",_markerAll=true );
    Omega.setMarker( _type="surface",_name="Omega",_markerAll=true );
#endif

    auto mesh = Omega.createMesh(_mesh=new mesh_type,
                                 _name ="omega_"+ mesh_type::shape_type::name() );
    return mesh;
}
boost::shared_ptr<Mesh< Simplex<3,1,3> > >
createTheMesh( Application_ptrtype & theApp, mpl::int_<3> /**/ )
{
    typedef Mesh< Simplex<3,1,3> > mesh_type;

    double meshSize = theApp->vm()["hsize"].as<double>();

    GeoTool::Sphere Omega(meshSize,"C1",
                          GeoTool::Node(0.,0,0),
                          GeoTool::Node(0.5,0,0) );
    Omega.setMarker(_type="surface",_name="Paroi",_markerAll=true);
    Omega.setMarker(_type="volume",_name="Sphere",_markerAll=true);

    auto mesh = Omega.createMesh(_mesh=new mesh_type,
                                 _name ="omega_"+ mesh_type::shape_type::name() );
    return mesh;
}
/*_________________________________________________*
 * Run
 *_________________________________________________*/

template <uint16_type Dim, uint16_type OrderPoly,typename Expr1,typename Expr2,typename Expr3>
void runGen( Application_ptrtype & theApp,const Expr1 & expr_f, const Expr2 & expr_g,const Expr3 & expr_gradg )
{
    const int nDim = Dim;//2;
    const int nOrderPoly = OrderPoly;
    const int nOrderGeo = 1;

    typedef Mesh< Simplex<nDim, 1,nDim> > mesh_type;
    typedef Mesh< Simplex<nDim-1, 1, nDim> > mesh_trace_type;
    //typedef Mesh< Simplex<1, 2, nDim> > mesh_1d_type2;

    double meshSize = theApp->vm()["hsize"].as<double>();
    double meshSize1d = theApp->vm()["1d-hsize"].as<double>();

    auto mesh = createTheMesh(theApp,mpl::int_<nDim>());


#if 1
    auto meshParoi = mesh->trace();//createSubmesh( mesh,markedfaces( mesh,"Paroi" ) );
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

    auto meshParoi = ( L1+L2+L3+L4 ).createMesh(_mesh=new mesh_trace_type,_name= "trace_"+ mesh_type::shape_type::name() );
#endif

    //--------------------------------------------------------------------------------------------------//

    typedef Lagrange<nOrderPoly, Scalar,Continuous,PointSetFekete> basis_u_type;
    typedef Lagrange<nOrderPoly, Scalar,Continuous,PointSetFekete> basis_l_type;

    typedef FunctionSpace<mesh_type, bases<basis_u_type> > space_u_type;
    typedef FunctionSpace<mesh_trace_type, bases<basis_l_type> > space_l_type;


    auto Xh_u = space_u_type::New( mesh );
    auto Xh_l = space_l_type::New( meshParoi );
    LOG(INFO) << "Xh_u " << Xh_u->nDof() << std::endl;
    LOG(INFO) << "Xh_l " << Xh_l->nDof() << std::endl;
    auto U_u = Xh_u->element();
    auto U_v = Xh_u->element();

    auto U_lambda = Xh_l->element();
    auto U_nu = Xh_l->element();

    //--------------------------------------------------------------------------------------------------//

    typedef Backend<double> backend_type;
    auto backend = backend_type::build( theApp->vm(),"laplacian" );

    boost::timer mytimer;
    auto A_uu = backend->newMatrix( _trial=Xh_u, _test=Xh_u, _diag_is_nonzero=false );
    LOG(INFO) << "time to build A_uu " << mytimer.elapsed() << std::endl;mytimer.restart();
    auto A_ul = backend->newMatrix( _trial=Xh_l, _test=Xh_u,  _diag_is_nonzero=false, _buildGraphWithTranspose=true,_collect_garbage=false );
    LOG(INFO) << "time to build A_ul " << mytimer.elapsed() << std::endl;mytimer.restart();
    auto A_lu = backend->newMatrix( _trial=Xh_u, _test=Xh_l, _diag_is_nonzero=false,_collect_garbage=false );
    LOG(INFO) << "time to build A_lu " << mytimer.elapsed() << std::endl;mytimer.restart();
    auto A_ll = backend->newZeroMatrix( _trial=Xh_l, _test=Xh_l);//, _diag_is_nonzero=false );
    LOG(INFO) << "time to build A_ll " << mytimer.elapsed() << std::endl;mytimer.restart();

    auto F_u = backend->newVector( Xh_u );
    auto F_l = backend->newVector( Xh_l );

    //--------------------------------------------------------------------------------------------------//

    auto solanal = vf::project( _space=Xh_u,
                                _range=elements( mesh ),
                                _expr=expr_g );


    // volume force
    //auto f = cst(1.);

    form2( Xh_u, Xh_u, A_uu ) =
        integrate( elements( mesh ), gradt( U_v )*trans( grad( U_v ) ) );
    //dsds form2( Xh_u, Xh_u, A_uu ) +=
    //    integrate( markedfaces( mesh,"Paroi" ), - gradt( U_v )*N()*id( U_v ) );
    A_uu->close();

    form2( Xh_u, Xh_l, A_ul ) =
        integrate( _range=elements(meshParoi),
                   _expr=id(U_u)*idt(U_lambda) );
        //integrate( markedfaces( mesh,"Paroi" ), id( U_u )*idt( U_lambda ) );
    A_ul->close();

    form2( Xh_l, Xh_u, A_lu ) =
        integrate( elements( meshParoi ), idt( U_u )*id( U_nu ) );
    //integrate( markedfaces(mesh,"Paroi"), idt(U_u)*id(U_nu) );
    A_lu->close();


    form1( Xh_l,  F_l ) =
        integrate( elements( meshParoi ), expr_g*id( U_nu ) ,_Q<15>() );
    //integrate( markedfaces(mesh,"Paroi"), g*id(U_nu) );
    F_l->close();

    //form2( Xh_l, Xh_l, A_ll, _init=true );
    //A_ll->close();

    //--------------------------------------------//

    form1( _test=Xh_u, _vector=F_u ) =
        integrate( _range=elements( mesh ), _expr=trans( expr_f )*id( U_v ),
                   _quad=_Q<15>() );
    F_u->close();


    //--------------------------------------------------------------------------------------------------//

    LOG(INFO) << "\n start create MatBlock" << std::endl;
    boost::timer time;

    auto myb = BlocksSparseMatrix<2,2>()<< A_uu << A_ul
                                        << A_lu << A_ll;
    auto AbB = backend->newBlockMatrix( myb );
    AbB->close();

    //form2( Xh_u, Xh_u, AbB ) +=
    //    on( markedfaces(mesh, "Paroi") , U_u, F_u, one()-one() );

    LOG(INFO) << "\n time elapsed for create MatBlock " << time.elapsed() << std::endl;

    //--------------------------------------------------------------------------------------------------//

    auto FbB = backend->newVector( F_u->size()+F_l->size(),F_u->size()+F_l->size() );
    auto UbB = backend->newVector( F_u->size()+F_l->size(),F_u->size()+F_l->size() );

    for ( size_type i = 0 ; i < F_u->size(); ++ i )
        FbB->set( i, ( *F_u )( i ) );

    for ( size_type i = 0 ; i < F_l->size(); ++ i )
        FbB->set( F_u->size()+i, ( *F_l )( i ) );

    //--------------------------------------------------------------------------------------------------//

    LOG(INFO) << "\n solve system start "<< std::endl;mytimer.restart();
    auto myprec = preconditioner( _matrix=AbB,
                                  _pc=PreconditionerType::LU_PRECOND,
                                  _backend=backend,
                                  _pcfactormatsolverpackage=MatSolverPackageType::MATSOLVER_UMFPACK );

    backend->solve( _matrix=AbB,
                    _solution=UbB,
                    _rhs=FbB,
                    _prec=myprec );
                    //_pcfactormatsolverpackage="umfpack" );
    LOG(INFO) << "\n solve system finish in"<< mytimer.elapsed() << "s"  << std::endl;

    //--------------------------------------------------------------------------------------------------//

    for ( size_type i = 0 ; i < U_u.size(); ++ i )
        U_u.set( i, ( *UbB )( i ) );

    for ( size_type i = 0 ; i < U_lambda.size(); ++ i )
        U_lambda.set( i, ( *UbB )( U_u.size()+i ) );

    //--------------------------------------------------------------------------------------------------//

    double L2error2 =integrate( _range=elements( mesh ),
                                _expr=( idv( U_u )-expr_g )*( idv( U_u )-expr_g ),
                                _quad=_Q<15>() ).evaluate()( 0,0 );
    double SemiH1error2 =integrate( _range=elements( mesh ),
                                    _expr=( gradv( U_u )-expr_gradg )*trans( gradv( U_u )-expr_gradg ),
                                    _quad=_Q<15>() ).evaluate()( 0,0 );
    double L2error =   math::sqrt( L2error2 );
    double H1error =   math::sqrt( L2error2+SemiH1error2 );

#if USE_BOOST_TEST
    BOOST_CHECK_SMALL( L2error,1e-5 );
    BOOST_CHECK_SMALL( H1error,1e-4 );
#else
    std::ostringstream ostr;
    ostr<<"testConv_P" << OrderPoly <<".data";
    std::ofstream out;
    out.open( ostr.str().c_str(), std::ios::out|std::ios::app);
    out << meshSize << " "
        << std::sqrt(L2error2) << " "
        << std::sqrt(SemiH1error2) << " "
        << std::sqrt(L2error2+SemiH1error2) << std::endl;
    out.close();
#endif
    //--------------------------------------------------------------------------------------------------//

#if 0 // export
#if 0 // oplagP1
    typedef FunctionSpace<mesh_type, bases<Lagrange<1, Scalar,Continuous,PointSetFekete> > > space_visu_u_type;
    typedef FunctionSpace<mesh_trace_type, bases<Lagrange<1, Scalar,Continuous,PointSetFekete> > > space_visu_l_type;

    OperatorLagrangeP1<space_u_type> opLagP1( Xh_u, backend );
    auto XhVisu_u = space_visu_u_type::New(opLagP1.mesh());
    auto u_visu = XhVisu_u->element();
    auto opI_u = opInterpolation(_domainSpace=Xh_u,
                                 _imageSpace=XhVisu_u,
                                 _range=elements(XhVisu_u->mesh()),
                                 _backend=backend,
                                 _type=InterpolationNonConforme() );
    opI_u->apply(U_u,u_visu);
    auto HOexporter = Exporter<mesh_type>::New( theApp->vm(), "ExportHO" );
    HOexporter->step( 0 )->setMesh( u_visu.mesh() );
    HOexporter->step( 0 )->add( "uHO", u_visu );
    HOexporter->save();

    OperatorLagrangeP1<space_l_type> opLagP1_l( Xh_l, backend );
    auto XhVisu_l = space_visu_l_type::New(opLagP1_l.mesh());
    auto l_visu = XhVisu_l->element();
    auto opI_l = opInterpolation(_domainSpace=Xh_l,
                                 _imageSpace=XhVisu_l,
                                 _range=elements(XhVisu_l->mesh()),
                                 _backend=backend,
                                 _type=InterpolationNonConforme() );
    opI_l->apply(U_lambda,l_visu);
    auto HOexporter_l = Exporter<mesh_trace_type>::New( theApp->vm(), "ExportHOTrace" );
    HOexporter_l->step( 0 )->setMesh( l_visu.mesh() );
    HOexporter_l->step( 0 )->add( "lHO", l_visu );
    HOexporter_l->save();
#else
    auto exporter = Exporter<mesh_type>::New( theApp->vm(), "Export" );
    exporter->step( 0 )->setMesh( mesh );
    exporter->step( 0 )->add( "u", U_u );
    exporter->step( 0 )->add( "solanal", solanal );
    exporter->save();
    auto exporterTrace = Exporter<mesh_trace_type>::New( theApp->vm(), "ExportTrace" );
    exporterTrace->step( 0 )->setMesh( meshParoi );
    exporterTrace->step( 0 )->add( "lambda", U_lambda );
    exporterTrace->save();
#endif
#endif

} // run


template <uint16_type OrderPoly>
void run( Application_ptrtype & theApp,mpl::int_<2> /**/ )
{

    auto pi = M_PI;
    double kk=10;
    auto g = cos(kk*pi*Px()/2.)*cos(kk*pi*Py()/2.);
    auto mygradg = trans(vec( -kk*pi/2.*sin(kk*pi*Px()/2.)*cos(kk*pi*Py()/2.),
                              -kk*pi/2.*cos(kk*pi*Px()/2.)*sin(kk*pi*Py()/2.) ));
    auto f = 2*pow(kk*pi/2.,2)*cos(kk*pi*Px()/2.)*cos(kk*pi*Py()/2.);

    runGen<2,OrderPoly>(theApp,f,g,mygradg);
}

template <uint16_type OrderPoly>
void run( Application_ptrtype & theApp,mpl::int_<3> /**/ )
{

    auto pi = M_PI;
    double kk=10;
    auto g = cos(kk*pi*Px()/2.)*cos(kk*pi*Py()/2.)*cos(kk*pi*Pz()/2.);
    auto mygradg = trans(vec( -kk*pi/2.*sin(kk*pi*Px()/2.)*cos(kk*pi*Py()/2.)*cos(kk*pi*Pz()/2.),
                              -kk*pi/2.*cos(kk*pi*Px()/2.)*sin(kk*pi*Py()/2.)*cos(kk*pi*Pz()/2.),
                              -kk*pi/2.*cos(kk*pi*Px()/2.)*cos(kk*pi*Py()/2.)*sin(kk*pi*Pz()/2.) ));
    auto f = 3*pow(kk*pi/2.,2)*cos(kk*pi*Px()/2.)*cos(kk*pi*Py()/2.)*cos(kk*pi*Pz()/2.);

    runGen<3,OrderPoly>(theApp,f,g,mygradg);
}

template <uint16_type Dim, uint16_type OrderPoly>
void run( Application_ptrtype & theApp )
{
    //using namespace Feel;
    /* change parameters below */
    const int nDim = Dim;
    const int nOrderPoly = OrderPoly;
    const int nOrderGeo = 1;

    theApp->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/%2%d/P%3%/" )
                              % theApp->about().appName()
                              % nDim
                              % OrderPoly );

    run<OrderPoly>(theApp,mpl::int_<Dim>());
}

} //namespace test_form_severaltrialtestmesh


/*_________________________________________________*
 * Main
 *_________________________________________________*/

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_form_severaltrialtestmesh::makeAbout(),
                                 test_form_severaltrialtestmesh::makeOptions() )

BOOST_AUTO_TEST_SUITE( form_severaltrialtestmesh )

typedef Feel::Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;


BOOST_AUTO_TEST_CASE( form_severaltrialtestmesh1 )
{
    auto theApp = Application_ptrtype( new Application_type );

    if ( theApp->vm().count( "help" ) )
    {
        std::cout << theApp->optionsDescription() << "\n";
        exit( 0 );
    }

    test_form_severaltrialtestmesh::run<2,4>( theApp );

}
BOOST_AUTO_TEST_SUITE_END()


#else
int
main( int argc, char** argv )
{
    typedef Feel::Application Application_type;
    typedef boost::shared_ptr<Application_type> Application_ptrtype;
    Feel::Environment env( argc,argv,
                           test_form_severaltrialtestmesh::makeAbout(),
                           test_form_severaltrialtestmesh::makeOptions() );
    auto theApp = Application_ptrtype( new Application_type );



    if ( theApp->vm().count( "help" ) )
    {
        std::cout << theApp->optionsDescription() << "\n";
        exit( 0 );
    }

    test_form_severaltrialtestmesh::run<2,4>( theApp );
    test_form_severaltrialtestmesh::run<3,2>( theApp );

}
#endif
