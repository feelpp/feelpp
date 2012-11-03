
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
test2dTo1d( Application_ptrtype test_app )
{
    typedef Backend<double> backend_type;

    typedef Mesh<Simplex<1,OrderGeo,2> > mesh_1d_type;
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_2d_type;

    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_2d_type;
    typedef FunctionSpace<mesh_2d_type, basis_2d_type> space_2d_type;
    typedef typename space_2d_type::element_type element_2d_type;

    typedef bases<Lagrange<5,Vectorial,Continuous,PointSetFekete> > basis_1d_type;
    typedef FunctionSpace<mesh_1d_type, basis_1d_type> space_1d_type;
    typedef typename space_1d_type::element_type element_1d_type;

    //-----------------------------------------------------------//

    auto meshSize = test_app->vm()["hsize"].as<double>();
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 2,1 );
    GeoTool::Rectangle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh2d = C.createMesh( _mesh=new mesh_2d_type,
                                _name="test2dTo1d_domain"+mesh_2d_type::shape_type::name() );

    GeoTool::Node x3( 0,0 );
    GeoTool::Node x4( 2,1 );
    GeoTool::Line L( meshSize, "Line",x3,x4 );
    L.setMarker( _type="point",_name="Sortie",_markerAll=true );
    L.setMarker( _type="line",_name="Omega1d",_markerAll=true );
    auto mesh1d = L.createMesh( _mesh=new mesh_1d_type,
                                _name="test2dTo1d_domain1d"+mesh_1d_type::shape_type::name() );

    //-----------------------------------------------------------//

    auto Xh1d = space_1d_type::New( _mesh=mesh1d );
    auto Xh2d = space_2d_type::New( _mesh=mesh2d );

    auto u1d = Xh1d->element();
    auto u2d = vf::project( _space=Xh2d,
                            _range=elements( mesh2d ),
                            _expr=vec( cos( M_PI*Px() ),sin( M_PI*Py() ) ) );

    //-----------------------------------------------------------//

    auto opI=opInterpolation( _domainSpace=Xh2d,
                              _imageSpace=Xh1d );
    opI->apply( u2d,u1d );

    auto s = integrate( _range=elements( mesh1d ),
                        _expr=trans( idv( u2d )-idv( u1d ) )*( idv( u2d )-idv( u1d ) ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s,1e-8 );
} // test2dTo1d

//---------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------//

template <uint32_type OrderGeo>
void
test2dTo2d( Application_ptrtype test_app )
{
    typedef Backend<double> backend_type;

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
    WorldComm myWorldComm;
    auto meshSize = test_app->vm()["hsize"].as<double>();

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

    u1 = vf::project( _space=Xh1,
                      _range=elements( mesh ),
                      _expr=vec( cos( M_PI*Px() ),sin( M_PI*Py() ) ) );
                      //_expr=cos(M_PI*Px() ) );

    auto mybackend = backend_type::build( test_app->vm() );

    auto opI=opInterpolation( _domainSpace=Xh1,
                              _imageSpace=Xh2,
                              _backend=mybackend );
    opI->apply( u1,u2 );

    auto s1 = integrate( _range=elements( mesh ),
                         _expr=trans( idv( u1 ) )*idv( u1 ) ).evaluate()( 0,0 );
    auto s2 = integrate( _range=elements( mesh ),
                         _expr=trans( idv( u2 ) )*idv( u2 ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s1-s2,1e-8 );

    auto opIa=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2,
                               _range=boundaryfaces( mesh ),
                               _backend=mybackend );
    opIa->apply( u1,u2a );

    auto s2a = integrate( _range=boundaryfaces( mesh ),
                          _expr=trans( idv( u1 )-idv( u2a ) )*( idv( u1 )-idv( u2a ) ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s2a,1e-8 );

    //-------------------------------------------------------//
    //case 2 : with interpolation tool
    //-------------------------------------------------------//
    //GeoTool::Node x3(0.4,0);
    GeoTool::Circle C2( meshSize/2.,"OMEGA",x1,x2);
    C2.setMarker(_type="line",_name="Boundary",_markerAll=true);
    C2.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh2 = C2.createMesh(_mesh=new mesh_type,
                               _name="test2dTo2d_domain2"+mesh_type::shape_type::name(),
                               _partitions=myWorldComm.localSize() );
    auto Xh2bis = space_2_type::New(_mesh=mesh2);
    auto u2bis = Xh2bis->element();
    auto u2bisbis = Xh2bis->element();

    auto opI2=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=elements( mesh2 ),
                               _backend=mybackend );
    opI2->apply( u1,u2bis );

#if 0
    auto myexporter = Exporter<mesh_type>::New( test_app->vm(), "test2dTo2d_MyExport" );
    myexporter->step(0)->setMesh( u2bis.mesh() );
    myexporter->step(0)->add( "test2dTo2d_u2bis", u2bis );
    myexporter->save();
    auto myexporter2 = Exporter<mesh_type>::New( test_app->vm(), "test2dTo2d_MyExport2" );
    myexporter2->step(0)->setMesh( u1.mesh() );
    myexporter2->step(0)->add( "test2dTo2d_u1", u1 );
    myexporter2->save();
#endif

    auto s3 = integrate(_range=elements(mesh2),
                        _expr=trans(idv(u1))*idv(u1) ).evaluate()(0,0);
    auto s4 = integrate(_range=elements(mesh2),
                        _expr=trans(idv(u2bis))*idv(u2bis) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s3-s4,1e-6);

    auto opI3=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=boundaryfaces( mesh2 ),
                               _backend=mybackend );
    opI3->apply( u1,u2bisbis );

    auto s7 = integrate( _range=boundaryfaces( mesh2 ),
                         _expr=trans( idv( u1 )-idv( u2bisbis ) )*( idv( u1 )-idv( u2bisbis ) ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( s7,1e-6 );

} // test2dTo2d


#if defined(FEELPP_HAS_VTK)
template <uint32_type OrderGeo>
void
test2dOpLagrangeP1( Application_ptrtype test_app )
{
    typedef Backend<double> backend_type;
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;
    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > basis_P1_type;
    typedef FunctionSpace<mesh_type, basis_P1_type> space_P1_type;

    WorldComm myWorldComm;
    auto meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.6,0 );
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name="test2dOpLagrangeP1_domain"+mesh_type::shape_type::name(),
                              _partitions=myWorldComm.localSize() );

    auto Xh = space_type::New( _mesh=mesh );
    auto u = vf::project( _space=Xh,
                          _range=elements( mesh ),
                          _expr=vec( cos( M_PI*Px() ),sin( M_PI*Py() ) ) );

    auto mybackend = backend_type::build(test_app->vm());

    auto opLagP1 = lagrangeP1(_space=Xh);
    auto meshLagP1 = opLagP1->mesh();

    auto XhLagP1 = space_P1_type::New( _mesh=meshLagP1 );
    auto uLagP1 = XhLagP1->element();

    auto opI=opInterpolation( _domainSpace=Xh,
                              _imageSpace=XhLagP1,
                              _range=elements( meshLagP1 ) );
    opI->apply( u,uLagP1 );

    auto s1 = integrate(_range=elements(mesh),
                        _expr=trans(idv(u)-idv(uLagP1))*(idv(u)-idv(uLagP1)) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s1,1e-6);

#if 0
    auto myexporter = Exporter<mesh_type>::New( test_app->vm(), "test2dOpLagrangeP1_MyExport" );
    myexporter->step(0)->setMesh( meshLagP1 );
    myexporter->step(0)->add( "test2dOpLagrangeP1_uLagP1", uLagP1 );
    myexporter->save();
#endif
}

template <uint32_type OrderGeo>
void
test2dOpLagrangeP1Composite( Application_ptrtype test_app )
{
    typedef Backend<double> backend_type;
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;

    typedef Lagrange<3,Vectorial,Continuous,PointSetFekete> basis_u_type;
    typedef Lagrange<2,Scalar,Continuous,PointSetFekete> basis_p_type;
    typedef FunctionSpace<mesh_type, bases<basis_u_type,basis_p_type> > space_type;

    typedef bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > basis_P1_type;
    typedef FunctionSpace<mesh_type, basis_P1_type> space_P1_type;

    //-------------------------------
    const int VelocityWorld=0;
    const int PressureWorld=1;
    std::vector<int> MapWorld(test_app->comm().size());
    WorldComm myWorldComm;
    if (test_app->comm().size()>1)
        {
            for (int proc = 0 ; proc < test_app->comm().size(); ++proc)
                {
                    if (proc < test_app->comm().size()/2 ) // if (proc%2==0 )
                        MapWorld[proc] = VelocityWorld;
                    else
                        MapWorld[proc] = PressureWorld;
                }
            myWorldComm = WorldComm(MapWorld);
        }

    //-------------------------------

    auto meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.6,0 );
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name="test2dOpLagrangeP1_domain"+mesh_type::shape_type::name(),
                              _partitions=myWorldComm.localSize(),
                              _worldcomm=myWorldComm);


    //-------------------------------

    std::vector<WorldComm> vecWorldComm(space_type::nSpaces);
    std::vector<WorldComm> vecLocWorldComm(1);

    int CurrentWorld=0;
    if (myWorldComm.globalRank() < myWorldComm.globalSize()/2 )
        CurrentWorld=VelocityWorld;
    else
        CurrentWorld=PressureWorld;

    if (myWorldComm.globalSize()>1)
        {
            vecWorldComm[0]=myWorldComm.subWorldComm(VelocityWorld);
            vecWorldComm[1]=myWorldComm.subWorldComm(PressureWorld);
            vecLocWorldComm[0]=myWorldComm.subWorldComm(CurrentWorld);
        }
    else
        {
            vecWorldComm[0]=WorldComm();
            vecWorldComm[1]=WorldComm();
            vecLocWorldComm[0]=WorldComm();
        }

    //-------------------------------
#if defined(FEELPP_ENABLE_MPI_MODE)
    auto Xh = space_type::New( _mesh=mesh, _worldscomm=vecWorldComm );
#else
    auto Xh = space_type::New( _mesh=mesh );
#endif
    auto U = Xh->element();
    auto u = U.template element<0>();
    u = vf::project( _space=Xh->template functionSpace<0>(),
                     _range=elements( mesh ),
                     _expr=vec( cos( M_PI*Px() ),sin( M_PI*Py() ) ) );

    auto mybackend = backend_type::build(test_app->vm());

    //OperatorLagrangeP1<typename space_type::template sub_functionspace<0>::type::element_type> opLagP1( Xh->template functionSpace<0>(), mybackend, vecLocWorldComm );
    auto opLagP1 = lagrangeP1(_space=Xh->template functionSpace<0>(),
                              _backend=mybackend,
                              _worldscomm=vecLocWorldComm);
    auto meshLagP1 = opLagP1->mesh();

#if defined(FEELPP_ENABLE_MPI_MODE)
    auto XhLagP1 = space_P1_type::New( _mesh=meshLagP1,_worldscomm=vecLocWorldComm );
#else
    auto XhLagP1 = space_P1_type::New( _mesh=meshLagP1 );
#endif
    auto uLagP1 = XhLagP1->element();

    auto opI=opInterpolation( _domainSpace=Xh->template functionSpace<0>(),
                              _imageSpace=XhLagP1,
                              _range=elements( meshLagP1 ),
                              _backend=mybackend );
    opI->apply( u,uLagP1 );

    auto s1 = integrate(_range=elements(mesh),
                        _expr=trans(idv(u)-idv(uLagP1))*(idv(u)-idv(uLagP1)) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s1,1e-6);

#if 0
    auto myexporter = Exporter<mesh_type>::New( test_app->vm(), "test2dOpLagrangeP1_MyExport", myWorldComm );
    myexporter->step(0)->setMesh( meshLagP1 );
    myexporter->step(0)->add( "test2dOpLagrangeP1_uHO", uLagP1 );
    myexporter->save();
#endif


}


#endif

} // namespace test_operatorinterpolation

/**
 * main code
 */
FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( interp_operatorinterpolation )

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation )
{
    using namespace Feel::vf;
    using namespace test_operatorinterpolation;

    Application_ptrtype test_app( new Application_type( boost::unit_test::framework::master_test_suite().argc,
                                  boost::unit_test::framework::master_test_suite().argv,
                                  test_operatorinterpolation::makeAbout(),
                                  test_operatorinterpolation::makeOptions()
                                                      ) );

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/" )
                                % test_app->about().appName() );

    test_operatorinterpolation::test2dTo2d<1>( test_app );
    test_operatorinterpolation::test2dTo2d<2>( test_app );
    //test_operatorinterpolation::test2dTo2d<3>(test_app);
    test_operatorinterpolation::test2dTo1d<1>( test_app );
    test_operatorinterpolation::test2dTo1d<2>( test_app );
#if defined(FEELPP_HAS_VTK)
    test_operatorinterpolation::test2dOpLagrangeP1<1>( test_app);
    test_operatorinterpolation::test2dOpLagrangeP1Composite<1>( test_app);
#endif
}

BOOST_AUTO_TEST_SUITE_END()

