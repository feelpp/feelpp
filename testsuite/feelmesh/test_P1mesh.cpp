
#define BOOST_TEST_MODULE test_P1mesh
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>



using namespace Feel;
using namespace Feel::vf;

namespace test_P1mesh
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
    po::options_description desc_options( "test_P1mesh options" );
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
    AboutData about( "Test_P1mesh" ,
                     "Test_P1mesh" ,
                     "0.1",
                     "test P1mesh",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2012 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}

template <uint32_type OrderGeo>
void
test2dP1mesh( Application_ptrtype test_app )
{
    typedef Backend<double> backend_type;
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;

    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;

    typedef Mesh<Simplex<2,1,2> > mesh_P1_type;
    typedef bases<Lagrange<1,Scalar,Continuous,PointSetFekete> > basis_P1_type;
    typedef FunctionSpace<mesh_P1_type, basis_P1_type> space_P1_type;

    WorldComm myWorldComm;
    auto meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.6,0 );
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name="test2dP1mesh_domain"+mesh_type::shape_type::name(),
                              _partitions=myWorldComm.localSize() );
    auto meshP1 = mesh->createP1mesh();
    auto XhP1 = space_P1_type::New( _mesh=meshP1 );
    auto u = XhP1->element();
    auto v = XhP1->element();

    auto mybackend = backend_type::build(test_app->vm());
    auto A = mybackend->newMatrix(_test=XhP1, _trial=XhP1);
    auto F = mybackend->newVector(XhP1);
    form2( _test=XhP1, _trial=XhP1, _matrix=A ) =
        integrate( _range=elements( meshP1 ),
                   _expr=gradt(u)*trans(grad(v)) );
    form1( _test=XhP1, _vector=F ) =
        integrate( _range=elements( meshP1 ),
                   _expr=id(v) );
    A->close();F->close();
    form2( _test=XhP1, _trial=XhP1, _matrix=A ) +=
        on ( _range=boundaryfaces(meshP1),
             _element=u,
             _rhs=F,
             _expr=cst(0.) );
    mybackend->solve(_matrix=A,_solution=u,_rhs=F);

    auto myexporter = Exporter<mesh_P1_type>::New( test_app->vm(), "test2dP1mesh_MyExport" );
    myexporter->step(0)->setMesh( meshP1 );
    myexporter->step(0)->add( "test2dP1mesh_uP1", u );
    myexporter->save();
}


template <uint32_type OrderGeo>
void
test2dP1meshComposite( Application_ptrtype test_app )
{
    typedef Backend<double> backend_type;
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;
    typedef Mesh<Simplex<2,1,2> > mesh_P1_type;

    typedef Lagrange<3,Vectorial,Continuous,PointSetFekete> basis_u_type;
    typedef Lagrange<2,Scalar,Continuous,PointSetFekete> basis_p_type;
    typedef FunctionSpace<mesh_type, bases<basis_u_type,basis_p_type> > space_type;

    typedef bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > basis_P1_type;
    typedef FunctionSpace<mesh_P1_type, basis_P1_type> space_P1_type;

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


    auto meshP1 = mesh->createP1mesh();
    auto XhP1 = space_P1_type::New( _mesh=meshP1,_worldscomm=vecLocWorldComm );
#if 0
    auto uP1 = XhP1->element();

    auto mybackend = backend_type::build(test_app->vm());

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
                              _backend=mytbackend );
    opI->apply( u,uLagP1 );

    auto s1 = integrate(_range=elements(mesh),
                        _expr=trans(idv(u)-idv(uLagP1))*(idv(u)-idv(uLagP1)) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s1,1e-6);
#endif
#if 0
    auto myexporter = Exporter<mesh_P1_type>::New( test_app->vm(), "test2dP1mesh_MyExport", myWorldComm );
    myexporter->step(0)->setMesh( meshP1 );
    myexporter->step(0)->add( "test2dP1mesh_uHO", uP1 );
    myexporter->save();
#endif


}


} // namespace test_P1mesh

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_P1mesh::makeAbout(), test_P1mesh::makeOptions() )
/**
 * main code
 */
BOOST_AUTO_TEST_SUITE( interp_P1mesh )

BOOST_AUTO_TEST_CASE( interp_P1mesh )
{
    using namespace Feel::vf;
    using namespace test_P1mesh;

    Application_ptrtype test_app( new Application_type( boost::unit_test::framework::master_test_suite().argc,
                                  boost::unit_test::framework::master_test_suite().argv,
                                  test_P1mesh::makeAbout(),
                                  test_P1mesh::makeOptions()
                                                      ) );

    test_app->changeRepository( boost::format( "/testsuite/feelmesh/%1%/" )
                                % test_app->about().appName() );

    test_P1mesh::test2dP1mesh<2>( test_app);
    //test_P1mesh::test2dP1meshComposite<2>( test_app);
}

BOOST_AUTO_TEST_SUITE_END()

