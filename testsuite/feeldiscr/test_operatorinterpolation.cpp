
#define BOOST_TEST_MODULE test_operatorinterpolation
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <feel/feeldiscr/operatorinterpolation.hpp>


using namespace Feel;
using namespace Feel::vf;



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
    po::options_description desc_options("test_operatorinterpolation options");
    desc_options.add_options()
        ("hsize", po::value<double>()->default_value( 0.1 ), "mesh size")
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
                     "Copyright (c) 2010 Universite Joseph Fourier");

    about.addAuthor("Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "");
    return about;
}


template <uint32_type OrderGeo>
void
test2dTo1d( Application_ptrtype test_app)
{


    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef Mesh<Simplex<1,OrderGeo,2> > mesh_1d_type;
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_2d_type;
    typedef boost::shared_ptr< mesh_2d_type > mesh_ptrtype;

    //typedef bases<Lagrange<2,Scalar,Continuous,PointSetFekete> > basis_2d_type;
    typedef bases<Lagrange<2,Vectorial,Continuous,PointSetFekete> > basis_2d_type;
    typedef FunctionSpace<mesh_2d_type, basis_2d_type> space_2d_type;
    typedef boost::shared_ptr<space_2d_type> space_2d_ptrtype;
    typedef typename space_2d_type::element_type element_2d_type;

    //typedef bases<Lagrange<3,Scalar,Continuous,PointSetFekete> > basis_1d_type;
    typedef bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > basis_1d_type;
    typedef FunctionSpace<mesh_1d_type, basis_1d_type> space_1d_type;
    typedef boost::shared_ptr<space_1d_type> space_1d_ptrtype;
    typedef typename space_1d_type::element_type element_1d_type;

    //-----------------------------------------------------------//

    double meshSize = test_app->vm()["hsize"].as<double>();
#if 0
    GeoTool::Node x1(0,0);
    GeoTool::Node x2(0.5,0);
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2);

    GeoTool::Node x3(-0.5,0);
    GeoTool::Node x4(0.5,0);
    GeoTool::Line L( meshSize, "Line",x3,x4);
#else
    GeoTool::Node x1(0,0);
    GeoTool::Node x2(2,1);
    GeoTool::Rectangle C( meshSize,"OMEGA",x1,x2);

    GeoTool::Node x3(0,0);
    GeoTool::Node x4(2,1);
    GeoTool::Line L( meshSize, "Line",x3,x4);
#endif
    GeoTool::Node x5(0.2,0.2);
    GeoTool::Node x6(0.8,0.8);
    GeoTool::Rectangle Rbis( meshSize,"koko",x5,x6);
    auto mesh2dBis = Rbis.createMesh<mesh_2d_type>("domainBis");
    auto Xh2dBis = space_2d_type::New(mesh2dBis);

    auto mesh2d = C.createMesh<mesh_2d_type>("domain");
    auto mesh1d = L.createMesh<mesh_1d_type>("domain1d");

    auto Xh1d = space_1d_type::New(mesh1d);
    auto Xh2d = space_2d_type::New(mesh2d);

    //auto u1d = Xh1d->element();
    //auto u1d = vf::project(Xh1d,elements(mesh1d),vec( cos(M_PI*Px()),sin(M_PI*Py()) ) );

    //auto u1d = vf::project(Xh1d,elements(mesh1d), cst(0.) );
    //auto u2d = vf::project(Xh2d,elements(mesh2d), cos(M_PI*Px()) );
    auto u1d = vf::project(Xh1d,elements(mesh1d),vec( cst(0.),cst(0.) ) );
    auto u2d = vf::project(Xh2d,elements(mesh2d),vec( cos(M_PI*Px()),sin(M_PI*Py()) ) );




    auto M_backend = backend_type::build( test_app->vm() );


    auto opI=opInterpolation( _domainSpace=Xh2d,
                              //_imageSpace=Xh2dBis,
                              _imageSpace=Xh1d,
                              //_range=elements(mesh),
                              _backend=M_backend );
    //size check
    std::cout << "\n size1 "<< opI->matPtr()->size1();
    std::cout << "\n size2 "<< opI->matPtr()->size2();
    std::cout << "\n size u2d "<< u2d.size();
    std::cout << "\n size u1d "<< u1d.size();
#if 1
    opI->apply(u2d,u1d);

    for (uint32_type i=0;i<u1d.size();++i)
        std::cout << "\n i = "<< i << " val : " << u1d(i);

#endif

    double s = integrate(elements(mesh1d), trans(idv(u2d)-idv(u1d))*(idv(u2d)-idv(u1d))).evaluate()(0,0);

    std::cout << "\ns = " << s ;
#if 1

    //BOOST_CHECK_SMALL( s1-s2,1e-8);
#if 0
    auto UNexporter2d = Exporter<mesh_2d_type>::New( test_app->vm(), "Export2d" );
    UNexporter2d->step( 0 )->setMesh( mesh2d );
    UNexporter2d->step( 0 )->add( "u2d", u2d );
    UNexporter2d->save();
#endif
    u1d = vf::project(Xh1d,elements(mesh1d),vec( cos(M_PI*Px()),sin(M_PI*Py()) ) );
    auto UNexporter1d = Exporter<mesh_1d_type>::New( test_app->vm(), "Export1d" );
    UNexporter1d->step( 0 )->setMesh( mesh1d );
    UNexporter1d->step( 0 )->add( "u1d", u1d );
    //    UNexporter1d->step( 0 )->add( "u2dbis", u2d );
    UNexporter1d->save();
#endif


}

//---------------------------------------------------------------------------------------------//



template <uint32_type OrderGeo>
void
test2d( Application_ptrtype test_app)
{

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_1_type;
    typedef FunctionSpace<mesh_type, basis_1_type> space_1_type;
    typedef boost::shared_ptr<space_1_type> space_1_ptrtype;
    typedef typename space_1_type::element_type element_1_type;

    typedef bases<Lagrange<4,Vectorial,Continuous,PointSetFekete> > basis_2_type;
    typedef FunctionSpace<mesh_type, basis_2_type> space_2_type;
    typedef boost::shared_ptr<space_2_type> space_2_ptrtype;
    typedef typename space_2_type::element_type element_2_type;

    //-------------------------------------------------------
    //case 1 : same mesh
    //-------------------------------------------------------

    auto meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1(0,0);
    GeoTool::Node x2(0.6,0);
    GeoTool::Circle R( meshSize,"OMEGA",x1,x2);

    auto mesh = R.createMesh(_mesh = new mesh_type,_name="domain");
    auto Xh1 = space_1_type::New( _mesh=mesh );
    auto Xh2 = space_2_type::New( _mesh=mesh );
    auto u1 = Xh1->element("u1");
    auto u2 = Xh2->element("u2");
    auto u2a = Xh2->element("u2a");

    u1 = vf::project(_space=Xh1,
                     _range=elements(mesh),
                     _expr=vec( cos(M_PI*Px()),sin(M_PI*Py()) ) );

    auto M_backend = backend_type::build( test_app->vm() );
    auto opI=opInterpolation( _domainSpace=Xh1,
                              _imageSpace=Xh2,
                              _backend=M_backend );
    opI->apply(u1,u2);

    auto s1 = integrate(_range=elements(mesh),
                        _expr=trans(idv(u1))*idv(u1) ).evaluate()(0,0);
    auto s2 = integrate(_range=elements(mesh),
                        _expr=trans(idv(u2))*idv(u2) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s1-s2,1e-8);

    auto opIa=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2,
                               _range=boundaryfaces(mesh),
                               _backend=M_backend );
    opIa->apply(u1,u2a);

    auto s2a = integrate(_range=boundaryfaces(mesh),
                         _expr=trans(idv(u1)-idv(u2a))*(idv(u1)-idv(u2a))).evaluate()(0,0);
    BOOST_CHECK_SMALL( s2a,1e-8);

    //-------------------------------------------------------//
    //case 2 : with interpolation tool
    //-------------------------------------------------------//

    GeoTool::Circle C2( meshSize/2.,"OMEGA",x1,x2);
    auto mesh2 = C2.createMesh(_mesh=new mesh_type,_name="domain2");
    auto Xh2bis = space_2_type::New(_mesh=mesh2);
    auto u2bis = Xh2bis->element("u2bis");
    auto u2bisbis = Xh2bis->element("u2bisbis");

    auto opI2=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=elements(mesh2),
                               _backend=M_backend );
    opI2->apply(u1,u2bis);

    auto s3 = integrate(_range=elements(mesh2),
                        _expr=trans(idv(u1))*idv(u1) ).evaluate()(0,0);
    auto s4 = integrate(_range=elements(mesh2),
                        _expr=trans(idv(u2bis))*idv(u2bis) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s3-s4,1e-6);

    auto opI3=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=boundaryfaces(mesh2),
                               _backend=M_backend );
    opI3->apply(u1,u2bisbis);

    auto s7 = integrate(_range=boundaryfaces(mesh2),
                        _expr=trans(idv(u1)-idv(u2bisbis))*(idv(u1)-idv(u2bisbis)) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s7,1e-6);
}



} // end test_operatorinterpolation

BOOST_AUTO_TEST_SUITE( interp_operatorinterpolation )

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation )
{

    using namespace test_operatorinterpolation;

    using namespace Feel::vf;

    auto test_app = Application_ptrtype( new Application_type(boost::unit_test::framework::master_test_suite().argc,
                                                              boost::unit_test::framework::master_test_suite().argv,
                                                              makeAbout(),
                                                              makeOptions()
                                                              ) );

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/" )
                                % test_app->about().appName()
        );
#if 1
    test_operatorinterpolation::test2d<1>(test_app);
    test_operatorinterpolation::test2d<2>(test_app);
    test_operatorinterpolation::test2d<3>(test_app);
#else
    test_operatorinterpolation::test2dTo1d<1>(test_app);
#endif
}

BOOST_AUTO_TEST_SUITE_END()



