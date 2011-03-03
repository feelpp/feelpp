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





template <uint OrderGeo>
void
test2d( Application_ptrtype test_app)
{

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    //typedef bases<Lagrange<1,Scalar,Continuous,PointSetFekete> > basis_1_type;
    typedef bases<Lagrange<2,Vectorial,Continuous,PointSetFekete> > basis_1_type;
    typedef FunctionSpace<mesh_type, basis_1_type> space_1_type;
    typedef boost::shared_ptr<space_1_type> space_1_ptrtype;
    typedef typename space_1_type::element_type element_1_type;

    //typedef bases<Lagrange<1,Scalar,Continuous,PointSetFekete> > basis_2_type;
    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_2_type;
    typedef FunctionSpace<mesh_type, basis_2_type> space_2_type;
    typedef boost::shared_ptr<space_2_type> space_2_ptrtype;
    typedef typename space_2_type::element_type element_2_type;

    //-----------------------------------------------------------//

    double meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1(0,0);
    GeoTool::Node x2(0.5,0);
    //GeoTool::Node x3(0,1);
    //GeoTool::Triangle R( meshSize,"OMEGA",x1,x2,x3);
    GeoTool::Circle R( meshSize,"OMEGA",x1,x2);

    auto mesh = R.createMesh<mesh_type>("domain");


    space_1_ptrtype Xh1 = space_1_type::New( mesh );
    space_2_ptrtype Xh2 = space_2_type::New( mesh );

    element_1_type u1( Xh1, "u1" );
    element_2_type u2( Xh2, "u2" );
    element_2_type u2a( Xh2, "u2" );

    u1 = vf::project(Xh1,elements(mesh),vec( cos(M_PI*Px()),sin(M_PI*Py()) ) );
    //u1 = vf::project(Xh1,elements(mesh),cos(M_PI*Px() ) );

    auto M_backend = backend_type::build( test_app->vm() );
    //OperatorInterpolation< space_1_type,space_2_type> opI( Xh1,Xh2, M_backend );
    //opI.apply(u1,u2);
    auto opI=opInterpolation( _domainSpace=Xh1,
                              _imageSpace=Xh2,
                              //_range=elements(mesh),
                              _backend=M_backend );
    opI->apply(u1,u2);

    double s1 = integrate(elements(mesh), trans(idv(u1))*idv(u1)).evaluate()(0,0);
    double s2 = integrate(elements(mesh), trans(idv(u2))*idv(u2)).evaluate()(0,0);

    //std::cout << "\ns1 = " << s1 ;
    //std::cout << "\ns2 = " << s2 <<"\n";

    BOOST_CHECK_SMALL( s1-s2,1e-8);

    auto opIa=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2,
                               _range=boundaryfaces(mesh),
                               _backend=M_backend );
    opIa->apply(u1,u2a);

    //double s1a = integrate(boundaryfaces(mesh), trans(idv(u1))*idv(u1)).evaluate()(0,0);
    //double s2a = integrate(boundaryfaces(mesh), trans(idv(u2a))*idv(u2a)).evaluate()(0,0);
    double s2a = integrate(boundaryfaces(mesh), trans(idv(u1)-idv(u2a))*(idv(u1)-idv(u2a))).evaluate()(0,0);
    //std::cout << "\ns1 = " << s1 ;
    //std::cout << "\ns2 = " << s2 <<"\n";

    BOOST_CHECK_SMALL( s2a,1e-8);

#if 0
    auto UNexporter = Exporter<mesh_type>::New( "gmsh"/*vm*/, "ExportOOOO" );
    UNexporter->step( 0 )->setMesh( mesh );
    UNexporter->step( 0 )->add( "u1", u1 );
    UNexporter->step( 0 )->add( "u2", u2 );
    UNexporter->step( 0 )->add( "u2a", u2a );
    UNexporter->save();
#endif

    //case 2 : with interpolation tool

    GeoTool::Circle C2( meshSize/2.,"OMEGA",x1,x2);
    auto mesh2 = C2.createMesh<mesh_type>("domain2");

    space_2_ptrtype Xh2bis = space_2_type::New( mesh2 );
    element_2_type u2bis( Xh2bis, "u2bis" );
    element_2_type u2bisbis( Xh2bis, "u2bisbis" );

    //OperatorInterpolation< space_1_type,space_2_type> opI2( Xh1,Xh2bis, M_backend );
    //opI2.apply(u1,u2bis);
    auto opI2=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=elements(mesh2),
                               _backend=M_backend );
    opI2->apply(u1,u2bis);

    double s3 = integrate(elements(mesh2), trans(idv(u1))*idv(u1)).evaluate()(0,0);
    double s4 = integrate(elements(mesh2), trans(idv(u2bis))*idv(u2bis)).evaluate()(0,0);


    BOOST_CHECK_SMALL( s3-s4,1e-6);

    auto opI3=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=boundaryfaces(mesh2),
                               _backend=M_backend );
    opI3->apply(u1,u2bisbis);

    //double s5 = integrate(boundaryfaces(mesh), trans(idv(u1))*idv(u1)).evaluate()(0,0);
    //double s6 = integrate(boundaryfaces(mesh), trans(idv(u2bisbis))*idv(u2bisbis)).evaluate()(0,0);
    double s7 = integrate(boundaryfaces(mesh2), trans(idv(u1)-idv(u2bisbis))*(idv(u1)-idv(u2bisbis)) ).evaluate()(0,0);

    BOOST_CHECK_SMALL( s7,1e-6);

#if 0
    auto UNexporter2 = Exporter<mesh_type>::New( "gmsh"/*vm*/, "Export2" );
    UNexporter2->step( 0 )->setMesh( mesh2 );
    UNexporter2->step( 0 )->add( "u2bis", u2bis );
    UNexporter2->step( 0 )->add( "u2bisbis", u2bisbis );
    UNexporter2->save();
#endif
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

    test_operatorinterpolation::test2d<1>(test_app);
    test_operatorinterpolation::test2d<2>(test_app);
    test_operatorinterpolation::test2d<3>(test_app);
}

BOOST_AUTO_TEST_SUITE_END()



