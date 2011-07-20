#define BOOST_TEST_MODULE test_normal3d
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <feel/options.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <iostream>

using namespace Feel;
using namespace Feel::vf;



namespace test_normal3d
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
    po::options_description desc_options("test_normal3d options");
    desc_options.add_options()
        ("hsize", po::value<double>()->default_value( 0.4 ), "mesh size")
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
    AboutData about( "Test_normal3d" ,
                     "Test_normal3d" ,
                     "0.1",
                     "test normal 3d",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier");

    about.addAuthor("Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "");
    return about;
}

template <uint orderGeo = 1>
void
runtest(Application_ptrtype test_app)
{
    typedef Simplex<3,orderGeo,3> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    typedef bases<Lagrange<2+orderGeo,Vectorial,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //typedef typename space_type::element_type element_type;

    //-----------------------------------------------------------//

    double meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node Centre(0,0,0);
    GeoTool::Node Rayon( 1);
    GeoTool::Node Dir(1,0,0);
    GeoTool::Node Lg(5,0,0);
    GeoTool::Cylindre C( meshSize,"UnCylindre",Centre,Dir,Rayon,Lg);
    C.setMarker(_type="surface",_name="Inlet",_marker1=true);
    C.setMarker(_type="surface",_name="Outlet",_marker2=true);
    C.setMarker(_type="surface",_name="Wall",_marker3=true);
    C.setMarker(_type="volume",_name="OmegaFluid",_markerAll=true);

    std::ostringstream __ostrData;
    __ostrData<<convex_type::name()<<"h"<<meshSize;

    auto mesh = C.createMesh<mesh_type>("domain"+__ostrData.str());

    //-----------------------------------------------------------//

    auto Xh = space_type::New( mesh );
    auto u = Xh->element();

    u = vf::project(Xh,elements(mesh),vf::vec( vf::cos(M_PI*Px()/5.),vf::cos(M_PI*Py()/5.),vf::cos(M_PI*Py()/5.)));


    auto value1 = integrate(elements(mesh),divv(u), _quad=_Q<15>()).evaluate()(0,0);
    auto value2 = integrate(boundaryfaces(mesh),trans(idv(u))*N(),_quad=_Q<15>()).evaluate()(0,0);

    std::cout << "\n value 1 =" << value1;
    std::cout << "\n value 2 =" << value2 <<"\n";

    //-----------------------------------------------------------//

    BOOST_CHECK_SMALL( value1-value2,1e-8);


#if 0
    auto nnn = Xh->element();

    nnn = vf::project(Xh,markedfaces(mesh,"Wall"),N());


    auto UNexporter = Exporter<mesh_type>::New( "gmsh"/*test_app->vm()*/, "ExportOOOO"+__ostrData.str() );
    UNexporter->step( 0 )->setMesh( mesh );
    UNexporter->step( 0 )->add( "u", u );
    UNexporter->step( 0 )->add( "n", nnn );
    UNexporter->save();
#endif
}

}

BOOST_AUTO_TEST_SUITE( normal3d )


BOOST_AUTO_TEST_CASE( normal3d )
{

    using namespace test_normal3d;

    using namespace Feel::vf;

    auto test_app = Application_ptrtype( new Application_type(boost::unit_test::framework::master_test_suite().argc,
                                                              boost::unit_test::framework::master_test_suite().argv,
                                                              makeAbout(),
                                                              makeOptions()
                                                              ) );

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/" )
                                % test_app->about().appName()
        );


    runtest<1>(test_app);
    runtest<2>(test_app);
    runtest<3>(test_app);
    runtest<4>(test_app);


}

BOOST_AUTO_TEST_SUITE_END()



