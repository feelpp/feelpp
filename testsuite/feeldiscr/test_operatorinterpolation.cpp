
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
    GeoTool::Node x1(0,0);
    GeoTool::Node x2(2,1);
    GeoTool::Rectangle C( meshSize,"OMEGA",x1,x2);
    C.setMarker(_type="line",_name="Sortie",_markerAll=true);
    C.setMarker(_type="surface",_name="OmegaFluide",_markerAll=true);
    auto mesh2d = C.createMesh(_mesh=new mesh_2d_type,
                               _name="test2dTo1d_domain"+mesh_2d_type::shape_type::name());

    GeoTool::Node x3(0,0);
    GeoTool::Node x4(2,1);
    GeoTool::Line L( meshSize, "Line",x3,x4);
    L.setMarker(_type="point",_name="Sortie",_markerAll=true);
    L.setMarker(_type="line",_name="Omega1d",_markerAll=true);
    auto mesh1d = L.createMesh(_mesh=new mesh_1d_type,
                               _name="test2dTo1d_domain1d"+mesh_1d_type::shape_type::name());

    //-----------------------------------------------------------//

    auto Xh1d = space_1d_type::New(_mesh=mesh1d);
    auto Xh2d = space_2d_type::New(_mesh=mesh2d);

    auto u1d = Xh1d->element();
    auto u2d = vf::project(_space=Xh2d,
                           _range=elements(mesh2d),
                           _expr=vec( cos(M_PI*Px()),sin(M_PI*Py()) ) );

    //-----------------------------------------------------------//

    auto opI=opInterpolation( _domainSpace=Xh2d,
                              _imageSpace=Xh1d);
    opI->apply(u2d,u1d);

    auto s = integrate(_range=elements(mesh1d),
                         _expr=trans(idv(u2d)-idv(u1d))*(idv(u2d)-idv(u1d)) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s,1e-8);
} // test2dTo1d

//---------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------//

template <uint32_type OrderGeo>
void
test2dTo2d( Application_ptrtype test_app)
{
    typedef Backend<double> backend_type;

    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;

    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_1_type;
    typedef FunctionSpace<mesh_type, basis_1_type> space_1_type;
    typedef typename space_1_type::element_type element_1_type;

    typedef bases<Lagrange<4,Vectorial,Continuous,PointSetFekete> > basis_2_type;
    typedef FunctionSpace<mesh_type, basis_2_type> space_2_type;
    typedef typename space_2_type::element_type element_2_type;

    //-------------------------------------------------------
    //case 1 : same mesh
    //-------------------------------------------------------
    WorldComm myWorldComm;
    auto meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1(0,0);
    GeoTool::Node x2(0.6,0);
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2);
    C.setMarker(_type="line",_name="Sortie",_markerAll=true);
    C.setMarker(_type="surface",_name="OmegaFluide",_markerAll=true);
    auto mesh = C.createMesh(_mesh = new mesh_type,
                             _name="test2dTo2d_domain"+mesh_type::shape_type::name(),
                             _partitions=myWorldComm.localSize() );

    auto Xh1 = space_1_type::New( _mesh=mesh );
    auto Xh2 = space_2_type::New( _mesh=mesh );
    auto u1 = Xh1->element("u1");
    auto u2 = Xh2->element("u2");
    auto u2a = Xh2->element("u2a");

    u1 = vf::project(_space=Xh1,
                     _range=elements(mesh),
                     _expr=vec( cos(M_PI*Px()),sin(M_PI*Py()) ) );

    auto mybackend = backend_type::build( test_app->vm() );

    auto opI=opInterpolation( _domainSpace=Xh1,
                              _imageSpace=Xh2,
                              _backend=mybackend );
    opI->apply(u1,u2);

    auto s1 = integrate(_range=elements(mesh),
                        _expr=trans(idv(u1))*idv(u1) ).evaluate()(0,0);
    auto s2 = integrate(_range=elements(mesh),
                        _expr=trans(idv(u2))*idv(u2) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s1-s2,1e-8);

    auto opIa=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2,
                               _range=boundaryfaces(mesh),
                               _backend=mybackend );
    opIa->apply(u1,u2a);

    auto s2a = integrate(_range=boundaryfaces(mesh),
                         _expr=trans(idv(u1)-idv(u2a))*(idv(u1)-idv(u2a))).evaluate()(0,0);
    BOOST_CHECK_SMALL( s2a,1e-8);

    //-------------------------------------------------------//
    //case 2 : with interpolation tool
    //-------------------------------------------------------//

    GeoTool::Circle C2( meshSize/2.,"OMEGA",x1,x2);
    C2.setMarker(_type="line",_name="Boundary",_markerAll=true);
    C2.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh2 = C2.createMesh(_mesh=new mesh_type,
                               _name="test2dTo2d_domain2"+mesh_type::shape_type::name());
    auto Xh2bis = space_2_type::New(_mesh=mesh2);
    auto u2bis = Xh2bis->element("u2bis");
    auto u2bisbis = Xh2bis->element("u2bisbis");

    auto opI2=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=elements(mesh2),
                               _backend=mybackend );
    opI2->apply(u1,u2bis);

    auto s3 = integrate(_range=elements(mesh2),
                        _expr=trans(idv(u1))*idv(u1) ).evaluate()(0,0);
    auto s4 = integrate(_range=elements(mesh2),
                        _expr=trans(idv(u2bis))*idv(u2bis) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s3-s4,1e-6);

    auto opI3=opInterpolation( _domainSpace=Xh1,
                               _imageSpace=Xh2bis,
                               _range=boundaryfaces(mesh2),
                               _backend=mybackend );
    opI3->apply(u1,u2bisbis);

    auto s7 = integrate(_range=boundaryfaces(mesh2),
                        _expr=trans(idv(u1)-idv(u2bisbis))*(idv(u1)-idv(u2bisbis)) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( s7,1e-6);
} // test2dTo2d



} // namespace test_operatorinterpolation


BOOST_AUTO_TEST_SUITE( interp_operatorinterpolation )
Environment env( boost::unit_test::framework::master_test_suite().argc,
                 boost::unit_test::framework::master_test_suite().argv );

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation )
{
    using namespace Feel::vf;
    using namespace test_operatorinterpolation;

    Application_ptrtype test_app( new Application_type(boost::unit_test::framework::master_test_suite().argc,
                                                       boost::unit_test::framework::master_test_suite().argv,
                                                       test_operatorinterpolation::makeAbout(),
                                                       test_operatorinterpolation::makeOptions()
                                                       ) );

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/" )
                                % test_app->about().appName() );

    test_operatorinterpolation::test2dTo2d<1>(test_app);
    test_operatorinterpolation::test2dTo2d<2>(test_app);
    //test_operatorinterpolation::test2dTo2d<3>(test_app);
    test_operatorinterpolation::test2dTo1d<1>(test_app);
    test_operatorinterpolation::test2dTo1d<2>(test_app);
}

BOOST_AUTO_TEST_SUITE_END()

