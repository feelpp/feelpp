/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#define BOOST_TEST_MODULE test_wirebasket
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;


#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>

namespace test_wirebasket
{

    using namespace Feel;
    using namespace Feel::vf;

    typedef Application Application_type;
    typedef boost::shared_ptr<Application_type> Application_ptrtype;

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_wire_basket options" );
    desc_options.add_options()
        ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
        ;
    return desc_options.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "Test_WireBasket" ,
                     "Test_WireBasket" ,
                     "0.1",
                     "test wire basket for three fields domain decomposition method ",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Abdoulaye Samake", "developer", "Abdoulaye.Samake@imag.fr", "" );
    return about;
}

template <uint16_type Order>
void run( Application_ptrtype & theApp )
{

    const int nDim = 3;
    const int nOrder = Order;
    double meshSize = theApp->vm()["hsize"].as<double>();

    theApp->changeRepository( boost::format( "testsuite/feeldiscr/%1%/P%2%/h_%3%/" )
                              % theApp->about().appName()
                              % Order
                              % meshSize );

    //--------------------------------------------------------------------------------------------------//

    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef Mesh< Simplex<nDim,1,nDim> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<Order,Scalar> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    // trace
    typedef typename mesh_type::trace_mesh_type trace_mesh_type;
    typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;
    typedef typename space_type::trace_functionspace_type trace_space_type;
    typedef typename boost::shared_ptr<trace_space_type> trace_space_ptrtype;
    typedef Exporter<trace_mesh_type> trace_export_type;
    typedef boost::shared_ptr<trace_export_type> trace_export_ptrtype;

    // trace_trace
    typedef typename trace_mesh_type::trace_mesh_type trace_trace_mesh_type;
    typedef typename trace_mesh_type::trace_mesh_ptrtype trace_trace_mesh_ptrtype;
    typedef typename trace_space_type::trace_functionspace_type trace_trace_space_type;
    typedef typename boost::shared_ptr<trace_trace_space_type> trace_trace_space_ptrtype;
    typedef Exporter<trace_trace_mesh_type> trace_trace_export_type;
    typedef boost::shared_ptr<trace_trace_export_type> trace_trace_export_ptrtype;

    //--------------------------------------------------------------------------------------------------//


    auto mesh1DFrom3D = createGMSHMesh( _mesh=new trace_trace_mesh_type,
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                        _desc=geo( _filename="domain.geo",_h=meshSize ) );


    auto mesh2DFrom3D = createGMSHMesh( _mesh=new trace_mesh_type,
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                        _desc=geo( _filename="domain.geo",_h=meshSize ) );

    Log() << "number of elements 1D/3D: " << mesh1DFrom3D->numElements() << "\n";
    Log() << "number of elements 2D/3D: " << mesh2DFrom3D->numElements() << "\n";
#if 1
#if 0
    auto mesh1DFrom2D = createGMSHMesh( _mesh=new trace_trace_mesh_type,
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                        _desc=geo( _filename="domain.geo",_h=meshSize ) );

    auto wirebasket_measure1F2 = integrate( elements( mesh1DFrom2D ),cst( 1. ) ).evaluate()( 0,0 );
#endif

    auto wirebasket_measure1F3 = integrate( elements( mesh1DFrom3D ),cst( 1. ) ).evaluate()( 0,0 );
    auto wirebasket_measure2F3 = integrate( elements( mesh2DFrom3D ),cst( 1. ) ).evaluate()( 0,0 );

#if 0
    BOOST_CHECK_SMALL( boundary_error,5e-5 );
    BOOST_CHECK_CLOSE( domain_measure, 1, 1e-10 );
    BOOST_CHECK_CLOSE( trace_measure, 1, 1e-12 );
    BOOST_CHECK_CLOSE( trace_trace_measure, 4, 1e-12 );
    BOOST_CHECK_SMALL( const_extention_error1,1e-10 );
    BOOST_CHECK_SMALL( const_extention_error2,5e-4 );
#endif

#endif

}

} //namespace test_wirebasket


BOOST_AUTO_TEST_SUITE( wire_basket )

typedef Feel::Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                       boost::unit_test::framework::master_test_suite().argv );

BOOST_AUTO_TEST_CASE( wire_basket1 )
{
    auto theApp = Application_ptrtype( new Application_type( boost::unit_test::framework::master_test_suite().argc,
                                                             boost::unit_test::framework::master_test_suite().argv,
                                                             test_wirebasket::makeAbout(),
                                                             test_wirebasket::makeOptions()
                                                             ) );

    if ( theApp->vm().count( "help" ) )
        {
            std::cout << theApp->optionsDescription() << "\n";
            exit( 0 );
        }

    test_wirebasket::run<2>( theApp );

}
BOOST_AUTO_TEST_SUITE_END()
