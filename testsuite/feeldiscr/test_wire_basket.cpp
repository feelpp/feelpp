/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#define BOOST_TEST_MODULE test_wire_basket
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;


#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feeldiscr/operatorlift.hpp>


namespace test_wire_basket
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
    po::options_description desc_options("test_wire_basket options");
    desc_options.add_options()
        ("hsize", po::value<double>()->default_value( 0.075 ), "mesh size")
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
    AboutData about( "Test_Wire_Basket" ,
                     "Test_Wire_Basket" ,
                     "0.1",
                     "test wire basket for three fields domain decomposition method ",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier");

    about.addAuthor("Abdoulaye Samake", "developer", "Abdoulaye.Samake@imag.fr", "");
    return about;
}

/*_________________________________________________*
 * Run
 *_________________________________________________*/

template <uint16_type OrderPoly>
void run(Application_ptrtype & theApp)
{
    //using namespace Feel;
    /* change parameters below */
    const int nDim = 3;
    const int nOrderPoly = OrderPoly;
    double meshSize = theApp->vm()["hsize"].as<double>();

    theApp->changeRepository( boost::format( "testsuite/feeldiscr/%1%/P%2%/h_%3%/" )
                              % theApp->about().appName()
                              % OrderPoly
                              % meshSize );

    //--------------------------------------------------------------------------------------------------//

    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef Mesh< Simplex<nDim,1,nDim> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<OrderPoly,Scalar> > > space_type;
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

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name=(boost::format( "hypercube-%1%" ) % nDim).str() ,
                                              _addmidpoint=false, _usenames=false, _shape="hypercube",
                                              _dim=nDim, _h=meshSize, _xmin=0., _xmax=1., _ymin=0.,
                                              _ymax=1., _zmin=0., _zmax=1.) );


    auto backend = backend_type::build( theApp->vm() );
    auto pi = M_PI;
    auto g = sin(pi*(2*Px()+Py()+1./4))*cos(pi*(Py()-1./4));


    auto Xh = space_type::New( mesh );

    auto domain_measure = integrate(elements(mesh),cst(1.) ).evaluate()(0,0);
    std::cout <<"domain_measure= " << domain_measure << std::endl;


    auto trace_mesh = createSubmesh(mesh,markedfaces(mesh,6));

    auto TXh = trace_space_type::New( trace_mesh );

    auto trace_measure = integrate(elements(trace_mesh),cst(1.) ).evaluate()(0,0);
    std::cout <<"trace_measure= " << trace_measure << std::endl;


    auto trace_trace_mesh = createSubmesh(trace_mesh,boundaryfaces(trace_mesh));

    auto TTXh = trace_trace_space_type::New( trace_trace_mesh );

    auto trace_trace_measure = integrate(elements(trace_trace_mesh),cst(1.) ).evaluate()(0,0);
    auto trace_trace_integrate_g = integrate(elements(trace_trace_mesh),g ).evaluate()(0,0);
    std::cout <<"trace_trace_measure= " << trace_trace_measure << std::endl;

    // projections

    auto projection_g = vf::project( Xh, elements(mesh), g );
    auto trace_projection_g = vf::project( TXh, elements(trace_mesh), g );
    auto trace_trace_projection_g = vf::project( TTXh, elements(trace_trace_mesh), g );

    // extensions

    auto trace_trace_integrate = integrate(elements(trace_trace_mesh),
                                           idv(trace_trace_projection_g) ).evaluate()(0,0);
    auto ttmean_g = trace_trace_integrate/trace_trace_measure;
    auto mean_g = trace_trace_integrate_g/trace_trace_measure;
    std::cout << "mean_g= " << mean_g << std::endl;
    std::cout << "ttmean_g= " << ttmean_g << std::endl;


    auto zero_extension = vf::project( TXh, boundaryfaces(trace_mesh), idv(trace_trace_projection_g) );
    auto const_extension = vf::project( TXh, boundaryfaces(trace_mesh), idv(trace_trace_projection_g)-ttmean_g );
    const_extension += vf::project( TXh, elements(trace_mesh), cst(ttmean_g) );
    auto op_lift = operatorLift(Xh,backend);
    auto glift = op_lift->lift( _range=markedfaces(mesh,6),_expr=idv(const_extension));


    auto boundary_error = integrate( markedfaces(mesh,6), idv(const_extension)-idv(glift) ).evaluate()(0,0);
    auto laplacian_error = integrate( elements(mesh), trace(hessv(glift)) ).evaluate()(0,0);

    auto const_extention_error1 = integrate( boundaryfaces(trace_mesh),
                                             idv(const_extension)-idv(trace_trace_projection_g) ).evaluate()(0,0);
    auto const_extention_error2 = integrate( elements(trace_mesh), idv(const_extension)-ttmean_g ).evaluate()(0,0);



    std::cout << "boundary_error= " << boundary_error << std::endl;
    std::cout << "laplacian_error= " << laplacian_error << std::endl;

    std::cout << "const_extention_error1= " << const_extention_error1 << std::endl;
    std::cout << "const_extention_error2= " << const_extention_error2 << std::endl;

    BOOST_CHECK_SMALL( boundary_error,1e-5);
    BOOST_CHECK_SMALL( laplacian_error,6e-3);
    BOOST_CHECK_CLOSE( domain_measure, 1, 1e-10 );
    BOOST_CHECK_CLOSE( trace_measure, 1, 1e-12 );
    BOOST_CHECK_CLOSE( trace_trace_measure, 4, 1e-12 );
    BOOST_CHECK_SMALL( const_extention_error1,1e-10);
    BOOST_CHECK_SMALL( const_extention_error2,5e-4);

    auto exporter = export_type::New( theApp->vm(), "Export" );

    auto trace_exporter = trace_export_type::New( theApp->vm(), "Trace_Export" ) ;

    auto trace_trace_exporter = trace_trace_export_type::New( theApp->vm(), "Trace_Trace_Export");


    exporter->step(0)->setMesh( mesh );
    exporter->step(0)->add( "g", projection_g );
    exporter->step(0)->add( "glift", glift );
    exporter->save();

    trace_exporter->step(0)->setMesh( trace_mesh );
    trace_exporter->step(0)->add( "traceg", trace_projection_g );
    trace_exporter->step(0)->add( "const_extension", const_extension );
    trace_exporter->save();

    trace_trace_exporter->step(0)->setMesh( trace_trace_mesh );
    trace_trace_exporter->step(0)->add( "tracetrace_g", trace_trace_projection_g );
    trace_trace_exporter->save();

} // run

} //namespace test_wire_basket


/*_________________________________________________*
 * Main
 *_________________________________________________*/

BOOST_AUTO_TEST_SUITE( wire_basket )

typedef Feel::Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                       boost::unit_test::framework::master_test_suite().argv );

BOOST_AUTO_TEST_CASE( wire_basket1 )
{
    auto theApp = Application_ptrtype( new Application_type(boost::unit_test::framework::master_test_suite().argc,
                                                            boost::unit_test::framework::master_test_suite().argv,
                                                            test_wire_basket::makeAbout(),
                                                            test_wire_basket::makeOptions()
                                                              ) );

    if ( theApp->vm().count( "help" ) )
        {
            std::cout << theApp->optionsDescription() << "\n";
            exit(0);
        }

    test_wire_basket::run<2>(theApp);

}
BOOST_AUTO_TEST_SUITE_END()


#if 0
int
main( int argc, char** argv )
{
    auto theApp = Application_ptrtype( new Application_type(argc,argv,
                                                            test_wire_basket::makeAbout(),
                                                            test_wire_basket::makeOptions() ) );


    if ( theApp->vm().count( "help" ) )
        {
            std::cout << theApp->optionsDescription() << "\n";
            exit(0);
        }

    test_wire_basket::run<2>(theApp);

}
#endif
