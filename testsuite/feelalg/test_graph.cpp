/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#define BOOST_TEST_MODULE test_graphcsr
#include <testsuite/testsuite.hpp>


#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>





namespace test_graphcsr
{

using namespace Feel;

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_graphcsr options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ;
    return desc_options.add( Feel::feel_options() ).add( backend_options( "graph" ) );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "Test_Graphcsr" ,
                     "Test_Graphcsr" ,
                     "0.1",
                     "test matrix graph and print the spy in python script",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}

/*_________________________________________________*
 * Run
 *_________________________________________________*/

void run( Application_ptrtype & theApp )
{
    //using namespace Feel;
    theApp->changeRepository( boost::format( "/testsuite/%1%/" )
                              % theApp->about().appName() );


    //--------------------------------------------------------------------------------------------------//
    // mesh
    const int nDim = 2;
    const int nOrderGeo = 1;
    double meshSize = theApp->vm()["hsize"].as<double>();

    typedef Mesh< Simplex<nDim, 1,nDim> > mesh_type;

    GeoTool::Node x1( 0,-0.5 );
    GeoTool::Node x2( 6,0.5 );
    GeoTool::Rectangle Omega( meshSize,"Omega",x1,x2 );
    Omega.setMarker( _type="line",_name="Wall",_markerAll=true );
    Omega.setMarker( _type="surface",_name="Omega",_markerAll=true );

    auto mesh = Omega.createMesh( _mesh=new mesh_type,
                                  _name="omega_"+ mesh_type::shape_type::name() );

    GeoTool::Node xc( 1.,0. );
    GeoTool::Node xr( 1.1,0. );
    GeoTool::Circle Circ( meshSize,"OmegaCirc",xc,xr );
    Circ.setMarker( _type="line",_name="Paroi",_markerAll=true );
    Circ.setMarker( _type="surface",_name="LocInterior",_markerAll=true );

    GeoTool::Node xr2( 1.3,0. );
    GeoTool::Circle Circ2( meshSize,"OmegaCirc2",xc,xr2 );
    Circ2.setMarker( _type="line",_name="FictitiousBoundary",_markerAll=true );
    Circ2.setMarker( _type="surface",_name="LocExterior",_markerAll=true );

    auto mesh2 = ( Circ2-Circ ).createMesh( _mesh=new mesh_type,
                                            _name="filename" );

    //--------------------------------------------------------------------------------------------------//
    // space
    const int nOrderPoly = 2;
    typedef Lagrange<nOrderPoly, Scalar,Continuous, PointSetFekete> basis_u_type;
    typedef Lagrange<nOrderPoly-1, Scalar, Continuous, PointSetFekete> basis_p_type;
    typedef FunctionSpace<mesh_type, bases<basis_u_type, basis_u_type, basis_p_type> > space_type;

    auto Xh1 = space_type::New( mesh );
    //std::cout << "nbDof Xh1 " << Xh1->nDof() << std::endl;

    typedef Lagrange<0, Scalar, Continuous> basis_l_type;
    typedef FunctionSpace<mesh_type, bases<basis_u_type, basis_u_type, basis_p_type, basis_l_type> > space2_type;
    auto Xh2 = space2_type::New( mesh2 );
    //std::cout << "nbDof Xh2 " << Xh2->nDof() << std::endl;

    //std::cout << "nbDof Xh1+Xh2 " << Xh1->nDof()+Xh2->nDof() << std::endl;

    //--------------------------------------------------------------------------------------------------//
    // matrix
    typedef Backend<double> backend_type;
    auto backend = backend_type::build( soption( _name="backend" ),"graph" );
    auto blockPattern = vf::Blocks<3,3,size_type>() << size_type( Pattern::COUPLED ) << size_type( Pattern::ZERO ) << size_type( Pattern::COUPLED )
                        << size_type( Pattern::ZERO ) << size_type( Pattern::COUPLED ) << size_type( Pattern::COUPLED )
                        << size_type( Pattern::COUPLED ) << size_type( Pattern::COUPLED ) << size_type( Pattern::ZERO );
    auto A_1_1 = backend->newMatrix( _test=Xh1, _trial=Xh1,
                                     _pattern_block=blockPattern,
                                     _diag_is_nonzero=false );
    A_1_1->graph()->printPython( "A_1_1.py" );

    auto A_1_2 = backend->newZeroMatrix( _test=Xh1, _trial=Xh2 );
    auto A_2_1 = backend->newZeroMatrix( _test=Xh2, _trial=Xh1 );
    auto blockPattern2 = vf::Blocks<4,4,size_type>() << size_type( Pattern::COUPLED ) << size_type( Pattern::ZERO )    << size_type( Pattern::COUPLED ) << size_type( Pattern::ZERO )
                         << size_type( Pattern::ZERO )    << size_type( Pattern::COUPLED ) << size_type( Pattern::COUPLED ) << size_type( Pattern::ZERO )
                         << size_type( Pattern::COUPLED ) << size_type( Pattern::COUPLED ) << size_type( Pattern::ZERO )    << size_type( Pattern::COUPLED )
                         << size_type( Pattern::ZERO )    << size_type( Pattern::ZERO )    << size_type( Pattern::COUPLED ) << size_type( Pattern::ZERO );
    auto A_2_2 = backend->newMatrix( _test=Xh2, _trial=Xh2,
                                     _pattern_block=blockPattern2,
                                     _diag_is_nonzero=false );
    A_2_2->graph()->printPython( "A_2_2.py" );

    auto BlockA = BlocksSparseMatrix<2,2>() << A_1_1 << A_1_2
                  << A_2_1 << A_2_2;
    auto A = backend->newBlockMatrix( _block=BlockA,
                                      _copy_values=false,
                                      _diag_is_nonzero=true );
    A->graph()->printPython( "A.py" );

} // run

} //namespace test_graphcsr


/*_________________________________________________*
 * Main
 *_________________________________________________*/

#if 1
FEELPP_ENVIRONMENT_WITH_OPTIONS( test_graphcsr::makeAbout(), test_graphcsr::makeOptions() )

BOOST_AUTO_TEST_SUITE( graphcsr )

typedef Feel::Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

BOOST_AUTO_TEST_CASE( graphcsr_case1 )
{
    using namespace Feel;

    auto theApp = Application_ptrtype( new Application_type );

    if ( theApp->vm().count( "help" ) )
    {
        std::cout << theApp->optionsDescription() << "\n";
        exit( 0 );
    }

    test_graphcsr::run( theApp );

}
BOOST_AUTO_TEST_SUITE_END()


#else

int
main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=boost::unit_test::framework::master_test_suite().argc,
                           _argv=boost::unit_test::framework::master_test_suite().argv,
                           _desc=makeOptions(),_about=makeAbout());

    typedef Feel::Application Application_type;
    typedef boost::shared_ptr<Application_type> Application_ptrtype;


    auto theApp = Application_ptrtype( new Application_type );


    if ( theApp->vm().count( "help" ) )
    {
        std::cout << theApp->optionsDescription() << "\n";
        exit( 0 );
    }

    test_graphcsr::run( theApp );

}
#endif
