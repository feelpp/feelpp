#define BOOST_TEST_MODULE test_vect_comp
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <iostream>

using namespace Feel;

namespace test_vect_comp
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
    po::options_description desc_options( "test_vect_comp options" );
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
    AboutData about( "Test_Vect_Comp" ,
                     "Test_Vect_Comp" ,
                     "0.1",
                     "test composant of a field vectorial",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}


}

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_vect_comp::makeAbout(), test_vect_comp::makeOptions() )

BOOST_AUTO_TEST_SUITE( interp_vect_comp )

BOOST_AUTO_TEST_CASE( interp_vect_comp )
{

    using namespace test_vect_comp;


    auto test_app = Application_ptrtype( new Application_type );

    test_app->changeRepository( boost::format( "/testsuite/feeldiscr/%1%/" )
                                % test_app->about().appName()
                              );

    typedef Mesh<Simplex<2,1,2> > mesh_type;
    typedef boost::shared_ptr<  mesh_type > mesh_ptrtype;

    typedef bases<Lagrange<2,Vectorial,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;

    //-----------------------------------------------------------//

    double meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );

    auto mesh = R.createMesh(_mesh=new mesh_type,_name= "domain" );


    space_ptrtype Xh = space_type::New( mesh );

    element_type u( Xh, "u" );

    u = vf::project( Xh,elements( mesh ),vf::vec( vf::cst( 1. ),vf::cst( -1. ) ) );

    auto ux=u.comp<X>();
    auto uy=u.comp<Y>();

    double s1 = integrate( elements( mesh ), trans( idv( u ) )*oneX() ).evaluate()( 0,0 );
    double sx = integrate( elements( mesh ), idv( ux ) ).evaluate()( 0,0 );
    double s2 = integrate( elements( mesh ), trans( idv( u ) )*oneY() ).evaluate()( 0,0 );
    double sy = integrate( elements( mesh ), idv( uy ) ).evaluate()( 0,0 );

    BOOST_CHECK_SMALL( s1-sx,1e-8 );
    BOOST_CHECK_SMALL( s2-sy,1e-8 );
}

BOOST_AUTO_TEST_SUITE_END()



