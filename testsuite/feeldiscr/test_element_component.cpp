#define BOOST_TEST_MODULE test_element_component
#include <testsuite/testsuite.hpp>

#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( element_component )

BOOST_AUTO_TEST_CASE( element_component_vectorial )
{
    typedef Mesh<Simplex<2,1,2> > mesh_type;
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle R( doption(_name="gmsh.hsize"),"OMEGA",x1,x2 );
    auto mesh = R.createMesh(_mesh=new mesh_type,_name= "domain" );

    auto Xh = Pchv<2>( mesh );
    auto u = Xh->element( vec( cst( 1. ),cst( 2. ) ) );
    auto ux = u[Component::X];
    auto uy = u[Component::Y];

    double sxRef = integrate( _range=elements( mesh ), _expr=cst(1.) ).evaluate()( 0,0 );
    double sx1 = integrate( _range=elements( mesh ), _expr=trans( idv( u ) )*oneX() ).evaluate()( 0,0 );
    double sx2 = integrate( _range=elements( mesh ), _expr=idv( ux ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sx1-sxRef,1e-12 );
    BOOST_CHECK_SMALL( sx2-sxRef,1e-12 );
    double syRef = integrate( _range=elements( mesh ), _expr=cst(2.) ).evaluate()( 0,0 );
    double sy1 = integrate( _range=elements( mesh ), _expr=trans( idv( u ) )*oneY() ).evaluate()( 0,0 );
    double sy2 = integrate( _range=elements( mesh ), _expr=idv( uy ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sy1-syRef,1e-12 );
    BOOST_CHECK_SMALL( sy2-syRef,1e-12 );

    auto sfull = integrate( _range=elements( mesh ), _expr=idv(u) ).evaluate();
    BOOST_CHECK_SMALL( sfull(0,0)-sxRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(1,0)-syRef, 1e-12 );
}


template<typename SpaceT>
void
test_tensor2()
{
    typedef Mesh<Simplex<2,1,2> > mesh_type;
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle R( doption(_name="gmsh.hsize"),"OMEGA",x1,x2 );
    auto mesh = R.createMesh(_mesh=new mesh_type,_name= "domain" );

    auto VhTensor2 = SpaceT::New( mesh );
    auto uTensor2 = VhTensor2->element();
    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( cst(1.),cst(2.),cst(3.),cst(4.) ) );
    auto uxx = uTensor2.comp( Component::X,Component::X );
    auto uxy = uTensor2.comp( Component::X,Component::Y );
    auto uyx = uTensor2.comp( Component::Y,Component::X );
    auto uyy = uTensor2.comp( Component::Y,Component::Y );

    double sxxRef = integrate( _range=elements( mesh ), _expr=cst(1.) ).evaluate()( 0,0 );
    double sxx1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneX(),oneX() ) ).evaluate()( 0,0 );
    double sxx2 = integrate( _range=elements( mesh ), _expr=idv(uxx) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sxx1-sxxRef, 1e-12 );
    BOOST_CHECK_SMALL( sxx2-sxxRef, 1e-12 );
    double sxyRef = integrate( _range=elements( mesh ), _expr=cst(2.) ).evaluate()( 0,0 );
    double sxy1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneY(),oneX() ) ).evaluate()( 0,0 );
    double sxy2 = integrate( _range=elements( mesh ), _expr=idv(uxy) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sxy1-sxyRef, 1e-12 );
    BOOST_CHECK_SMALL( sxy2-sxyRef, 1e-12 );
    double syxRef = 0;
    if ( SpaceT::is_tensor2symm )
        syxRef = integrate( _range=elements( mesh ), _expr=cst(2.) ).evaluate()( 0,0 );
    else
        syxRef = integrate( _range=elements( mesh ), _expr=cst(3.) ).evaluate()( 0,0 );

    double syx1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneX(),oneY() ) ).evaluate()( 0,0 );
    double syx2 = integrate( _range=elements( mesh ), _expr=idv(uyx) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( syx1-syxRef, 1e-12 );
    BOOST_CHECK_SMALL( syx2-syxRef, 1e-12 );
    double syyRef = integrate( _range=elements( mesh ), _expr=cst(4.) ).evaluate()( 0,0 );
    double syy1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneY(),oneY() ) ).evaluate()( 0,0 );
    double syy2 = integrate( _range=elements( mesh ), _expr=idv(uyy) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( syy1-syyRef, 1e-12 );
    BOOST_CHECK_SMALL( syy2-syyRef, 1e-12 );

    auto sfull = integrate( _range=elements( mesh ), _expr=idv(uTensor2) ).evaluate();
    BOOST_CHECK_SMALL( sfull(0,0)-sxxRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(0,1)-sxyRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(1,0)-syxRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(1,1)-syyRef, 1e-12 );

}
BOOST_AUTO_TEST_CASE( element_component_tensor2_continuous )
{
    test_tensor2<Pchm_type<Mesh<Simplex<2>>,2>>();
}
BOOST_AUTO_TEST_CASE( element_component_tensor2_discontinuous )
{
    test_tensor2<Pdhm_type<Mesh<Simplex<2>>,2>>();
}

BOOST_AUTO_TEST_CASE( element_component_tensor2symm_continuous )
{
    test_tensor2<Pchms_type<Mesh<Simplex<2>>,2>>();
}


BOOST_AUTO_TEST_CASE( element_component_tensor2symm_discontinuous )
{
    test_tensor2<Pdhms_type<Mesh<Simplex<2>>,2>>();
}

BOOST_AUTO_TEST_SUITE_END()
