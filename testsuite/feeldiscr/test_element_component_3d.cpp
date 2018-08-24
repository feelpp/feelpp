#define BOOST_TEST_MODULE test_element_component_3d
#include <testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( element_component_3d )

BOOST_AUTO_TEST_CASE( element_component_vectorial )
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);

    auto Xh = Pchv<2>( mesh );
    auto u = Xh->element( vec( cst( 1. ),cst( 2. ),cst( 3. ) ) );
    auto ux = u[Component::X];
    auto uy = u[Component::Y];
    auto uz = u[Component::Z];

    double sxRef = integrate( _range=elements( mesh ), _expr=cst(1.) ).evaluate()( 0,0 );
    double sx1 = integrate( _range=elements( mesh ), _expr=inner( idv( u ), oneX() ) ).evaluate()( 0,0 );
    double sx2 = integrate( _range=elements( mesh ), _expr=idv( ux ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sx1-sxRef,1e-12 );
    BOOST_CHECK_SMALL( sx2-sxRef,1e-12 );
    double syRef = integrate( _range=elements( mesh ), _expr=cst(2.) ).evaluate()( 0,0 );
    double sy1 = integrate( _range=elements( mesh ), _expr=inner( idv( u ),oneY() ) ).evaluate()( 0,0 );
    double sy2 = integrate( _range=elements( mesh ), _expr=idv( uy ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sy1-syRef,1e-12 );
    BOOST_CHECK_SMALL( sy2-syRef,1e-12 );
    double szRef = integrate( _range=elements( mesh ), _expr=cst(3.) ).evaluate()( 0,0 );
    double sz1 = integrate( _range=elements( mesh ), _expr=inner( idv( u ), oneZ() ) ).evaluate()( 0,0 );
    double sz2 = integrate( _range=elements( mesh ), _expr=idv( uz ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sz1-szRef,1e-12 );
    BOOST_CHECK_SMALL( sz2-szRef,1e-12 );

    auto sfull = integrate( _range=elements( mesh ), _expr=idv(u) ).evaluate();
    BOOST_CHECK_SMALL( sfull(0,0)-sxRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(1,0)-syRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(2,0)-szRef, 1e-12 );
}

template <typename SpaceType>
void
test_tensor2( std::shared_ptr<SpaceType> const& Xh )
{
    auto mesh = Xh->mesh();
    auto u = Xh->element( mat<3,3>( cst( 1. ),cst( 2. ),cst( 3. ),
                                    cst( 4. ),cst( 5. ),cst( 6. ),
                                    cst( 7. ),cst( 8. ),cst( 9. ) ) );

    auto uxx = u.comp( Component::X,Component::X );
    auto uxy = u.comp( Component::X,Component::Y );
    auto uxz = u.comp( Component::X,Component::Z );
    auto uyx = u.comp( Component::Y,Component::X );
    auto uyy = u.comp( Component::Y,Component::Y );
    auto uyz = u.comp( Component::Y,Component::Z );
    auto uzx = u.comp( Component::Z,Component::X );
    auto uzy = u.comp( Component::Z,Component::Y );
    auto uzz = u.comp( Component::Z,Component::Z );

    double sxxRef = integrate( _range=elements( mesh ), _expr=cst(1.) ).evaluate()( 0,0 );
    double sxx1 = integrate( _range=elements( mesh ), _expr=idv( u )(0,0) ).evaluate()( 0,0 );
    double sxx2 = integrate( _range=elements( mesh ), _expr=idv( uxx ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sxx1-sxxRef,1e-12 );
    BOOST_CHECK_SMALL( sxx2-sxxRef,1e-12 );
    double sxyRef = integrate( _range=elements( mesh ), _expr=cst(2.) ).evaluate()( 0,0 );
    double sxy1 = integrate( _range=elements( mesh ), _expr=idv( u )(0,1) ).evaluate()( 0,0 );
    double sxy2 = integrate( _range=elements( mesh ), _expr=idv( uxy ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sxy1-sxyRef,1e-12 );
    BOOST_CHECK_SMALL( sxy2-sxyRef,1e-12 );
    double sxzRef = integrate( _range=elements( mesh ), _expr=cst(3.) ).evaluate()( 0,0 );
    double sxz1 = integrate( _range=elements( mesh ), _expr=idv( u )(0,2) ).evaluate()( 0,0 );
    double sxz2 = integrate( _range=elements( mesh ), _expr=idv( uxz ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( sxz1-sxzRef,1e-12 );
    BOOST_CHECK_SMALL( sxz2-sxzRef,1e-12 );

    double syxRef = integrate( _range=elements( mesh ), _expr=cst(4.) ).evaluate()( 0,0 );
    double syx1 = integrate( _range=elements( mesh ), _expr=idv( u )(1,0) ).evaluate()( 0,0 );
    double syx2 = integrate( _range=elements( mesh ), _expr=idv( uyx ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( syx1-syxRef,1e-12 );
    BOOST_CHECK_SMALL( syx2-syxRef,1e-12 );
    double syyRef = integrate( _range=elements( mesh ), _expr=cst(5.) ).evaluate()( 0,0 );
    double syy1 = integrate( _range=elements( mesh ), _expr=idv( u )(1,1) ).evaluate()( 0,0 );
    double syy2 = integrate( _range=elements( mesh ), _expr=idv( uyy ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( syy1-syyRef,1e-12 );
    BOOST_CHECK_SMALL( syy2-syyRef,1e-12 );
    double syzRef = integrate( _range=elements( mesh ), _expr=cst(6.) ).evaluate()( 0,0 );
    double syz1 = integrate( _range=elements( mesh ), _expr=idv( u )(1,2) ).evaluate()( 0,0 );
    double syz2 = integrate( _range=elements( mesh ), _expr=idv( uyz ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( syz1-syzRef,1e-12 );
    BOOST_CHECK_SMALL( syz2-syzRef,1e-12 );

    double szxRef = integrate( _range=elements( mesh ), _expr=cst(7.) ).evaluate()( 0,0 );
    double szx1 = integrate( _range=elements( mesh ), _expr=idv( u )(2,0) ).evaluate()( 0,0 );
    double szx2 = integrate( _range=elements( mesh ), _expr=idv( uzx ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( szx1-szxRef,1e-12 );
    BOOST_CHECK_SMALL( szx2-szxRef,1e-12 );
    double szyRef = integrate( _range=elements( mesh ), _expr=cst(8.) ).evaluate()( 0,0 );
    double szy1 = integrate( _range=elements( mesh ), _expr=idv( u )(2,1) ).evaluate()( 0,0 );
    double szy2 = integrate( _range=elements( mesh ), _expr=idv( uzy ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( szy1-szyRef,1e-12 );
    BOOST_CHECK_SMALL( szy2-szyRef,1e-12 );
    double szzRef = integrate( _range=elements( mesh ), _expr=cst(9.) ).evaluate()( 0,0 );
    double szz1 = integrate( _range=elements( mesh ), _expr=idv( u )(2,2) ).evaluate()( 0,0 );
    double szz2 = integrate( _range=elements( mesh ), _expr=idv( uzz ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( szz1-szzRef,1e-12 );
    BOOST_CHECK_SMALL( szz2-szzRef,1e-12 );

    auto sfull = integrate( _range=elements( mesh ), _expr=idv(u) ).evaluate();
    BOOST_CHECK_SMALL( sfull(0,0)-sxxRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(0,1)-sxyRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(0,2)-sxzRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(1,0)-syxRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(1,1)-syyRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(1,2)-syzRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(2,0)-szxRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(2,1)-szyRef, 1e-12 );
    BOOST_CHECK_SMALL( sfull(2,2)-szzRef, 1e-12 );

}
template <typename SpaceType>
void
test_tensor2symm( std::shared_ptr<SpaceType> const& Xh )
{
    auto mesh = Xh->mesh();
    auto u = Xh->element();

    auto uxx = u.comp( Component::X,Component::X );
    auto uxy = u.comp( Component::X,Component::Y );
    auto uxz = u.comp( Component::X,Component::Z );
    auto uyx = u.comp( Component::Y,Component::X );
    auto uyy = u.comp( Component::Y,Component::Y );
    auto uyz = u.comp( Component::Y,Component::Z );
    auto uzx = u.comp( Component::Z,Component::X );
    auto uzy = u.comp( Component::Z,Component::Y );
    auto uzz = u.comp( Component::Z,Component::Z );

    uxx.on(_range=elements(mesh),_expr=cst(1.));
    uxy.on(_range=elements(mesh),_expr=cst(2.));
    uxz.on(_range=elements(mesh),_expr=cst(3.));
    uyy.on(_range=elements(mesh),_expr=cst(4.));
    uyz.on(_range=elements(mesh),_expr=cst(5.));
    uzz.on(_range=elements(mesh),_expr=cst(6.));

    BOOST_CHECK_CLOSE( uxx.max(), 1., 1e-12 );
    BOOST_CHECK_CLOSE( uxy.max(), 2., 1e-12 );
    BOOST_CHECK_CLOSE( uyx.max(), 2., 1e-12 );
    BOOST_CHECK_CLOSE( uxz.max(), 3., 1e-12 );
    BOOST_CHECK_CLOSE( uzx.max(), 3., 1e-12 );
    BOOST_CHECK_CLOSE( uyy.max(), 4., 1e-12 );
    BOOST_CHECK_CLOSE( uyz.max(), 5., 1e-12 );
    BOOST_CHECK_CLOSE( uzy.max(), 5., 1e-12 );
    BOOST_CHECK_CLOSE( uzz.max(), 6., 1e-12 );

    BOOST_CHECK_CLOSE( uxx.min(), 1., 1e-12 );
    BOOST_CHECK_CLOSE( uxy.min(), 2., 1e-12 );
    BOOST_CHECK_CLOSE( uyx.min(), 2., 1e-12 );
    BOOST_CHECK_CLOSE( uxz.min(), 3., 1e-12 );
    BOOST_CHECK_CLOSE( uzx.min(), 3., 1e-12 );
    BOOST_CHECK_CLOSE( uyy.min(), 4., 1e-12 );
    BOOST_CHECK_CLOSE( uyz.min(), 5., 1e-12 );
    BOOST_CHECK_CLOSE( uzy.min(), 5., 1e-12 );
    BOOST_CHECK_CLOSE( uzz.min(), 6., 1e-12 );

    uyx.on(_range=elements(mesh),_expr=cst(7.));
    uzx.on(_range=elements(mesh),_expr=cst(8.));
    uzy.on(_range=elements(mesh),_expr=cst(9.));
    BOOST_CHECK_CLOSE( uxy.min(), 7., 1e-12 );
    BOOST_CHECK_CLOSE( uxz.min(), 8., 1e-12 );
    BOOST_CHECK_CLOSE( uyz.min(), 9., 1e-12 );
    BOOST_CHECK_CLOSE( uxy.max(), 7., 1e-12 );
    BOOST_CHECK_CLOSE( uxz.max(), 8., 1e-12 );
    BOOST_CHECK_CLOSE( uyz.max(), 9., 1e-12 );
}

BOOST_AUTO_TEST_CASE( element_component_tensor2 )
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
    test_tensor2( Pchm<2>( mesh ) );
    test_tensor2( Pchm<2>( mesh, true ) );
    test_tensor2( Pdhm<2>( mesh ) );
    test_tensor2( Pdhm<2>( mesh, true ) );
}

BOOST_AUTO_TEST_CASE( element_component_tensor2symm )
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
    test_tensor2symm( Pchms<2>( mesh ) );
    test_tensor2symm( Pchms<2>( mesh, true ) );
    test_tensor2symm( Pdhms<2>( mesh ) );
    test_tensor2symm( Pdhms<2>( mesh, true ) );
}


BOOST_AUTO_TEST_SUITE_END()
