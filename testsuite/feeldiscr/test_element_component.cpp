#define BOOST_TEST_MODULE test_element_component
#include <testsuite/testsuite.hpp>

#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdhv.hpp>
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
    using Feel::cout;
    typedef Mesh<Simplex<2,1,2> > mesh_type;
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Node x3( 0,1 );
    GeoTool::Rectangle R( doption(_name="gmsh.hsize"),"OMEGA",x1,x2 );
    //GeoTool::Triangle R( 3,"OMEGA",x1,x2,x3 );
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
    BOOST_CHECK_CLOSE( sxx1, sxxRef, 1e-12 );
    BOOST_CHECK_CLOSE( sxx2, sxxRef, 1e-12 );
    double sxyRef = 0;
    if ( SpaceT::is_tensor2symm )
        sxyRef = integrate( _range=elements( mesh ), _expr=cst(3.) ).evaluate()( 0,0 );
    else
        sxyRef = integrate( _range=elements( mesh ), _expr=cst(2.) ).evaluate()( 0,0 );
    double sxy1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneY(),oneX() ) ).evaluate()( 0,0 );
    double sxy2 = integrate( _range=elements( mesh ), _expr=idv(uxy) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( sxy1, sxyRef, 1e-12 );
    BOOST_CHECK_CLOSE( sxy2, sxyRef, 1e-12 );
    double syxRef = 0;
    syxRef = integrate( _range=elements( mesh ), _expr=cst(3.) ).evaluate()( 0,0 );

    double syx1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneX(),oneY() ) ).evaluate()( 0,0 );
    double syx2 = integrate( _range=elements( mesh ), _expr=idv(uyx) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( syx1, syxRef, 1e-12 );
    BOOST_CHECK_CLOSE( syx2, syxRef, 1e-12 );
    double syyRef = integrate( _range=elements( mesh ), _expr=cst(4.) ).evaluate()( 0,0 );
    double syy1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneY(),oneY() ) ).evaluate()( 0,0 );
    double syy2 = integrate( _range=elements( mesh ), _expr=idv(uyy) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( syy1, syyRef, 1e-12 );
    BOOST_CHECK_CLOSE( syy2, syyRef, 1e-12 );

    auto sfull = integrate( _range=elements( mesh ), _expr=idv(uTensor2) ).evaluate();
    BOOST_CHECK_CLOSE( sfull(0,0),sxxRef, 1e-12 );
    BOOST_CHECK_CLOSE( sfull(0,1),sxyRef, 1e-12 );
    BOOST_CHECK_CLOSE( sfull(1,0),syxRef, 1e-12 );
    BOOST_CHECK_CLOSE( sfull(1,1),syyRef, 1e-12 );

    double area = integrate( _range=elements( mesh ), _expr=cst(1.) ).evaluate()(0,0);
    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Px(),cst(0.),cst(3.),2.3*Py() ) );
    auto sdiv = integrate( _range=elements( mesh ), _expr=divv(uTensor2) ).evaluate();
    BOOST_CHECK_CLOSE( sdiv(0,0),area, 1e-12 );
    BOOST_CHECK_CLOSE( sdiv(1,0),2.3*area, 1e-12 );

    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Py(),Px(),Px(),2.3*Px() ) );
    sdiv = integrate( _range=elements( mesh ), _expr=divv(uTensor2) ).evaluate();
    BOOST_CHECK_SMALL( sdiv(0,0), 1e-12 );
    BOOST_CHECK_CLOSE( sdiv(1,0),area, 1e-12 );

    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Px(),Py(),Py(),2.3*Py()+Px() ) );
    sdiv = integrate( _range=elements( mesh ), _expr=divv(uTensor2) ).evaluate();
    BOOST_CHECK_CLOSE( sdiv(0,0),2*area, 1e-12 );
    BOOST_CHECK_CLOSE( sdiv(1,0),2.3*area, 1e-12 );


#if 0
    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Px()*Py(),Py()+Px(),Py(),2.3*Py()+Px() ) );
    sdiv = integrate( _range=elements( mesh ), _expr=divv(uTensor2) ).evaluate();
    BOOST_CHECK_CLOSE( sdiv(0,0),2*area, 1e-12 );
    BOOST_CHECK_CLOSE( sdiv(1,0),2.3*area, 1e-12 );
#endif

    /*
     * Test bilinear forms for HDG linear elasticity
     */

    auto Wh = Pdhv<3>( mesh );
    auto w = Wh->element();

    w.on(_range=elements(mesh), _expr=ones<2,1>());
    uTensor2.on(_range=elements(mesh),_expr= ones<2,2>() );

    // a11

    auto a11a = form2( _trial=VhTensor2, _test=VhTensor2 );
    a11a += integrate( _range=elements( mesh ), _expr=inner(idt(uTensor2),id(uTensor2)));
    a11a.close();
    auto a11a_v = a11a(uTensor2, uTensor2);
    BOOST_CHECK_CLOSE( a11a_v, 4*area, 1e-11 );
    cout << "a11a(ones, ones) = " << a11a_v << std::endl;

    auto a11b = form2( _trial=VhTensor2, _test=VhTensor2 );
    a11b += integrate( _range=elements( mesh ), _expr=trace(idt(uTensor2))*trace(id(uTensor2)));
    a11b.close();
    auto a11b_v = a11b(uTensor2, uTensor2);
    BOOST_CHECK_CLOSE( a11b_v, 4*area, 1e-11 );
    cout << "a11b(ones, ones) = " << a11b_v << std::endl;

    cout << "a11 works" << std::endl;

    // a12

    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Px(),Py(),Px(),Py() ) );

    auto a12 = form2( _trial=Wh, _test=VhTensor2 );
    a12 += integrate( _range=elements(mesh), _expr=trans(idt(w))*div(uTensor2));
    a12.close();
    auto a12v = a12(uTensor2,w);
    //    BOOST_CHECK_CLOSE( a12v, 4*area, 1e-11 );
    cout << "a12(ones, ones) = " << a12v << std::endl;

    cout << "a12 works" << std::endl;

    // a13

    // mesh for the Lagrange multiplier
    typedef Simplex<1,1,2> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    using Mh_t = Pdhv_type<face_mesh_type, 2>;
    using Mh_ptr_t = Pdhv_ptrtype<face_mesh_type, 2>;

    auto face_mesh = createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
    Mh_ptr_t Mh = Pdhv<2>( face_mesh, true );
    auto l = Mh->element();
    l.on( _range=elements(face_mesh), _expr=ones<2,1>());

    uTensor2.on(_range=elements(mesh),_expr= eye<2>() );

    auto a13a = form2( _trial=Mh, _test=VhTensor2 );
    a13a += integrate(_range=internalfaces(mesh),
                      _expr=-( trans(idt(l))*leftface(id(uTensor2)*N())+
                               trans(idt(l))*rightface(id(uTensor2)*N())) );
    a13a.close();
    auto a13a_v = a13a(uTensor2,l);
    BOOST_CHECK_SMALL( a13a_v, 1e-11 );
    cout << "a13a(eye, ones) = " << a13a_v << std::endl;

    l.on( _range=elements(face_mesh), _expr=oneX());

    auto a13b1 = form2( _trial=Mh, _test=VhTensor2 );
    a13b1 += integrate(_range=boundaryfaces(mesh),
                       _expr=-( trans(idt(l))*(id(uTensor2)*N())));
    a13b1.close();
    auto a13b1_v = a13b1(uTensor2,l);
    BOOST_CHECK_SMALL( a13b1_v, 1e-11 );
    cout << "a13b1(eye, oneX) = " << a13b1_v << std::endl;

    l.on( _range=elements(face_mesh), _expr=oneY());

    auto a13b2 = form2( _trial=Mh, _test=VhTensor2 );
    a13b2 += integrate(_range=boundaryfaces(mesh),
                       _expr=-( trans(idt(l))*(id(uTensor2)*N())));
    a13b2.close();
    auto a13b2_v = a13b2(uTensor2,l);
    BOOST_CHECK_SMALL( a13b2_v, 1e-11 );
    cout << "a13b2(eye, oneY) = " << a13b2_v << std::endl;

    cout << "a13 works" << std::endl;

    // a21

    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Px(),Py(),Px(),Py() ) );
    w.on(_range=elements(mesh), _expr=ones<2,1>());

    auto a21 = form2( _trial=VhTensor2, _test=Wh );
    a21 += integrate( _range=elements(mesh), _expr=trans(id(w))*divt(uTensor2));
    a21.close();
    auto a21v = a21(w,uTensor2);
    cout << "a21(ones,ones) = " << a12v << std::endl;

    cout << "a21 works" << std::endl;

    // a22

    auto a22 = form2( _trial=Wh, _test=Wh );

    a22 += integrate(_range=internalfaces(mesh),
                     _expr=
                     -( leftfacet( pow(h(),0)*trans(idt(w)))*leftface(id(w)) +
                       rightfacet( pow(h(),0)*trans(idt(w)))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(mesh),
                     _expr=-(pow(h(),0)*trans(idt(w))*id(w)));

    auto I = integrate( _range=internalfaces(mesh), _expr=cst(4.)).evaluate()(0,0)
        +integrate(_range=boundaryfaces(mesh),
                   _expr=cst(2.)).evaluate()(0,0);

    auto a22v = a22(w,w);
    BOOST_CHECK_CLOSE( a22v, -I, 1e-11 );
    cout << "a22(ones, ones) = " << a22v << std::endl;

    cout << "a22 works fine" << std::endl;

    // a23

    l.on( _range=elements(face_mesh), _expr=ones<2,1>());

    auto a23 = form2( _trial=Mh, _test=Wh );
    a23 += integrate(_range=internalfaces(mesh),
                     _expr=trans(idt(l)) *
                     ( leftface( pow(h(),0)*id(w) )+
                       rightface( pow(h(),0)*id(w) )));
    a23 += integrate(_range=boundaryfaces(mesh),
                     _expr=trans(idt(l)) * pow(h(),0)*id(w) );

    auto a23v = a23(w,l);
    cout << "a23(ones, ones) = " << a23v << std::endl;

    cout << "a23 works fine" << std::endl;

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
    test_tensor2<Pdhms_type<Mesh<Simplex<2>>,1>>();
}

BOOST_AUTO_TEST_SUITE_END()
