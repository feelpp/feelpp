#define BOOST_TEST_MODULE test_element_component
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( element_component )

template<int O=1>
std::shared_ptr<Mesh<Simplex<2,O,2>>>
getReferenceTriangleMesh()
{
    GeoTool::Node x1( -1,-1 );
    GeoTool::Node x2( 1,-1 );
    GeoTool::Node x3( -1,1 );
    GeoTool::Triangle R( 3,"OMEGA",x1,x2,x3 );
    R.setMarker(_type="line",_name="Neumann",_marker3=true);
    R.setMarker(_type="line",_name="Dirichlet",_marker1=true,_marker2=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    return R.createMesh(_mesh=new Mesh<Simplex<2,O,2>>(Environment::worldCommSeqPtr()),
                        _name= (boost::format("domainRef_p%1%")%Environment::worldComm().rank()).str() );
}

template<int O=1>
std::shared_ptr<Mesh<Simplex<2,O,2>>>
getTriangleMesh()
{
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Node x3( 0,1 );
    GeoTool::Triangle R( doption(_name="gmsh.hsize"),"OMEGA",x1,x2,x3 );
    R.setMarker(_type="line",_name="Neumann",_marker3=true);
    R.setMarker(_type="line",_name="Dirichlet",_marker1=true,_marker2=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    return R.createMesh(_mesh=new Mesh<Simplex<2,O,2>>,_name= "domain" );
}

BOOST_AUTO_TEST_CASE( element_component_vectorial )
{
    auto mesh = getTriangleMesh();

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
test_tensor2symm_basic(std::string const& name, std::shared_ptr<SpaceT> const& Vh )
{
    auto mesh = Vh->mesh();
    auto uTensor2 = Vh->element();
    auto uxx = uTensor2.comp( Component::X,Component::X );
    auto uxy = uTensor2.comp( Component::X,Component::Y );
    auto uyx = uTensor2.comp( Component::Y,Component::X );
    auto uyy = uTensor2.comp( Component::Y,Component::Y );
    uxx.on(_range=elements(mesh),_expr=cst(1.));
    uxy.on(_range=elements(mesh),_expr=cst(2.));
    uyy.on(_range=elements(mesh),_expr=cst(3.));

    BOOST_CHECK_CLOSE( uxx.max(), 1., 1e-12 );
    BOOST_CHECK_CLOSE( uxx.min(), 1., 1e-12 );
    BOOST_CHECK_CLOSE( uxy.max(), 2., 1e-12 );
    BOOST_CHECK_CLOSE( uxy.min(), 2., 1e-12 );
    BOOST_CHECK_CLOSE( uyx.max(), 2., 1e-12 );
    BOOST_CHECK_CLOSE( uyx.min(), 2., 1e-12 );
    BOOST_CHECK_CLOSE( uyy.max(), 3., 1e-12 );
    BOOST_CHECK_CLOSE( uyy.min(), 3., 1e-12 );

    uyx.on(_range=elements(mesh),_expr=cst(5.));
    BOOST_CHECK_CLOSE( uxy.max(), 5., 1e-12 );
    BOOST_CHECK_CLOSE( uxy.min(), 5., 1e-12 );
    BOOST_CHECK_CLOSE( uyx.max(), 5., 1e-12 );
    BOOST_CHECK_CLOSE( uyx.min(), 5., 1e-12 );
}

BOOST_AUTO_TEST_CASE( element_component_tensor2symm_continuous_basic )
{
    auto m1 = getReferenceTriangleMesh();
    test_tensor2symm_basic( "tensor2_cs_basic", Pchms<2>(  m1 ) );
    auto m2 = getTriangleMesh();
    test_tensor2symm_basic( "tensor2_cs_basic", Pchms<2>(  m2, false ) );
    test_tensor2symm_basic( "tensor2_cs_basic_edt", Pchms<2>(  m2, true ) );
}

template</*typename MeshT,*/ typename SpaceT>
void
test_tensor2(std::string const& name, std::shared_ptr<typename SpaceT::mesh_type/*MeshT*/> const& mesh )
{
    using Feel::cout;
    using mesh_type = typename SpaceT::mesh_type;
#if 0
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name = "domain",
                                              _shape = "hypercube",
                                              _usenames = true,
                                              _dim = 2,
                                              _h = doption(_name="gmsh.hsize"),
                                              _xmin=0,_xmax=4,
                                              _ymin=0,_ymax=1 ) );
#endif
    auto VhTensor2 = SpaceT::New( _mesh=mesh, _worldscomm=makeWorldsComm( 1,mesh->worldCommPtr() ),
                                  _extended_doftable=std::vector<bool>( 1,true ));
    auto uTensor2 = VhTensor2->element();
    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( cst(1.),cst(2.),cst(3.),cst(4.) ) );
    auto uxx = uTensor2.comp( Component::X,Component::X );
    auto uxy = uTensor2.comp( Component::X,Component::Y );
    auto uyx = uTensor2.comp( Component::Y,Component::X );
    auto uyy = uTensor2.comp( Component::Y,Component::Y );

    double sxxRef = integrate( _range=elements( mesh ), _expr=cst(1.) ).evaluate(parallelEvaluation)( 0,0 );
    double sxx1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneX(),oneX() ) ).evaluate(parallelEvaluation)( 0,0 );
    double sxx2 = integrate( _range=elements( mesh ), _expr=idv(uxx) ).evaluate(parallelEvaluation)( 0,0 );
    BOOST_CHECK_CLOSE( sxx1, sxxRef, 1e-12 );
    BOOST_CHECK_CLOSE( sxx2, sxxRef, 1e-12 );
    double sxyRef = 0;
    if ( SpaceT::is_tensor2symm )
        sxyRef = integrate( _range=elements( mesh ), _expr=cst(3.) ).evaluate(parallelEvaluation)( 0,0 );
    else
        sxyRef = integrate( _range=elements( mesh ), _expr=cst(2.) ).evaluate(parallelEvaluation)( 0,0 );
    double sxy1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneY(),oneX() ) ).evaluate(parallelEvaluation)( 0,0 );
    double sxy2 = integrate( _range=elements( mesh ), _expr=idv(uxy) ).evaluate(parallelEvaluation)( 0,0 );
    BOOST_CHECK_CLOSE( sxy1, sxyRef, 1e-12 );
    BOOST_CHECK_CLOSE( sxy2, sxyRef, 1e-12 );
    double syxRef = 0;
    syxRef = integrate( _range=elements( mesh ), _expr=cst(3.) ).evaluate(parallelEvaluation)( 0,0 );

    double syx1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneX(),oneY() ) ).evaluate(parallelEvaluation)( 0,0 );
    double syx2 = integrate( _range=elements( mesh ), _expr=idv(uyx) ).evaluate(parallelEvaluation)( 0,0 );
    BOOST_CHECK_CLOSE( syx1, syxRef, 1e-12 );
    BOOST_CHECK_CLOSE( syx2, syxRef, 1e-12 );
    double syyRef = integrate( _range=elements( mesh ), _expr=cst(4.) ).evaluate(parallelEvaluation)( 0,0 );
    double syy1 = integrate( _range=elements( mesh ), _expr=inner( idv( uTensor2 )*oneY(),oneY() ) ).evaluate(parallelEvaluation)( 0,0 );
    double syy2 = integrate( _range=elements( mesh ), _expr=idv(uyy) ).evaluate(parallelEvaluation)( 0,0 );
    BOOST_CHECK_CLOSE( syy1, syyRef, 1e-12 );
    BOOST_CHECK_CLOSE( syy2, syyRef, 1e-12 );

    auto sfull = integrate( _range=elements( mesh ), _expr=idv(uTensor2) ).evaluate(parallelEvaluation);
    BOOST_CHECK_CLOSE( sfull(0,0),sxxRef, 1e-12 );
    BOOST_CHECK_CLOSE( sfull(0,1),sxyRef, 1e-12 );
    BOOST_CHECK_CLOSE( sfull(1,0),syxRef, 1e-12 );
    BOOST_CHECK_CLOSE( sfull(1,1),syyRef, 1e-12 );

    double area = integrate( _range=elements( mesh ), _expr=cst(1.) ).evaluate(parallelEvaluation)(0,0);
    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Px(),cst(0.),cst(3.),2.3*Py() ) );
    auto sdiv = integrate( _range=elements( mesh ), _expr=divv(uTensor2) ).evaluate(parallelEvaluation);
    BOOST_CHECK_CLOSE( sdiv(0,0),area, 1e-12 );
    BOOST_CHECK_CLOSE( sdiv(1,0),2.3*area, 1e-12 );

    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Py(),3.4*Px(),4.5*Py(),2.3*Py() ) );
    sdiv = integrate( _range=elements( mesh ), _expr=divv(uTensor2) ).evaluate(parallelEvaluation);
    if ( SpaceT::is_tensor2symm )
        BOOST_CHECK_CLOSE( sdiv(0,0),4.5*area, 1e-12 );
    else
        BOOST_CHECK_SMALL( sdiv(0,0), 1e-12 );
    BOOST_CHECK_CLOSE( sdiv(1,0),2.3*area, 1e-12 );

    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Px(),3.4*Py(),3.8*Px(),2.3*Py()+Px() ) );
    sdiv = integrate( _range=elements( mesh ), _expr=divv(uTensor2) ).evaluate(parallelEvaluation);
    if ( SpaceT::is_tensor2symm )
        BOOST_CHECK_CLOSE( sdiv(0,0),area, 1e-12 );
    else
        BOOST_CHECK_CLOSE( sdiv(0,0),4.4*area, 1e-12 );
    BOOST_CHECK_CLOSE( sdiv(1,0),6.1*area, 1e-12 );

#if 0
    uTensor2.on(_range=elements(mesh),_expr= mat<2,2>( Px()*Py(),Py()+Px(),Py(),2.3*Py()+Px() ) );
    sdiv = integrate( _range=elements( mesh ), _expr=divv(uTensor2) ).evaluate(parallelEvaluation);
    BOOST_CHECK_CLOSE( sdiv(0,0),2*area, 1e-12 );
    BOOST_CHECK_CLOSE( sdiv(1,0),2.3*area, 1e-12 );
#endif

    /*
     * Test bilinear forms for HDG linear elasticity
     */

    auto Wh = Pdhv<1>( mesh, true );
    auto w = Wh->element();
    auto w1 = Wh->element();

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
    //auto f = expr<2,2>( soption( "functions.f" ) );
    //auto f = expr<2,2>( "{2*x,2,3,2*y}:x:y" );
    auto f = expr<2,2>( "{2*x,5*y,3*x,2*y}:x:y" );

    uTensor2.on(_range=elements(mesh),_expr= f );
    w1.on(_range=elements(mesh),_expr= divv(uTensor2));

    uxx = uTensor2.comp( Component::X,Component::X );
    uxy = uTensor2.comp( Component::X,Component::Y );
    uyx = uTensor2.comp( Component::Y,Component::X );
    uyy = uTensor2.comp( Component::Y,Component::Y );
    auto w1x=w1.comp( Component::X );
    auto w1y=w1.comp( Component::Y );

#if 0
    std::cout << "Tensor2: " << uTensor2 << std::endl;
    std::cout << "uxx: " << uxx << std::endl;
    std::cout << "uyy: " << uyy << std::endl;
    std::cout << "uxy: " << uxy << std::endl;
    std::cout << "uyx: " << uyx << std::endl;
    std::cout << "w1x: " << w1x << std::endl;
    std::cout << "w1y: " << w1y << std::endl;

    auto e  = exporter(_mesh=mesh,_name=name);
    e->add( "u", uTensor2 );
    e->add( "w1", w1 );
    e->add( "uxx", uxx );
    e->add( "uyy", uyy );
    e->add( "uxy", uxy );
    e->add( "uyx", uyx );
    e->save();
#endif
    auto a12 = form2( _trial=Wh, _test=VhTensor2 );
    a12 = integrate( _range=elements(mesh), _expr=trans(idt(w))*div(uTensor2));
    a12.close();
    auto a12v = a12(uTensor2,w);
    if ( SpaceT::is_tensor2symm )
        BOOST_CHECK_CLOSE( a12v, 7*area, 1e-11 );
    else
        BOOST_CHECK_CLOSE( a12v, 12*area, 1e-11 );
    double a12vv = integrate( _range=elements(mesh), _expr=trans(idv(w))*divv(uTensor2)).evaluate(parallelEvaluation)(0,0);
    double a12vvv = integrate( _range=elements(mesh), _expr=trans(idv(w))*(2*ones<2,1>())).evaluate(parallelEvaluation)(0,0);
    cout << "a12vv=" << a12vv << " a12v=" << a12v << " a12vvv=" << a12vvv << std::endl;
    BOOST_CHECK_CLOSE( a12v, a12vv, 1e-11 );
    cout << "a12(ones, ones) = " << a12v << std::endl;

    cout << "a12 works" << std::endl;

    // a13
    // mesh for the Lagrange multiplier
    //typedef Simplex<1,1,2> face_convex_type;
    typedef /*Mesh<face_convex_type>*/typename mesh_type::trace_mesh_type face_mesh_type;
    using Mh_t = Pdhv_type<face_mesh_type, 2>;
    using Mh_ptr_t = Pdhv_ptrtype<face_mesh_type, 2>;

    auto face_mesh = createSubmesh( _mesh=mesh, _range=faces(mesh),_update=0 );
    Mh_ptr_t Mh = Pdhv<2>( face_mesh, true );
    auto l = Mh->element();
    l.on( _range=elements(face_mesh), _expr=ones<2,1>());

    uTensor2.on(_range=elements(mesh),_expr= eye<2>() );

    auto a13aMat = backend()->newMatrix(_test=VhTensor2,_trial=Mh,_buildGraphWithTranspose=true);
    auto a13a = form2( _trial=Mh, _test=VhTensor2,_matrix=a13aMat );
    a13a += integrate(_range=internalfaces(mesh),
                      _expr=-( trans(idt(l))*leftface(id(uTensor2)*N())+
                               trans(idt(l))*rightface(id(uTensor2)*N())) );
    a13a.close();
    auto a13a_v = a13a(uTensor2,l);

    BOOST_CHECK_SMALL( a13a_v, 1e-11 );
    cout << "a13a(eye, ones) = " << a13a_v << std::endl;
#if 1
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

    auto I = integrate( _range=internalfaces(mesh), _expr=cst(4.)).evaluate(parallelEvaluation)(0,0)
        +integrate(_range=boundaryfaces(mesh),
                   _expr=cst(2.)).evaluate(parallelEvaluation)(0,0);

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
    a23.close();
    auto a23v = a23(w,l);
    cout << "a23(ones, ones) = " << a23v << std::endl;

    auto a23b = form2( _trial=Mh, _test=Wh );
    a23b += integrate(_range=boundaryfaces(mesh),
                      _expr=trans(idt(l)) * pow(h(),0)*id(w) );
    a23b.close();
    auto a23b_v = a23b(w,l);
    cout << "a23b(ones, ones) = " << a23b_v << std::endl;

    cout << "a23 works fine" << std::endl;

    // a31

    l.on( _range=elements(face_mesh), _expr=ones<2,1>());

    uTensor2.on(_range=elements(mesh),_expr= eye<2>() );

    auto a31a = form2( _trial=VhTensor2, _test=Mh );
    a31a += integrate(_range=internalfaces(mesh),
                      _expr=-trans(id(l))*(leftfacet(idt(uTensor2)*N())+
                                           rightfacet(idt(uTensor2)*N())) );
    a31a.close();
    auto a31a_v = a31a(l,uTensor2);
    BOOST_CHECK_SMALL( a31a_v, 1e-11 );
    cout << "a31a(ones, eye) = " << a31a_v << std::endl;

    auto a31b = form2( _trial=VhTensor2, _test=Mh );
    a31b += integrate(_range=markedfaces(mesh,"Neumann"),
                      _expr=trans(id(l))*(idt(uTensor2)*N()));
    a31b.close();
    double a31b_v_res = integrate(_range=markedfaces(mesh,"Neumann"),
                                  _expr=trans(one())*N()).evaluate(parallelEvaluation)(0,0);
    auto a31b_v = a31b(l,uTensor2);
    BOOST_CHECK_SMALL( a31b_v - a31b_v_res, 1e-11 );
    cout << "a31b(ones, eye) = " << a31b_v << std::endl;

    cout << "a31 works fine" << std::endl;

    // a32

    auto a32a = form2( _trial=Wh, _test=Mh );

    a32a += integrate(_range=internalfaces(mesh),
                     _expr=-trans(id(l)) * (leftfacet( pow(h(),0)*idt(w) )+
                                            rightfacet( pow(h(),0)*idt(w) )));
    a32a.close();
    auto a32a_v = a32a(l,w);
    BOOST_CHECK_CLOSE( a32a_v, -a23v, 1e-11 );
    cout << "a32a(ones, ones) = " << a32a_v << std::endl;

    auto a32b = form2( _trial=Wh, _test=Mh );

    a32b += integrate(_range=markedfaces(mesh,"Neumann"),
                      _expr=-trans(id(l)) * ( pow(h(),0)*idt(w) ) );
    a32b.close();
    auto a32b_v = a32b(l, w);
    cout << "a32b(ones, ones) = " << a32b_v << std::endl;

    cout << "a32 works fine" << std::endl;

    // a33

    auto I1 = integrate( _range=internalfaces(mesh), _expr=cst(2.)).evaluate(parallelEvaluation)(0,0);
    auto I2 = integrate( _range=elements(face_mesh), _expr=cst(2.)).evaluate(parallelEvaluation)(0,0);

    auto a33a = form2(_trial=Mh, _test=Mh );
    a33a += integrate(_range=internalfaces(mesh),
                      _expr=leftfacet(trans(idt(l))) * leftface(id(l)));//( leftface( pow(h(),0) )+
    //                                                      rightface( pow(h(),0) )));
    a33a.close();
    auto a33a_v = a33a(l,l);
    BOOST_CHECK_CLOSE( a33a_v, I1, 1e-11 );

    a33a = integrate(_range=elements(face_mesh),
                     _expr=trans(idt(l)) * id(l));
    a33a.close();
    a33a_v = a33a(l,l);
    BOOST_CHECK_CLOSE( a33a_v, I2, 1e-11 );

    auto a33b = form2(_trial=Mh, _test=Mh );
    a33b += integrate(_range=markedfaces(mesh,"Neumann"),
                      _expr=trans(idt(l)) * id(l) * ( pow(h(),0) ) );
    a33b.close();
    auto a33b_v = a33b(l,l);

    auto a33c = form2(_trial=Mh, _test=Mh);
    a33c += integrate(_range=markedfaces(mesh,"Dirichlet"),
                      _expr=trans(idt(l)) * id(l) );
    a33c.close();
    auto a33c_v = a33c(l,l);
    double lenghtDirichlet = integrate(_range=markedfaces(mesh,"Dirichlet"),
                                     _expr=cst(1.) ).evaluate(parallelEvaluation)(0,0);
    double lenghtNeumann = integrate(_range=markedfaces(mesh,"Neumann"),
                                     _expr=cst(1.) ).evaluate(parallelEvaluation)(0,0);
    BOOST_CHECK_CLOSE(a33b_v + a33c_v, 2*lenghtDirichlet+2*lenghtNeumann, 1e-11);

    cout << "a33 works fine" << std::endl;
#endif
}

BOOST_AUTO_TEST_CASE( element_component_tensor2_continuous_ref )
{
    BOOST_TEST_MESSAGE( "element_component_tensor2_continuous_ref starts" );
    auto m1 = getReferenceTriangleMesh();
    test_tensor2<Pchm_type<Mesh<Simplex<2>>,1>>("tensor2_c_ref",m1);
    BOOST_TEST_MESSAGE( "element_component_tensor2_continuous_ref ends" );
}

BOOST_AUTO_TEST_CASE( element_component_tensor2_continuous )
{
    BOOST_TEST_MESSAGE( "element_component_tensor2_continuous starts" );
    auto m2 = getTriangleMesh();
    test_tensor2<Pchm_type<Mesh<Simplex<2>>,1>>("tensor2_c",m2);
    BOOST_TEST_MESSAGE( "element_component_tensor2_continuous ends" );
}

BOOST_AUTO_TEST_CASE( element_component_tensor2_discontinuous )
{
    auto m1 = getReferenceTriangleMesh();
    test_tensor2<Pdhm_type<Mesh<Simplex<2>>,2>>("tensor2_c",m1);
    auto m2 = getTriangleMesh();
    test_tensor2<Pdhm_type<Mesh<Simplex<2>>,2>>("tensor2_d",m2);
}

BOOST_AUTO_TEST_CASE( element_component_tensor2symm_continuous )
{
    auto m1 = getReferenceTriangleMesh();
    test_tensor2<Pchms_type<Mesh<Simplex<2>>,2>>("tensor2_cs",m1);
    auto m2 = getTriangleMesh();
    test_tensor2<Pchms_type<Mesh<Simplex<2>>,2>>("tensor2_cs",m2);
}

BOOST_AUTO_TEST_CASE( element_component_tensor2symm_discontinuous )
{
    auto m1 = getReferenceTriangleMesh();
    test_tensor2<Pdhms_type<Mesh<Simplex<2>>,2>>("tensor2_ds",m1);
    auto m2 = getTriangleMesh();
    test_tensor2<Pdhms_type<Mesh<Simplex<2>>,2>>("tensor2_ds",m2);
}

BOOST_AUTO_TEST_SUITE_END()
