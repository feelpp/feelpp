
#define BOOST_TEST_MODULE operatorlagrangep1 testsuite
#include <feel/feelcore/testsuite.hpp>

#include <boost/mp11/utility.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>


using namespace Feel;

template <int OrderLagP1,template<class, int, class> class PointSetType,typename MeshType>
void run_test_oplagp1( std::shared_ptr<MeshType> meshBase )
{
    auto VhBase = Pch<OrderLagP1,double,PointSetType>( meshBase );
    auto opLagP1 = lagrangeP1( _space = VhBase );
    auto mesh = opLagP1->mesh();
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements( mesh ),
                   _expr=inner( gradt( u ), grad( u ) ) );
    auto l = form1( _test = Vh );
    auto solution = Px()*Py();
    a += on( _range=boundaryfaces( mesh ), _rhs=l, _element=u, _expr=solution );
    a.solve( _rhs=l, _solution=u, _rebuild=true );

    double l2 = normL2( _range=elements( mesh ), _expr=idv( u )-solution );
    BOOST_TEST_MESSAGE( "l2 Norm" + std::to_string( l2 ) );
    BOOST_CHECK_SMALL( l2, 1e-7 );
}

FEELPP_ENVIRONMENT_NO_OPTIONS

using order_geo_types = boost::mp11::mp_list_c<int, 1, 2>;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_oplagp1_2d, T, order_geo_types )
{
    static const int orderGeo = T::value;
    BOOST_TEST_MESSAGE( "test_oplagp1_2d_geo" + std::to_string(orderGeo) );
    typedef Mesh<Simplex<2,orderGeo>> mesh_type;
    auto meshBase = loadMesh( _mesh = new mesh_type );

    BOOST_TEST_MESSAGE( "run P2 PointSetEquiSpaced" );
    run_test_oplagp1<2,PointSetEquiSpaced>( meshBase );
    BOOST_TEST_MESSAGE( "run P3 PointSetEquiSpaced" );
    run_test_oplagp1<3,PointSetEquiSpaced>( meshBase );
    BOOST_TEST_MESSAGE( "run P4 PointSetEquiSpaced" );
    run_test_oplagp1<4,PointSetEquiSpaced>( meshBase );
    BOOST_TEST_MESSAGE( "run P5 PointSetEquiSpaced" );
    run_test_oplagp1<5,PointSetEquiSpaced>( meshBase );

    BOOST_TEST_MESSAGE( "run P3 PointSetFekete" );
    run_test_oplagp1<3,PointSetFekete>( meshBase );
    BOOST_TEST_MESSAGE( "run P4 PointSetFekete" );
    run_test_oplagp1<4,PointSetFekete>( meshBase );
    BOOST_TEST_MESSAGE( "run P5 PointSetFekete" );
    run_test_oplagp1<5,PointSetFekete>( meshBase );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_oplagp1_3d, T, order_geo_types )
{
    static const int orderGeo = T::value;
    BOOST_TEST_MESSAGE( "test_oplagp1_3d_geo" + std::to_string(orderGeo) );
    typedef Mesh<Simplex<3,orderGeo>> mesh_type;
    auto meshBase = loadMesh( _mesh = new mesh_type );
    BOOST_TEST_MESSAGE( "run P2 PointSetEquiSpaced" );
    run_test_oplagp1<2,PointSetEquiSpaced>( meshBase );
    BOOST_TEST_MESSAGE( "run P3 PointSetEquiSpaced" );
    run_test_oplagp1<3,PointSetEquiSpaced>( meshBase );

    BOOST_TEST_MESSAGE( "run P3 PointSetWarpBlend" );
    run_test_oplagp1<3,PointSetWarpBlend>( meshBase );
}

BOOST_AUTO_TEST_CASE( test_operatorinterpolation )
{
    BOOST_TEST_MESSAGE( "test_operatorinterpolation" );

    static const int OrderGeo = 1;
    typedef Mesh<Simplex<2, OrderGeo, 2>> mesh_type;
    double meshSize = doption(_name="gmsh.hsize" );
    GeoTool::Node x1( 0, 0 );
    GeoTool::Node x2( 0.6, 0 );
    GeoTool::Circle C( meshSize, "OMEGA", x1, x2 );
    C.setMarker( _type = "line", _name = "Sortie", _markerAll = true );
    C.setMarker( _type = "surface", _name = "OmegaFluide", _markerAll = true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name = "test2dOpLagrangeP1_domain" + mesh_type::shape_type::name() );

    auto Xh = Pchv<3, PointSetFekete>( mesh );
    auto exprProj = vec( cos( M_PI * Px() ), sin( M_PI * Py() ) );
    auto u = vf::project( _space = Xh,
                          _range = elements( mesh ),
                          _expr = exprProj );

    auto mybackend = backend( _rebuild = true );

    auto opLagP1 = lagrangeP1( _space = Xh );
    auto meshLagP1 = opLagP1->mesh();

    auto XhLagP1 = Pchv<1>( meshLagP1 );
    auto uLagP1interp = XhLagP1->element();

    auto uLagP1proj = vf::project( _space = XhLagP1,
                                   _range = elements( meshLagP1 ),
                                   _expr = exprProj );

    auto opI = opInterpolation( _domainSpace = Xh,
                                _imageSpace = XhLagP1,
                                _range = elements( meshLagP1 ) );
    opI->apply( u, uLagP1interp );

    auto s1 = integrate( _range = elements( meshLagP1 ),
                         _expr = inner( idv( uLagP1interp ) - idv( uLagP1proj ), idv( uLagP1interp ) - idv( uLagP1proj ) ) )
                  .evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( s1, 1e-6 );

    BOOST_TEST_MESSAGE( "test_operatorinterpolation done" );
}
