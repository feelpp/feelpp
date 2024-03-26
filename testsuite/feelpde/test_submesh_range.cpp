#define BOOST_TEST_MODULE test_laplacian
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
inline po::options_description
makeOptions()
{
    po::options_description opts( "test_submesh" );
    opts.add_options()( "mu", po::value<double>()->default_value( 1. ), "");
    opts.add_options()( "marker.name", po::value<std::string>()->default_value( "Omega1" ), "marker name" );
    return opts.add( Feel::feel_options() );
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "test_laplacian" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( test_laplacian )

typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_scalar, T, dim_types )
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<T::value>>, _filename=T::value==2?"bidomain_square.geo":"bidomain_cube.geo");
    double I1 = 0, I0 = 0, I2 = 0, I3 = 0;
    if ( T::value == 2 )
    {
        I1 = 4*0.4; // perimeter
        I2 = 0.4*0.4; // area
        I3 = 4*1+4*0.4; // perimenter external
    }
    if ( T::value == 3 )
    {
        I1 = 6*0.4*0.4;
        I2 = 0.4*0.4*0.4;
        I3 = 6*1+6*0.4*0.4;
    }
    auto Vh = Pch<1>( mesh, markedelements(mesh, "Omega1" ) );
    auto u = Vh->element();
    u.on( _range=elements( support(Vh) ), _expr=constant(1.) );
    auto I = integrate( _range=elements( support(Vh) ), _expr=idv(u) ).evaluate();
    std::cout << fmt::format( "meas Omega1 = {}", I ) << std::endl;
    BOOST_CHECK_CLOSE( I.norm(), I2, 1.e-10 );
    I = integrate( _range=boundaryfaces( support(Vh) ), _expr=idv(u) ).evaluate();
    std::cout << fmt::format( "meas boundary Omega1 = {}", I ) << std::endl;
    BOOST_CHECK_CLOSE( I.norm(), I1, 1.e-10 );

    auto Xh = Pch<1>( mesh, markedelements( mesh, "Omega2" ) );
    auto v = Xh->element();
    v.on( _range = elements( support( Xh ) ), _expr = constant( 1. ) );
    I = integrate( _range = elements( support( Xh ) ), _expr = idv( v )  ).evaluate();
    std::cout << fmt::format( "meas Omega2 = {}", I ) << std::endl;
    BOOST_CHECK_CLOSE( I.norm(), 1 - I2, 1.e-10 );
    I = integrate( _range=boundaryfaces( support(Xh) ), _expr=idv(v) ).evaluate();
    std::cout << fmt::format( "measure boundary Omega2 = {}", I ) << std::endl;
    BOOST_CHECK_CLOSE( I.norm(), I3, 1.e-10 );
}

using dim_types = boost::mpl::list< boost::mpl::int_<2>, boost::mpl::int_<3> >;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_vectorial, T, dim_types )
{
    auto mesh = loadMesh( _mesh = new Mesh<Simplex<T::value>>, _filename = T::value == 2 ? "bidomain_square.geo" : "bidomain_cube.geo" );
    double I1 = 0, I0 = 0, I2 = 0, I3 = 0;
    if ( T::value == 2 )
    {
        I1 = 4 * 0.4;         // perimeter
        I2 = 0.4 * 0.4;       // area
        I3 = 4 * 1 + 4 * 0.4; // perimenter external
    }
    if ( T::value == 3 )
    {
        I1 = 6 * 0.4 * 0.4;
        I2 = 0.4 * 0.4 * 0.4;
        I3 = 6 * 1 + 6 * 0.4 * 0.4;
    }
    auto Vh = Pchv<1>( mesh, markedelements( mesh, "Omega1" ) );
    auto u = Vh->element();
    u.on( _range = elements( support( Vh ) ), _expr = one() );
    auto I = integrate( _range = elements( support( Vh ) ), _expr = trans(idv( u ))*one()/T::value ).evaluate();
    std::cout << fmt::format( "meas Omega1 = {}", I ) << std::endl;
    BOOST_CHECK_CLOSE( I.norm(), I2, 1.e-10 );
    auto Iv = integrate( _range = boundaryfaces( support( Vh ) ), _expr = idv( u ) ).evaluate();
    std::cout << fmt::format( "meas boundary Omega1 = {}", Iv ) << std::endl;
    BOOST_CHECK_CLOSE( Iv.sum()/T::value, I1, 1.e-10 );
    I = integrate( _range = boundaryfaces( support( Vh ) ), _expr = trans(idv( u ))*N() ).evaluate();
    std::cout << fmt::format( "int 1 . N = {}", I ) << std::endl;
    BOOST_CHECK_SMALL( I.norm(), 1.e-10 );
    u.on( _range = elements( support( Vh ) ), _expr = Px()*oneX() );
    I = integrate( _range = boundaryfaces( support( Vh ) ), _expr = trans( idv( u ) ) * N() ).evaluate();
    std::cout << fmt::format( "int X . N = {}, exact={}", I, I2 ) << std::endl;
    BOOST_CHECK_CLOSE( I.norm(), I2, 1.e-10 );

    auto Xh = Pchv<1>( mesh, markedelements( mesh, "Omega2" ) );
    auto v = Xh->element();
    v.on( _range = elements( support( Xh ) ), _expr = one() );
    I = integrate( _range = elements( support( Xh ) ), _expr = trans( idv( v ) ) * one() / T::value ).evaluate();
    std::cout << fmt::format( "meas Omega2 = {}", I ) << std::endl;
    BOOST_CHECK_CLOSE( I.norm(), 1 - I2, 1.e-10 );
    Iv = integrate( _range = boundaryfaces( support( Xh ) ), _expr = idv( v ) ).evaluate();
    std::cout << fmt::format( "measure boundary Omega2 = {}", Iv ) << std::endl;
    BOOST_CHECK_CLOSE( Iv.sum()/T::value, I3, 1.e-10 );
    I = integrate( _range = boundaryfaces( support( Xh ) ), _expr = trans( idv( v ) ) * N() ).evaluate();
    std::cout << fmt::format( "int 1 . N = {}", I ) << std::endl;
    BOOST_CHECK_SMALL( I.norm(), 1.e-10 );

    v.on( _range = elements( support( Xh ) ), _expr = Px()*oneX() );
    I = integrate( _range = boundaryfaces( support( Xh ) ), _expr = trans( idv( v ) ) * N() ).evaluate();
    std::cout << fmt::format( "int X . N = {}, exact = {}", I, 1-I2 ) << std::endl;
    BOOST_CHECK_CLOSE( I.norm(), 1-I2, 1.e-10 );

}
BOOST_AUTO_TEST_SUITE_END()
