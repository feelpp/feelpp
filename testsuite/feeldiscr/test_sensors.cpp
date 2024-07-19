#define BOOST_TEST_MODULE sensors testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/sensors.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sensorssuite )

typedef boost::mpl::list<boost::mpl::integral_c<int, 2>, boost::mpl::integral_c<int, 3>> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_create, T, test_types)
{
    using mesh_type = Mesh<Simplex<T::value> >;
    using space_type = Pch_type<mesh_type, 2>;
    using sensormap_type = SensorMap<space_type>;

    auto filename = "markerhdf5_"+std::to_string(T::value)+"D.geo";
    auto mesh = loadMesh( _mesh=new mesh_type, _filename=filename );
    auto Xh = Pch<2>(mesh);
    auto i = std::ifstream(Environment::expand("$cfgdir/sensorsdesc.json"));
    auto sm = sensormap_type(Xh, json::parse(i));

    for( auto e : std::vector{ "x:x:y", "y:x:y", "x+y:x:y", "(x+y)*(x+y):x:y" } )
    {
        BOOST_TEST_MESSAGE(fmt::format("[test_create] check expression: {}\n", e ) );
        auto ex = expr(e);
        auto u = Xh->element(ex);
        for( auto const& [name, f] : sm )
        {
            if( auto ff = std::dynamic_pointer_cast<SensorPointwise<space_type>>(f) )
            {
                auto n = ff->position();
                auto vv = ex.evaluate({{"x", n(0)}, {"y", n(1)}})(0,0);
                auto v = (*ff)(u);
                BOOST_TEST_MESSAGE(fmt::format("[test_create] sensor::pointwise: exact : {} SensorPointwise : {}\n",vv,v));
                BOOST_CHECK_SMALL(std::abs(v-vv), 1e-10 );
            } else if( auto ff = std::dynamic_pointer_cast<SensorGaussian<space_type>>(f) )
            {
                auto n = ff->position();
                auto vv = ex.evaluate({{"x", n(0)}, {"y", n(1)}})(0,0);
                auto v = (*ff)(u);
                BOOST_CHECK_SMALL(std::abs(v-vv), 1e-2 );
            } else if( auto ff = std::dynamic_pointer_cast<SensorSurface<space_type>>(f) )
            {
                auto m = ff->markers();
                auto vv = mean(_range=markedfaces(mesh,m), _expr=ex)(0,0);
                auto v = (*ff)(u);
                BOOST_TEST_MESSAGE(fmt::format("[test_create] sensor::sensorsurface: markers: {}, exact : {} SensorPointwise : {}\n",m,vv,v));
                BOOST_CHECK_SMALL(std::abs(v-vv), 1e-10 );
            }
        }
    }

    auto ofa = std::ofstream("sensors.db");
    boost::archive::binary_oarchive oa( ofa );
    oa << sm;

    json j = sm.to_json();
    BOOST_CHECK_EQUAL(j.size(), 3);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_load, T, test_types)
{
    using mesh_type = Mesh<Simplex<T::value> >;
    using space_type = Pch_type<mesh_type, 2>;
    using sensormap_type = SensorMap<space_type>;

    auto filename = "markerhdf5_"+std::to_string(T::value)+"D.geo";
    auto mesh = loadMesh( _mesh=new mesh_type, _filename=filename );
    auto Xh = Pch<2>(mesh);
    auto sm = sensormap_type(Xh);

    auto ifa = std::ifstream("sensors.db");
    boost::archive::binary_iarchive ia( ifa );
    ia >> sm;

    std::cout << "loaded with " << sm.size() << " sensors" << std::endl;

    for( auto e : std::vector{ "x:x:y", "y:x:y", "x+y:x:y", "(x+y)*(x+y):x:y" } )
    {
        BOOST_TEST_MESSAGE(fmt::format("[test_load] check expression: {}\n", e ) );
#if 0 
        auto ex = expr(e);
       
        auto u = Xh->element(ex);
        for( auto const& [name, f] : sm )
        {
            if( auto ff = std::dynamic_pointer_cast<SensorPointwise<space_type>>(f) )
            {
                auto n = ff->position();
                auto vv = ex.evaluate({{"x", n(0)}, {"y", n(1)}})(0,0);
                auto v = (*ff)(u);
                BOOST_TEST_MESSAGE(fmt::format("-- sensor::pointwise: exact : {} SensorPointwise : {}",vv,v));
                BOOST_CHECK_SMALL(std::abs(v-vv), 1e-10 );
            } else if( auto ff = std::dynamic_pointer_cast<SensorGaussian<space_type>>(f) )
            {
                auto n = ff->position();
                auto vv = ex.evaluate({{"x", n(0)}, {"y", n(1)}})(0,0);
                auto v = (*ff)(u);
                BOOST_CHECK_SMALL(std::abs(v-vv), 1e-2 );
            } else if( auto ff = std::dynamic_pointer_cast<SensorSurface<space_type>>(f) )
            {
                auto m = ff->markers();
                auto vv = mean(_range=markedfaces(mesh,m), _expr=ex)(0,0);
                auto v = (*ff)(u);
                BOOST_CHECK_SMALL(std::abs(v-vv), 1e-10 );
            }
        }
#endif        
    }
    json j = sm.to_json();
    BOOST_CHECK_EQUAL(j.size(), 3);
}

BOOST_AUTO_TEST_SUITE_END()

