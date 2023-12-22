#define BOOST_TEST_MODULE bvh_tests
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/bvh.hpp>
using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "test_bvh" ,
                     "test_bvh" ,
                     "0.2",
                     "nD(n=2,3) test bvh",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2023 Feel++ Consortium" );

    about.addAuthor( "Luca Berti", "developer", "luca.berti@cemosis.fr", "" );
    return about;

}

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()
        ( "mesh2D.filename", po::value<std::string>(), "mesh2D.filename" )
        ( "mesh3D.filename", po::value<std::string>(), "mesh3D.filename" )
        ;
    return opts;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( bvh_intersection_tests )


template <typename BvhType,typename RayIntersectionResultType>
void printRayIntersectionResults( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    if ( bvh->worldComm().isMasterRank() )
        BOOST_TEST_MESSAGE( "Number of intersection: " << rirs.size() );
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            BOOST_TEST_MESSAGE( " --  ProcessId: " << rir.processId() );
            BOOST_TEST_MESSAGE( " --  PrimitiveId: " << rir.primitiveId() );
            BOOST_TEST_MESSAGE( " --  Distance: " << rir.distance() );

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            BOOST_TEST_MESSAGE( " --  Mesh entity id: " << prim.meshEntity().id() );
            BOOST_TEST_MESSAGE( " --  Mesh entity barycenter: " << prim.meshEntity().barycenter() );
            BOOST_TEST_MESSAGE( " ---------------------------------" );
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}
template <typename RangeType>
void test2D( RangeType const& range )
{
    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    Eigen::VectorXd origin(2);
    origin<< 0.5,0.5;
    Eigen::VectorXd direction_perp_1(2);
    direction_perp_1 << 1.,0.;
    Eigen::VectorXd direction_perp_2(2);
    direction_perp_2 << 0.,1.;
    bvh_ray_type ray_1(origin,direction_perp_1);
    bvh_ray_type ray_2(origin,direction_perp_2);

    auto bvhInHouse = boundingVolumeHierarchy(_range=range,_kind="in-house");
    auto rayIntersectionResult1 = bvhInHouse->intersect(_ray=ray_1) ;
    BOOST_CHECK_MESSAGE(!rayIntersectionResult1.empty(), fmt::format("Intersection between ray1 and BVH tree has been found"));
    printRayIntersectionResults( bvhInHouse,rayIntersectionResult1 );
    auto rayIntersectionResult2 = bvhInHouse->intersect(_ray=ray_2) ;
    BOOST_CHECK_MESSAGE(!rayIntersectionResult2.empty(), fmt::format("Intersection between ray2 and BVH tree has been found"));
    printRayIntersectionResults( bvhInHouse,rayIntersectionResult2 );
}

template <typename RangeType>
void test3D( RangeType const& range )
{
    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    Eigen::VectorXd origin(3);
    origin<< 0.5,0.5,0.5;
    Eigen::VectorXd origin2(3);
    origin2<< -10.0,0.5,0.5;
    Eigen::VectorXd direction_perp_1(3);
    direction_perp_1 << 1.,0.,0.;
    Eigen::VectorXd direction_perp_2(3);
    direction_perp_2 << 0.,1.,0.;
    std::vector<bvh_ray_type> rays;
    rays.push_back( bvh_ray_type(origin,direction_perp_1) );
    rays.push_back( bvh_ray_type(origin,direction_perp_2) );
    rays.push_back( bvh_ray_type(origin2,direction_perp_1) );

    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
    for (int k=0;k<rays.size();++k)
    {
        raysDistributed.push_back( rays[k] );
    }

    auto bvhInHouse = boundingVolumeHierarchy(_range=range,_kind="in-house");
    auto bvhThirdParty = boundingVolumeHierarchy(_range=range,_kind="third-party");
    auto bvhThirdPartyLow = boundingVolumeHierarchy(_range=range,_kind="third-party",_quality=BVHEnum::Quality::Low);
    for ( auto bvhCurrent : {bvhInHouse.get(),bvhThirdParty.get(),bvhThirdPartyLow.get()} )
    {
        for ( auto const& ray : rays )
        {
            auto rayIntersectionResult1 = bvhCurrent->intersect(_ray=ray) ;
            BOOST_CHECK_MESSAGE(!rayIntersectionResult1.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            printRayIntersectionResults( bvhCurrent,rayIntersectionResult1 );
        }

        auto multiRayIntersectionResult = bvhCurrent->intersect(_ray=rays);
        for ( auto const& rayIntersectionResult : multiRayIntersectionResult)
        {
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            printRayIntersectionResults( bvhCurrent,rayIntersectionResult );
        }

        auto multiRayDistributedIntersectionResult = bvhCurrent->intersect(_ray=raysDistributed);
        for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            printRayIntersectionResults( bvhCurrent,rayIntersectionResult );
        }
    }
}

BOOST_AUTO_TEST_CASE( intersection_bvh_2D )
{
    using mesh_type = Mesh<Simplex<2,1,2>>;
    auto mesh = loadMesh(_mesh = new mesh_type, _filename=soption(_name="mesh2D.filename" ) );

    auto submesh = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,{"RequiredBoundaryOfRequiredElements"}));
    test2D( elements(submesh) );
}

BOOST_AUTO_TEST_CASE( intersection_bvh_3D )
{
    using mesh_type = Mesh<Simplex<3,1,3>>;
    auto mesh = loadMesh(_mesh = new mesh_type, _filename=soption(_name="mesh3D.filename" ) );
    auto rangeFaces = markedfaces(mesh,{"CavityBottom","CavitySides","CavityTop"});
    auto submesh = createSubmesh(_mesh=mesh,_range=rangeFaces);

    test3D( rangeFaces );
    test3D( elements(submesh) );
}

BOOST_AUTO_TEST_SUITE_END()
