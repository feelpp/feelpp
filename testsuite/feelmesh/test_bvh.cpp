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
void printRayIntersectionResults( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs, std::vector<typename std::remove_pointer_t<BvhType>::vector_realdim_type  /*Eigen::Vector3d*/> const& pointIntersection )
{
    if ( bvh->worldComm().isMasterRank() )
        BOOST_TEST_MESSAGE( "Number of intersection: " << rirs.size() );
    BOOST_CHECK_MESSAGE( pointIntersection.size() == rirs.size(), fmt::format("Number of intersection between ray and BVH tree is not correct : {} vs {}", pointIntersection.size(), rirs.size() ) );

    int counter = 0;
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            BOOST_TEST_MESSAGE( " --  ProcessId: " << rir.processId() );
            BOOST_TEST_MESSAGE( " --  PrimitiveId: " << rir.primitiveId() );
            BOOST_TEST_MESSAGE( " --  Distance: " << rir.distance() );

            if ( rir.hasCoordinates() )
            {
                if constexpr ( std::decay_t<decltype(*bvh)>::nRealDim == 3 )
                    BOOST_TEST_MESSAGE( " --  Coordinates: " << rir.coordinates()[0] << "," << rir.coordinates()[1] << "," << rir.coordinates()[2] );
                BOOST_CHECK_SMALL( (rir.coordinates()-pointIntersection[counter]).norm(), 1e-8 );
            }

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            BOOST_TEST_MESSAGE( " --  Mesh entity id: " << prim.meshEntity().id() );
            BOOST_TEST_MESSAGE( " --  Mesh entity barycenter: " << prim.meshEntity().barycenter() );
            BOOST_TEST_MESSAGE( " ---------------------------------" );
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        counter++;
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

    std::vector<std::vector<Eigen::Vector2d>> pointIntersections;
    bvh_ray_type ray_1(origin,direction_perp_1);
    pointIntersections.push_back( { Eigen::Vector2d( 0.6, origin.y() ) } );
    bvh_ray_type ray_2(origin,direction_perp_2);
    pointIntersections.push_back( { Eigen::Vector2d( origin.x(), 0.6 ) } );

    auto bvhInHouse = boundingVolumeHierarchy(_range=range,_kind="in-house");
    auto rayIntersectionResult1 = bvhInHouse->intersect(_ray=ray_1) ;
    printRayIntersectionResults( bvhInHouse.get(),rayIntersectionResult1,pointIntersections[0] );
    auto rayIntersectionResult2 = bvhInHouse->intersect(_ray=ray_2) ;
    printRayIntersectionResults( bvhInHouse.get(),rayIntersectionResult2,pointIntersections[1] );
}

template <typename RangeType>
void test3D( RangeType const& range )
{
    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    Eigen::Vector3d origin = { 0.5,0.5,0.5 };
    Eigen::Vector3d origin2 = { -10.0,0.5,0.5 };

    Eigen::Vector3d direction_perp_1 = { 1.,0.,0. };
    Eigen::Vector3d direction_perp_2 = { 0.,1.,0. };

    std::vector<bvh_ray_type> rays;
    std::vector<std::vector<Eigen::Vector3d>> pointIntersections;
    rays.push_back( bvh_ray_type(origin,direction_perp_1) );
    pointIntersections.push_back( { Eigen::Vector3d( 1, origin.y(), origin.z() ) } );
    rays.push_back( bvh_ray_type(origin,direction_perp_2) );
    pointIntersections.push_back( { Eigen::Vector3d( origin.x(), 1., origin.z() ) } );
    rays.push_back( bvh_ray_type(origin2,direction_perp_1) );
    pointIntersections.push_back( { Eigen::Vector3d( -1, origin2.y(), origin2.z() ) } );

    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
    std::vector<std::size_t> raysDistributedIndices;
    for (int k=0;k<rays.size();++k)
    {
        raysDistributed.push_back( rays[k] );
        raysDistributedIndices.push_back( k );
    }

    auto bvhInHouse = boundingVolumeHierarchy(_range=range,_kind="in-house");
    auto bvhThirdParty = boundingVolumeHierarchy(_range=range,_kind="third-party");
    auto bvhThirdPartyLow = boundingVolumeHierarchy(_range=range,_kind="third-party",_quality=BVHEnum::Quality::Low);
    for ( auto bvhCurrent : {bvhInHouse.get(),bvhThirdParty.get(),bvhThirdPartyLow.get()} )
    {
        std::size_t counter = 0;
        for ( auto const& ray : rays )
        {
            auto rayIntersectionResult1 = bvhCurrent->intersect(_ray=ray) ;
            printRayIntersectionResults( bvhCurrent,rayIntersectionResult1, pointIntersections[counter++] );
        }
        auto multiRayDistributedIntersectionResult = bvhCurrent->intersect(_ray=raysDistributed);
        counter = 0;
        for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            printRayIntersectionResults( bvhCurrent,rayIntersectionResult, pointIntersections[raysDistributedIndices[counter++]] );
        }
    }
}

BOOST_AUTO_TEST_CASE( intersection_bvh_2D )
{
    if ( Environment::isParallel() )
        return;

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
