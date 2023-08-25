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
    return opts;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( bvh_intersection_tests )


BOOST_AUTO_TEST_CASE( intersection_bvh )
{

    auto mesh = loadMesh(_mesh = new Mesh<Simplex<FEELPP_DIM,1,FEELPP_DIM>>);
#if FEELPP_DIM==2
    auto submesh = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,{"RequiredBoundaryOfRequiredElements"}));
#elif FEELPP_DIM==3
    auto submesh = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,{"CavityBottom","CavitySides","CavityTop"}));
#endif

    auto rangeEntity = elements(submesh);
    using mesh_entity_type = std::remove_const_t<entity_range_t<std::decay_t<decltype(rangeEntity)>>>;

    //using bvh_ray_type = typename std::decay_t<decltype(*bvhInHouse)>::ray_type;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    if constexpr ( mesh_entity_type::nRealDim == 2 )
    {
        Eigen::VectorXd origin(2);
        origin<< 0.5,0.5;
        Eigen::VectorXd direction_perp_1(2);
        direction_perp_1 << 1.,0.;
        Eigen::VectorXd direction_perp_2(2);
        direction_perp_2 << 0.,1.;
        bvh_ray_type ray_1(origin,direction_perp_1);
        bvh_ray_type ray_2(origin,direction_perp_2);

        auto bvhInHouse = boundingVolumeHierarchy(_range=elements(submesh),_kind="in-house");
        auto rayIntersectionResult1 = bvhInHouse->intersect(ray_1) ;
        BOOST_CHECK_MESSAGE(!rayIntersectionResult1.empty(), fmt::format("Intersection between ray1 and BVH tree has been found"));
        auto rayIntersectionResult2 = bvhInHouse->intersect(ray_2) ;
        BOOST_CHECK_MESSAGE(!rayIntersectionResult2.empty(), fmt::format("Intersection between ray2 and BVH tree has been found"));
    }
    else if constexpr ( mesh_entity_type::nRealDim == 3 )
    {
        Eigen::VectorXd origin(3);
        origin<< 0.5,0.5,0.5;
        Eigen::VectorXd direction_perp_1(3);
        direction_perp_1 << 1.,0.,0.;
        Eigen::VectorXd direction_perp_2(3);
        direction_perp_2 << 0.,1.,0.;
        bvh_ray_type ray_1(origin,direction_perp_1);
        bvh_ray_type ray_2(origin,direction_perp_2);

        auto bvhInHouse = boundingVolumeHierarchy(_range=elements(submesh),_kind="in-house");
        auto bvhThirdParty = boundingVolumeHierarchy(_range=elements(submesh),_kind="third-party");
        for ( auto bvhCurrent : {bvhInHouse.get(),bvhThirdParty.get()} )
        {
            auto rayIntersectionResult1 = bvhCurrent->intersect(ray_1) ;
            BOOST_CHECK_MESSAGE(!rayIntersectionResult1.empty(), fmt::format("Intersection between ray1 and BVH tree has been found"));
            auto rayIntersectionResult2 = bvhCurrent->intersect(ray_2) ;
            BOOST_CHECK_MESSAGE(!rayIntersectionResult2.empty(), fmt::format("Intersection between ray2 and BVH tree has been found"));
        }
    }

}

BOOST_AUTO_TEST_SUITE_END()
