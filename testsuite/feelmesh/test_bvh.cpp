#define BOOST_TEST_MODULE bvh_tests
#include <feel/feelcore/testsuite.hpp>


#include <feel/feelmesh/bvh.hpp>
using namespace Feel;
using Feel::project;

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

    using bvh_type = BVHTree<FEELPP_DIM,FEELPP_DIM>;
    using bvh_ray_type = typename bvh_type::ray_type;
    BVHTree<FEELPP_DIM,FEELPP_DIM> bvh_tree;

    bvh_tree.buildPrimitivesInfo(submesh);
#if 0
    for(auto info: bvh_tree.M_primitiveInfo)
    {
        std::cout << "primitive number" << info.M_primitiveNumber <<std::endl;
        std::cout << "Min_bound" << info.M_bound_min <<std::endl;
        std::cout << "Max_bound" << info.M_bound_max <<std::endl;
        std::cout << "Centroid" << info.M_centroid <<std::endl;
    }
#endif

    bvh_tree.buildRootTree();

    // for(auto i: bvh_tree.orderedPrims )
    // {
    //     std::cout << i << std::endl;
    //     std::cout << bvh_tree.M_primitiveInfo[i].M_centroid << std::endl;
    // }
#if FEELPP_DIM==2
    Eigen::VectorXd origin(2);
    origin<< 0.5,0.5;
    Eigen::VectorXd direction_perp_1(2);
    direction_perp_1 << 1.,0.;
    Eigen::VectorXd direction_perp_2(2);
    direction_perp_2 << 0.,1.;
    bvh_ray_type ray_1(origin,direction_perp_1);
    bvh_ray_type ray_2(origin,direction_perp_2);

    int element_number = bvh_tree.raySearch(ray_1);
    BOOST_CHECK_MESSAGE(element_number > 0, fmt::format("Intersection between ray1 and BVH tree has been found"));
    element_number = bvh_tree.raySearch(ray_2);
    BOOST_CHECK_MESSAGE(element_number > 0, fmt::format("Intersection between ray2 and BVH tree has been found"));
#elif FEELPP_DIM==3
    Eigen::VectorXd origin(3);
    origin<< 0.5,0.5,0.5;
    Eigen::VectorXd direction_perp_1(3);
    direction_perp_1 << 1.,0.,0.;
    Eigen::VectorXd direction_perp_2(3);
    direction_perp_2 << 0.,1.,0.;
    bvh_ray_type ray_1(origin,direction_perp_1);
    bvh_ray_type ray_2(origin,direction_perp_2);

    int element_number = bvh_tree.raySearch(ray_1);
    BOOST_CHECK_MESSAGE(element_number > 0, fmt::format("Intersection between ray1 and BVH tree has been found"));
    element_number = bvh_tree.raySearch(ray_2);
    BOOST_CHECK_MESSAGE(element_number > 0, fmt::format("Intersection between ray2 and BVH tree has been found"));
#endif
}

BOOST_AUTO_TEST_SUITE_END()
