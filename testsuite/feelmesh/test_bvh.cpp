// #define BOOST_TEST_MODULE bvh tests
// #include <feel/feelcore/testsuite.hpp>


#include <feel/feelmesh/bvh.hpp>
using namespace Feel;

int main(int argc , char **argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="bvh",                                 
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    auto mesh = loadMesh(_mesh = new Mesh<Simplex<FEELPP_DIM,1>>);
#if FEELPP_DIM==2    
    auto submesh = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,{"RequiredBoundaryOfRequiredElements"}));
#elif FEELPP_DIM==3
    auto submesh = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,{"CavityBottom","CavitySides","CavityTop"}));
#endif

    BVH_tree<FEELPP_DIM> bvh_tree;

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
    Ray_bvh ray_1(origin,direction_perp_1);
    Ray_bvh ray_2(origin,direction_perp_2);

    bvh_tree.ray_search(ray_1,"");
#elif FEELPP_DIM==3  
    Eigen::VectorXd origin(3);
    origin<< 0.5,0.5,0.5;   
    Eigen::VectorXd direction_perp_1(3); 
    direction_perp_1 << 1.,0.,0.;  
    Eigen::VectorXd direction_perp_2(3);
    direction_perp_2 << 0.,1.,0.;      
    Ray_bvh ray_1(origin,direction_perp_1);
    Ray_bvh ray_2(origin,direction_perp_2);

    bvh_tree.ray_search(ray_1,"");
    std::cout <<"============================"<< std::endl;
    bvh_tree.ray_search(ray_2,"");
#endif
    return 0;
}
