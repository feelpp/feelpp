#define BOOST_TEST_MODULE bvhgpu_tests

// NOTA : Objective: Ray tracing from inside a cube using BVH Ray Tracing with a CPU and a GPU method. Compare the performances of the two methods.


#define COMPILE_WITH_HIP

#include <feel/feelmesh/ranges.hpp>

#include <feel/feelcore/enumerate.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/bvh.hpp>
#include <fmt/chrono.h>

#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelvf/vf.hpp>


#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"

#include "thrust/device_vector.h"
#include "thrust/transform.h"
#include "thrust/functional.h"
#include "thrust/execution_policy.h"
#include "thrust/random.h"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "thrust/sort.h"

#include "thrust/generate.h"
#include "thrust/sort.h"
#include "thrust/copy.h"
#include "thrust/count.h"

#include <hwloc.h>



using namespace Feel;


// Section to use if you want to make meshes, otherwise it doesn't work.
inline
AboutData
makeAbout()
{
    AboutData about( "test_bvhgpu" ,
                     "test_bvhgpu" ,
                     "0.2",
                     "nD(n=2,3) test bvh",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2024 Feel++ Consortium" );

    about.addAuthor( "Noname", "developer", "Noname@cemosis.fr", "" );
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



template <typename BvhType,typename RayIntersectionResultType>
void printRayIntersectionResults2( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs,std::vector<typename std::remove_pointer_t<BvhType>::vector_realdim_type  /*Eigen::Vector3d*/> const& pointIntersection )
{
    if ( bvh->worldComm().isMasterRank() )
        std::cout << "Number of intersection: " << rirs.size()<< "\n";

    int counter = 0;
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            std::cout << " --  ProcessId: " << rir.processId()<< "\n";
            std::cout << " --  PrimitiveId: " << rir.primitiveId()<< "\n";
            std::cout << " --  Distance: " << rir.distance()<< "\n";

            if ( rir.hasCoordinates() )
            {
                if constexpr ( std::decay_t<decltype(*bvh)>::nRealDim == 3 )
                    std::cout <<  " --  Coordinates: " << rir.coordinates()[0] << "," << rir.coordinates()[1] << "," << rir.coordinates()[2]<< "\n";;
            }

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            std::cout << " --  Mesh entity id: " << prim.meshEntity().id()<< "\n";
            std::cout << " --  Mesh entity barycenter: " << prim.meshEntity().barycenter()<< "\n";
            std::cout <<" ---------------------------------"<< "\n" ;
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        counter++;
    }
}

template <typename BvhType,typename RayIntersectionResultType>
void printRayIntersectionResults( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    if ( bvh->worldComm().isMasterRank() )
        std::cout << "Number of intersection: " << rirs.size()<< "\n";
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            std::cout << " --  ProcessId: " << rir.processId()<< "\n";
            std::cout << " --  PrimitiveId: " << rir.primitiveId()<< "\n";
            std::cout << " --  Distance: " << rir.distance()<< "\n";

            if ( rir.hasCoordinates() )
            {
                if constexpr ( std::decay_t<decltype(*bvh)>::nRealDim == 3 )
                    std::cout <<  " --  Coordinates: " << rir.coordinates()[0] << "," << rir.coordinates()[1] << "," << rir.coordinates()[2]<< "\n";;
            }

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            std::cout << " --  Mesh entity id: " << prim.meshEntity().id()<< "\n";
            std::cout << " --  Mesh entity barycenter: " << prim.meshEntity().barycenter()<< "\n";
            //std::cout << " --  Mesh entity barycenter: " << prim.meshEntity()
            std::cout <<" ---------------------------------"<< "\n" ;
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}
 
// Function that detects if there is GPU in the environment in order to switch the calculations to the GPU.
bool isThereAnyGPUhere() {
    hwloc_topology_t topology;
    hwloc_obj_t obj = nullptr;   
    bool isGPU      = false;
    unsigned n, i;  

    if (hwloc_topology_init(&topology) < 0) {
        std::cerr << "Error Num" << std::endl;
        return false;
    }

    hwloc_topology_set_io_types_filter(topology, HWLOC_TYPE_FILTER_KEEP_IMPORTANT);

    if (hwloc_topology_load(topology) < 0) {
        std::cerr << "Error Num" << std::endl;
        hwloc_topology_destroy(topology);
        return false;
    }

    n = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_OS_DEVICE);
    std::cout<<"n="<<n<<"\n";

    for (i = 0; i < n ; i++) {
        obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_OS_DEVICE, i);
        printf("[INFO]: %s:\n", obj->name);
        const char *s;
        s = hwloc_obj_get_info_by_name(obj, "Backend");
        printf("[INFO]: %s\n",s);
        if (s && !strcmp(s, "OpenCL")) { isGPU=true; };
    }

    hwloc_topology_destroy(topology);
    return isGPU;
}

 
template <typename RangeType>
void test3DWithHybrid( RangeType const& range )
{
    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    Eigen::Vector3d origin1={1000.0,0.0,0.0};
    Eigen::Vector3d direction_perp_1={1.,0.,0.};

    Eigen::Vector3d origin2={-10.0,-0.25,-0.25};
    Eigen::Vector3d direction_perp_2={1.,0.,0.};

    Eigen::Vector3d origin3={-5.0,-0.20,-0.20};
    Eigen::Vector3d direction_perp_3={1.,0.,0.};


    std::vector<bvh_ray_type> rays;

    rays.push_back( bvh_ray_type(origin1,direction_perp_1) );
    //rays.push_back( bvh_ray_type(origin2,direction_perp_2) );
    //rays.push_back( bvh_ray_type(origin3,direction_perp_3) );

    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;

    for (int k=0;k<rays.size();++k)
    {
        raysDistributed.push_back( rays[k] );
    }
    

    std::cout<<"\n";

    // In normal CPU mode
    std::cout<<"[INFO]: In normal CPU mode\n";
    auto bvhThirdPartyLow = boundingVolumeHierarchy(_range=range,_kind="third-party");
    auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            std::cout << "multiRayDistributedIntersectionResult parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhThirdPartyLow,rayIntersectionResult);
        }

    // In GPU mode with AMD HIP
    std::cout<<"[INFO]: In GPU mode with AMD HIP\n";
    auto bvhHIPParty       = boundingVolumeHierarchy(_range=range,_kind="hip-party");   
    auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
        {
            std::cout << "Hip parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhHIPParty,rayIntersectionResult);
        }
}



template <typename BvhType,typename RayIntersectionResultType>
std::vector<double> getAllDistanceRayIntersections( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    std::vector<double> distance;
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            //std::cout << " --  Distance: " << rir.distance()<< "\n";     
            distance.push_back(rir.distance());     
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    return distance;
}



Eigen::Vector3d sphericalToCartesian(double r, double theta, double alpha) {
    Eigen::Vector3d position;
    position[0] = r * sin(theta) * cos(alpha);
    position[1] = r * sin(theta) * sin(alpha);
    position[2] = r * cos(theta);
    return position;
}

template <typename RangeType>
void test3DInsideObjectWithHybrid( RangeType const& range )
{
    std::chrono::steady_clock::time_point t_begin_cpu,t_begin_gpu;
    std::chrono::steady_clock::time_point t_end_cpu,t_end_gpu;

    std::chrono::steady_clock::time_point t_begin_raytracing_cpu,t_begin_raytracing_gpu;
    std::chrono::steady_clock::time_point t_end_raytracing_cpu,t_end_raytracing_gpu;
    std::chrono::steady_clock::time_point t_end_bvh_cpu,t_end_bvh_gpu;

    long int t_laps;

    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;


    // Save all informations
    std::string filename= "results.txt";
    std::ofstream myfile (filename);

    //for (int kkk=1; kkk<=20; kkk++) 
    for (int kkk=1; kkk<=1; kkk++) 
    {

        Eigen::Vector3d ray_origin={0.0f,0.0f,0.0f};
        double thetaStart = 0.0f;        
        double thetaEnd   = M_PI;       
        double alphaStart = 0.0f;        
        double alphaEnd   = 2.0f * M_PI;   
        double thetaStep  = M_PI / (18.0f*0.125f*float(kkk)); 
        double alphaStep  = M_PI / (18.0f*0.125f*float(kkk)); 
        bool   isViewInfo = false;

        std::vector<bvh_ray_type> rays;
        for (double theta = thetaStart; theta <= thetaEnd; theta += thetaStep) {
            for (double alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep) {
                Eigen::Vector3d ray_direction = sphericalToCartesian(1.0f, theta, alpha);

                if (isViewInfo)
                {
                    std::cout << "Origin: <" << ray_origin[0] << ", " << ray_origin[1] << ", " << ray_origin[2] << ">" << " ";
                    std::cout << "Direction: <" << ray_direction[0] << ", " << ray_direction[1] << ", " << ray_direction[2] << ">" << std::endl;              
                }
                rays.push_back( bvh_ray_type(ray_origin,ray_direction) );
            }
        }

        BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
        for (int k=0;k<rays.size();++k)
        {
            raysDistributed.push_back( rays[k] );
        }

        std::vector<double> dist;



        // In normal CPU mode
        std::cout<<"[INFO]: In normal CPU mode\n";

        t_begin_cpu                                = std::chrono::steady_clock::now();
        auto bvhThirdPartyLow                      = boundingVolumeHierarchy(_range=range,_kind="third-party");
        t_end_bvh_cpu                              = std::chrono::steady_clock::now();

        t_begin_raytracing_cpu                     = std::chrono::steady_clock::now();
        auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect(_ray=raysDistributed);
        t_end_raytracing_cpu                       = std::chrono::steady_clock::now();


        std::vector<double> distance_CPU_mode;
        for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            dist=getAllDistanceRayIntersections(bvhThirdPartyLow,rayIntersectionResult);
            //printRayIntersectionResults(bvhThirdPartyLow,rayIntersectionResult);
            distance_CPU_mode.insert(distance_CPU_mode.end(), dist.begin(), dist.end());
        }
        t_end_cpu = std::chrono::steady_clock::now();



        // In GPU mode with AMD HIP
        std::cout<<"[INFO]: In GPU mode with AMD HIP\n";

        t_begin_gpu                                   = std::chrono::steady_clock::now();
        auto bvhHIPParty                              = boundingVolumeHierarchy(_range=range,_kind="hip-party");   
        t_end_bvh_gpu                                 = std::chrono::steady_clock::now();

        t_begin_raytracing_gpu                        = std::chrono::steady_clock::now();
        auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect(_ray=raysDistributed);
        t_end_raytracing_gpu                          = std::chrono::steady_clock::now();

        std::vector<double> distance_GPU_mode;
        for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
        {
            dist=getAllDistanceRayIntersections(bvhHIPParty,rayIntersectionResult);
            //printRayIntersectionResults(bvhHIPParty,rayIntersectionResult);
            distance_GPU_mode.insert(distance_GPU_mode.end(), dist.begin(), dist.end());
        }
        t_end_gpu = std::chrono::steady_clock::now();

        // Distance comparison
        double deltaError = 0.00001f;
        double sumErrors  = 0.0f;
        bool   isError = false;



        std::cout<<"[INFO]: Size vector distance CPU="<<distance_CPU_mode.size()<<" GPU="<<distance_GPU_mode.size()<<"\n";

        if (distance_GPU_mode.size()!=distance_CPU_mode.size())
        {
            isError = true;
            std::cout << "[INFO]: Error size vector distance CPU vs GPU\n"; 
        }

        for ( int k; k<distance_GPU_mode.size(); ++k)
        {
            double e=fabs(distance_GPU_mode[k]-distance_CPU_mode[k]);
            sumErrors = sumErrors + e;
            if (e>deltaError) 
            { 
                std::cout << "[INFO]: ERROR "<<k<<" Distance CPU = "<<distance_CPU_mode[k]<<" Dist GPU = "<<distance_GPU_mode[k]<<" e="<<e<<"\n"; 
                isError = true;
            }
        }
        if (!isError) {  std::cout << "[INFO]: WELL DONE :-) No error (same distance). \n";  } 


        // Save all informations
        //std::string filename= "results.txt";
        //std::ofstream myfile (filename);

        if (myfile.is_open())
        {
            std::cout << "[INFO]: Nb Rays : "<<rays.size()<<"\n";
            //myfile <<  "Nb Rays,"<<rays.size()<< "\n";

            // Time laps comparison
            std::cout << "[INFO]: Build BVH\n";

            if (1==0) {
                t_laps= std::chrono::duration_cast<std::chrono::microseconds>(t_end_bvh_cpu - t_begin_cpu).count();
                std::cout << "[INFO]: Elapsed microseconds inside BVH CPU : "<<t_laps<< " us\n";
                //myfile <<  "Elapsed microseconds inside BVH CPU,"<<t_laps<< "\n";

                t_laps= std::chrono::duration_cast<std::chrono::microseconds>(t_end_bvh_gpu - t_begin_gpu).count();
                std::cout << "[INFO]: Elapsed microseconds inside BVH GPU : "<<t_laps<< " us (+load mesh in GPU)\n";
                //myfile <<  "Elapsed microseconds inside BVH GPU, "<<t_laps<< "\n";

                std::cout << "[INFO]: Ray Tracing\n";
                t_laps= std::chrono::duration_cast<std::chrono::microseconds>(t_end_raytracing_cpu  - t_begin_raytracing_cpu).count();
                std::cout << "[INFO]: Elapsed microseconds inside Ray Tracing CPU : "<<t_laps<< " us\n";
                //myfile << "Elapsed microseconds inside Ray Tracing CPU," <<t_laps<< "\n";

                t_laps= std::chrono::duration_cast<std::chrono::microseconds>(t_end_raytracing_gpu  - t_begin_raytracing_gpu).count();
                std::cout << "[INFO]: Elapsed microseconds inside Ray Tracing GPU : "<<t_laps<< " us (+load rays data in GPU)\n";
                //myfile <<  "Elapsed microseconds inside Ray Tracing GPU,"<<t_laps<< "\n";
            }
 

            std::cout << "[INFO]: Elapse all\n";

            t_laps= std::chrono::duration_cast<std::chrono::microseconds>(t_end_cpu - t_begin_cpu).count();
            std::cout << "[INFO]: Elapsed microseconds inside BVH Ray Tracing CPU : "<<t_laps<< " us\n";
            //myfile <<  "Elapsed microseconds inside  BVH Ray Tracing CPU,"<<t_laps<< "\n";

            t_laps= std::chrono::duration_cast<std::chrono::microseconds>(t_end_gpu - t_begin_gpu).count();
            std::cout << "[INFO]: Elapsed microseconds inside BVH Ray Tracing GPU : "<<t_laps<< " us\n";
            //myfile <<  "Elapsed microseconds inside  BVH Ray Tracing GPU,"<<t_laps<< "\n";

            myfile <<  "Nb Rays,"<<rays.size()<< ",";
            t_laps= std::chrono::duration_cast<std::chrono::microseconds>(t_end_cpu - t_begin_cpu).count();
            myfile <<t_laps<< ",";
            t_laps= std::chrono::duration_cast<std::chrono::microseconds>(t_end_gpu - t_begin_gpu).count();
            myfile <<t_laps<< "\n";


        }
        else std::cout << "Unable to open file";
        
    }

    myfile.close();

}




BOOST_AUTO_TEST_SUITE( bvh_intersection_gpu_tests )


BOOST_AUTO_TEST_CASE( test_load_mesh3 )
{

    //TEST Thrust
    thrust::device_vector<float> d_x(10, 1.0f);
    using namespace Feel;
    using Feel::cout;   
    //typedef Mesh<Simplex<2> > mesh_type;

    using mesh_type = Mesh<Simplex<3,1,3>>; //<Dim,Order,RDim>

    auto mesh = loadMesh(_mesh=new  mesh_type,_filename="/nvme0/lemoinep/feelppGPUBeta/feelpp/testsuite/feelbvhgpu/cases/mesh/Test.geo");

    std::cout << "[INFO] maxNumElement : "<< mesh->maxNumElements() << std::endl;
    std::cout << "[INFO] maxNumFace    : "<< mesh->maxNumFaces() << std::endl;
    std::cout << "[INFO] maxNumPoints  : "<< mesh->maxNumPoints() << std::endl;
    std::cout << "[INFO] maxNumVerices : "<< mesh->maxNumVertices() << std::endl;

    auto rangeFaces = markedfaces(mesh,{"CavityBottom","CavitySides","CavityTop","Up3","Down3","Back3","Left3","Rigth3"});
    auto submesh = createSubmesh(_mesh=mesh,_range=rangeFaces);

    auto Xhd0 = Pdh<0>(mesh);
    auto measures = Xhd0->element();
    measures.on(_range=elements(mesh),_expr=vf::meas());
    double measMin = measures.min();
    double measMax = measures.max();
    size_type nbdyfaces = nelements(boundaryfaces(mesh));

        std::cout << "[INFO]: mesh entities" << std::endl;
        std::cout << "[INFO]:    number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << "[INFO]:    number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << "[INFO]:    number of boundary faces : " << nbdyfaces << std::endl;
        std::cout << "[INFO]:    number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << "[INFO]:    number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << "[INFO]: mesh sizes" << std::endl;
        std::cout << "[INFO]:    h max : " << mesh->hMax() << std::endl;
        std::cout << "[INFO]:    h min : " << mesh->hMin() << std::endl;
        std::cout << "[INFO]:    h avg : " << mesh->hAverage() << std::endl;
        std::cout << "[INFO]:    measure : " << mesh->measure() << "\t" << measMin<<" : " << measMax << std::endl;
        std::cout << "[INFO]: Number of Partitions : " << mesh->numberOfPartitions() << std::endl ;
        std::cout << "[INFO]: nbdyfaces : " <<nbdyfaces<<"\n";
        std::cout<<"\n";

    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    // GPU or not GPU
    bool isGPUdetected; //put elsewhere
    isGPUdetected=isThereAnyGPUhere();

    if (isGPUdetected)
    {
        std::cout << "[INFO]: GPU detected" << std::endl;
    }
    else
    {
        std::cout << "[INFO]: No GPU detected" << std::endl;
    }
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    test3DWithHybrid( rangeFaces );
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    test3DInsideObjectWithHybrid( rangeFaces );
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    std::cout<<"\n";
}


BOOST_AUTO_TEST_SUITE_END()




