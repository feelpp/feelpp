#define BOOST_TEST_MODULE bvhgpu_tests


#define COMPILE_WITH_HIP

#include <feel/feelmesh/ranges.hpp>

#include <feel/feelcore/enumerate.hpp>

#include <feel/feelcore/environment.hpp>

#include <feel/feelcore/kokkos.hpp>

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

#include <feel/feelfilters/unitcube.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feells/distancetorange.hpp>


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



inline
AboutData
makeAbout()
{
    AboutData about("test_distance BVH RT CPU AND GPU",
        "test_distance",
        "0.1",
        "nD(n=3)",
        Feel::AboutData::License_GPL,
        "Copyright (c) 2024 Feel++ Consortium");

    about.addAuthor("Noname", "developer", "Noname@cemosis.fr", "");
    return about;
}

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts("Test Environment options");
    opts.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "number_rays_desired", po::value<int>()->default_value( 703 ), "mesh size" )
        ;
    return opts;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS(makeAbout(), makeOptions());

struct DataDistanceErrTime {
    int id;
    double distanceMinREAL;
    double distanceFastMarching;
    double errFastMarching;
    double distanceMinCPU;
    double errCPU;
    double distanceMinGPU;
    double errGPU;
    long int t_laps_CPU;
    long int t_laps_GPU;
};

struct DataTimeLapsConfig {
    int  nbRays; // <= uniform distribution
    int  nbRaysDesired;
    double hsize;
    long int t_laps_BVH_CPU;
    long int t_laps_BVH_GPU;
    long int t_laps_RT_CPU;
    long int t_laps_RT_GPU;
    long int t_laps_FastMarching;
};




template <typename BvhType, typename RayIntersectionResultType>
std::vector<double> getAllDistanceRayIntersections(BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs)
{
    std::vector<double> distance;
    for (auto const& rir : rirs)
    {
        if (rir.processId() == bvh->worldComm().rank())
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


double calculateStepRays(int n, double start, double end) {
    double n_total = std::sqrt(n);
    return (end - start) / (n_total - 1);
}


template <typename RangeType2>
void distToBoundaryBVHpu(
    RangeType2 const& range,
    DataTimeLapsConfig& allDataPU,
    std::vector<std::vector<double>>& allNodeCoordinates,
    std::vector<DataDistanceErrTime>& allDataDistanceBVHRT,
    bool isViewInfo)
{
    std::chrono::steady_clock::time_point t_begin_cpu, t_begin_gpu;
    std::chrono::steady_clock::time_point t_end_cpu, t_end_gpu;

    std::chrono::steady_clock::time_point t_begin_raytracing_cpu, t_begin_raytracing_gpu;
    std::chrono::steady_clock::time_point t_end_raytracing_cpu, t_end_raytracing_gpu;
    std::chrono::steady_clock::time_point t_end_bvh_cpu, t_end_bvh_gpu;

    long int t_laps_CPU = 0;
    long int t_laps_GPU = 0;

    long int t_laps_CPU_Total = 0;
    long int t_laps_GPU_Total = 0;

    int nbRays = 0;

    int number_rays_desired = allDataPU.nbRaysDesired;

    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType2>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    // BEGIN::Build BVH CPU
    t_begin_cpu = std::chrono::steady_clock::now();
    auto bvhThirdPartyLow = boundingVolumeHierarchy(_range = range, _kind = "third-party");
    t_end_bvh_cpu = std::chrono::steady_clock::now();
    // END::Build BVH CPU

    sleep(1);

    // BEGIN::Build BVH GPU
    t_begin_gpu = std::chrono::steady_clock::now();
    auto bvhHIPParty = boundingVolumeHierarchy(_range = range, _kind = "hip-party");
    t_end_bvh_gpu = std::chrono::steady_clock::now();
    // END::Build BVH GPU

    sleep(1);


    t_laps_CPU = std::chrono::duration_cast<std::chrono::microseconds>(t_end_bvh_cpu - t_begin_cpu).count();
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside BVH CPU : " << t_laps_CPU << " us\n";

    t_laps_GPU = std::chrono::duration_cast<std::chrono::microseconds>(t_end_bvh_gpu - t_begin_gpu).count();
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside BVH GPU : " << t_laps_GPU << " us\n";

    const double epsilon = 0.00001f;
    double distanceMinCPU = 0.0f;
    double distanceMinGPU = 0.0f;
    double distanceMinREAL = 0.0f;
    int nbValues = 0;

    for (int k = 0; k < allNodeCoordinates.size(); ++k)
    {
        // Build ray
        std::vector<bvh_ray_type> rays;
        std::vector<double> distance_Real_mode;
        BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
        Eigen::Vector3d ray_origin = { allNodeCoordinates[k][0], allNodeCoordinates[k][1], allNodeCoordinates[k][2] };

        bool ok = false;
        //if ( (ray_origin[0]>0.0f) && (ray_origin[0]<1.0f) && (ray_origin[1]>0.0f) && (ray_origin[1]<1.0f) && (ray_origin[2]>0.0f) && (ray_origin[2]<1.0f)) { ok = true; }

        ok = true; // All points

        if (ok) {
            nbValues++;
            distanceMinREAL = std::min(ray_origin[0], std::min(1.0 - ray_origin[0], std::min(ray_origin[1], std::min(1.0 - ray_origin[1], std::min(ray_origin[2], 1.0 - ray_origin[2])))));

            double thetaStart = 0.0f;
            double thetaEnd = M_PI;
            double alphaStart = 0.0f;
            double alphaEnd = 2.0f * M_PI;
            double thetaStep = calculateStepRays(number_rays_desired, thetaStart, thetaEnd);
            double alphaStep = calculateStepRays(number_rays_desired, alphaStart, alphaEnd);

            for (double theta = thetaStart; theta <= thetaEnd; theta += thetaStep) {
                for (double alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep) {
                    Eigen::Vector3d ray_direction = sphericalToCartesian(1.0f, theta, alpha);
                    Eigen::Vector3d ray_origin_Epsilon = ray_origin;
                    rays.push_back(bvh_ray_type(ray_origin_Epsilon, ray_direction));
                }
            }

            nbRays = rays.size();

            for (int i = 0; i < rays.size(); ++i)
            {
                raysDistributed.push_back(rays[i]);
            }

            std::vector<double> dist;
            std::vector<double> distance_CPU_mode;
            std::vector<double> distance_GPU_mode;

            // In normal CPU mode
            // Ray Tracing
            t_begin_raytracing_cpu = std::chrono::steady_clock::now();
            auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect(_ray = raysDistributed);
            for (auto const& rayIntersectionResult : multiRayDistributedIntersectionResult)
            {
                dist = getAllDistanceRayIntersections(bvhThirdPartyLow, rayIntersectionResult);
                distance_CPU_mode.insert(distance_CPU_mode.end(), dist.begin(), dist.end());
            }
            t_end_raytracing_cpu = std::chrono::steady_clock::now();


            t_begin_raytracing_gpu = std::chrono::steady_clock::now();
            auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect(_ray = raysDistributed);
            for (auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult)
            {
                dist = getAllDistanceRayIntersections(bvhHIPParty, rayIntersectionResult);
                distance_GPU_mode.insert(distance_GPU_mode.end(), dist.begin(), dist.end());
            }
            t_end_raytracing_gpu = std::chrono::steady_clock::now();


            distanceMinCPU = INFINITY;
            for (int k = 0; k < distance_CPU_mode.size(); ++k)
            {
                distanceMinCPU = fmin(distanceMinCPU, distance_CPU_mode[k]);
            }

            distanceMinGPU = INFINITY;

            for (int k = 0; k < distance_GPU_mode.size(); ++k)
            {
                distanceMinGPU = fmin(distanceMinGPU, distance_GPU_mode[k]);
            }


            t_laps_CPU = std::chrono::duration_cast<std::chrono::microseconds>(t_end_raytracing_cpu - t_begin_raytracing_cpu).count();
            t_laps_GPU = std::chrono::duration_cast<std::chrono::microseconds>(t_end_raytracing_gpu - t_begin_raytracing_gpu).count();

            t_laps_CPU_Total = t_laps_CPU_Total + t_laps_CPU;
            t_laps_GPU_Total = t_laps_GPU_Total + t_laps_GPU;

            double errCPU = abs(distanceMinCPU - distanceMinREAL);
            double errGPU = abs(distanceMinGPU - distanceMinREAL);

            //double errFastMarching = abs(distanceFastMarching[k] - distanceMinREAL);

            if (isViewInfo) std::cout << "[INFO] [" << k
                << "]"
                << "<" << std::setprecision(5)<< ray_origin[0]
                << "," << std::setprecision(5)<< ray_origin[1]
                << "," << std::setprecision(5)<< ray_origin[2]
                << ">"
                << " Distance Min REAL=" << distanceMinREAL
                //<< " FastMarching=" << distanceFastMarching[k]
                //<< " err=" << errFastMarching
                << " CPU=" << std::setprecision(5)<< distanceMinCPU
                << " err=" << std::setprecision(5)<< errCPU
                << " GPU=" << std::setprecision(5)<< distanceMinGPU
                << " err=" << std::setprecision(5)<< errGPU
                << " t_laps_CPU=" << t_laps_CPU
                << " t_laps_GPU=" << t_laps_GPU
                << "\n";

            DataDistanceErrTime data = {
                k,
                distanceMinREAL,
                -1,
                -1,
                distanceMinCPU,
                errCPU,
                distanceMinGPU,
                errGPU,
                t_laps_CPU,
                t_laps_GPU
            };

            allDataDistanceBVHRT.push_back(data);

            // Memory cleaning
            distance_CPU_mode.clear();
            distance_GPU_mode.clear();
        }
    } // END for k

    // Elapse Time BVH - RT - CPU - GPU
    allDataPU.t_laps_BVH_CPU = std::chrono::duration_cast<std::chrono::microseconds>(t_end_bvh_cpu - t_begin_cpu).count();
    allDataPU.t_laps_RT_CPU = t_laps_CPU_Total / nbValues;
    allDataPU.t_laps_BVH_GPU = std::chrono::duration_cast<std::chrono::microseconds>(t_end_bvh_gpu - t_begin_gpu).count();
    allDataPU.t_laps_RT_GPU = t_laps_GPU_Total / nbValues;
    allDataPU.nbRays = nbRays;
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside BVH CPU : " << allDataPU.t_laps_BVH_CPU << " us\n";
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside Ray Tracing CPU : " << allDataPU.t_laps_RT_CPU << " us\n";
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside BVH GPU : " << allDataPU.t_laps_BVH_GPU << " us \n";
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside Ray Tracing GPU : " << allDataPU.t_laps_RT_GPU << " us\n";
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside FastMarching : " << allDataPU.t_laps_FastMarching << " us\n";
}


BOOST_AUTO_TEST_SUITE(distance_bvh_cpu_gpu_gpu_tests)


BOOST_AUTO_TEST_CASE(all_distance)
{

    bool isViewInfo = true;  //isViewInfo = false;;
    double hsize = 1.0f / 2.0f;
    int number_rays_desired = 703;


    using namespace Feel;
    using Feel::cout;
    using mesh_type = Mesh<Simplex<3, 1, 3>>; //<Dim,Order,RDim>
    auto mesh = unitCube();
    if (isViewInfo)
    {
        std::cout << "[INFO] maxNumElement : " << mesh->maxNumElements() << std::endl;
        std::cout << "[INFO] maxNumFace    : " << mesh->maxNumFaces() << std::endl;
        std::cout << "[INFO] maxNumPoints  : " << mesh->maxNumPoints() << std::endl;
        std::cout << "[INFO] maxNumVerices : " << mesh->maxNumVertices() << std::endl;
    }

    auto rangeFaces = markedfaces(mesh);
    auto submeshFaces  = boundaryfaces( mesh );
    auto rangeElements = markedelements(mesh);


    auto Vh = Pch<1>(mesh);
#if 0
    auto const& nodes = Vh->mesh()->points();

    for (auto const& pointPair :nodes)
    {
        auto const& point = pointPair.second;
        auto const& coords = point.node();
        allNodeCoordinates.push_back({ coords[0], coords[1], coords[2] });
    }
#endif

    std::vector<std::vector<double>> allNodeCoordinates;
    DataTimeLapsConfig allDataPU;
    allDataPU.nbRaysDesired=number_rays_desired;
    allDataPU.hsize=hsize;
    std::vector<DataDistanceErrTime> allDataDistanceBVHRT;
    
    // List of node coordinates
    for (size_type k=0;k<Vh->nLocalDofWithGhost();++k)
    {
        auto const& dofPt = Vh->dof()->dofPoint(k).template get<0>();
        allNodeCoordinates.push_back({ dofPt[0], dofPt[1], dofPt[2] });
    }

    int nbNode=allNodeCoordinates.size();


    // Calculates Node points to Surface distances by the method FastMarching
    std::chrono::steady_clock::time_point t_begin_FastMarching, t_end_FastMarching;
    t_begin_FastMarching = std::chrono::steady_clock::now();
        auto distToBoundary = distanceToRange( _space=Vh, _range=submeshFaces);
    t_end_FastMarching = std::chrono::steady_clock::now();
    long int t_laps_FastMarching = std::chrono::duration_cast<std::chrono::microseconds>(t_end_FastMarching - t_begin_FastMarching).count();

    allDataPU.t_laps_FastMarching=t_laps_FastMarching;


    // Calculates Node points to Surface distances by the method BVH RT CPU and GPU
    distToBoundaryBVHpu(rangeFaces,allDataPU,allNodeCoordinates,allDataDistanceBVHRT,isViewInfo);

    //We fill the Fast Marching distance data into the allDataDistanceBVHRT data structure
    for (int i = 0; i < nbNode; ++i)
    {
        allDataDistanceBVHRT[i].distanceFastMarching = distToBoundary[i];
        allDataDistanceBVHRT[i].errFastMarching = abs(distToBoundary[i]-allDataDistanceBVHRT[i].distanceMinREAL);
    }

    // Debriefing Save all data
    std::string filenameDataDistanceErrTime = "all_results_per_vertex.csv";
    if (remove(filenameDataDistanceErrTime.c_str()) != 0) {
        std::cerr << "Error delete file." << std::endl;
    }
    std::ofstream myfileB(filenameDataDistanceErrTime);
    myfileB << "Num Vertex,PosX,PosY,PosZ,distanceMinREAL,distanceFastMarching,errFastMarching,distanceMinCPU,errCPU,distanceMinGPU,errGPU,time_RT_CPU,time_RT_GPU\n";
    for (int i = 0; i < nbNode; ++i)
    {
        myfileB << allDataDistanceBVHRT[i].id << ","
                << std::setprecision(5)<< allNodeCoordinates[i][0] << ","
                << std::setprecision(5)<< allNodeCoordinates[i][1] << ","
                << std::setprecision(5)<< allNodeCoordinates[i][2] << ","
                << std::setprecision(5)<< allDataDistanceBVHRT[i].distanceMinREAL << ","
                << std::setprecision(5)<< allDataDistanceBVHRT[i].distanceFastMarching << ","
                << std::setprecision(5)<< allDataDistanceBVHRT[i].errFastMarching << ","
                << std::setprecision(5)<< allDataDistanceBVHRT[i].distanceMinCPU << ","
                << std::setprecision(5)<< allDataDistanceBVHRT[i].errCPU << ","
                << std::setprecision(5)<< allDataDistanceBVHRT[i].distanceMinGPU << ","
                << std::setprecision(5)<< allDataDistanceBVHRT[i].errGPU << ","
                << allDataDistanceBVHRT[i].t_laps_CPU << ","
                << allDataDistanceBVHRT[i].t_laps_GPU << "\n";
    }
    myfileB.close();

    std::string filenameA = "results.csv";
    if (remove(filenameA.c_str()) != 0) {
        std::cerr << "Error delete file." << std::endl;
    }
    std::ofstream myfileA(filenameA);
    myfileA << "hsize=" << allDataPU.hsize<< "\n";
    myfileA << "maxNumElement= " << mesh->maxNumElements() << "\n";
    myfileA << "maxNumFace=" << mesh->maxNumFaces() << "\n";
    myfileA << "maxNumPoints=" << mesh->maxNumPoints() << "\n";
    myfileA << "maxNumVerices=" << mesh->maxNumVertices() << "\n";
    myfileA << "nbRaysDesired=" << allDataPU.nbRaysDesired<< "\n";
    myfileA << "nbRays=" << allDataPU.nbRays<< "\n";
    myfileA << "timeBVHcpu=" << allDataPU.t_laps_BVH_CPU<< "\n";
    myfileA << "timeMeanRTcpu=" << allDataPU.t_laps_RT_CPU << "\n";
    myfileA << "timeBVHgpu=" << allDataPU.t_laps_BVH_GPU << "\n";
    myfileA << "timeMeanRTgpu=" << allDataPU.t_laps_RT_GPU << "\n";
    myfileA << "timeFastMarching=" << allDataPU.t_laps_FastMarching<< "\n";
    myfileA << "totalTimeBVHRTcpu=" << allDataPU.t_laps_BVH_CPU+allDataPU.t_laps_RT_CPU<< "\n";
    myfileA << "totalTimeBVHRTgpu=" << allDataPU.t_laps_BVH_GPU+allDataPU.t_laps_RT_GPU<< "\n";
    myfileA.close();

    // Data backup file distances for paraview
    auto exp = exporter( _mesh = mesh, _name = fmt::format( "distance_{}d_o{}", 3, 1 ) );
    exp->addRegions();
    exp->add( "distToBoundary", distToBoundary );
    exp->save();


    // Memory cleaning
    allNodeCoordinates.clear();
    allDataDistanceBVHRT.clear();

    
}



BOOST_AUTO_TEST_SUITE_END()




