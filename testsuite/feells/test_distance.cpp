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
    AboutData about("test_distance",
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
        //( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        //( "number_rays_desired", po::value<int>()->default_value( 703 ), "mesh size" )
        ("mesh2D.filename", po::value<std::string>(), "mesh2D.filename")
        ("mesh3D.filename", po::value<std::string>(), "mesh3D.filename")
        ;
    return opts;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS(makeAbout(), makeOptions());

// BEGIN::Global data
std::vector<std::vector<double>> allNodeCoordinates;
std::vector<double> distanceFastMarching;
long int t_laps_FastMarching;
// END::Global data


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

template <typename MeshEntityType>
struct MeshPrimitiveInfo
{
    using mesh_entity_type = std::decay_t<MeshEntityType>;
    static constexpr uint16_type nDim = mesh_entity_type::nDim;
    static constexpr uint16_type nRealDim = mesh_entity_type::nRealDim;
    using vector_realdim_type = Eigen::Matrix<double, nRealDim, 1>;

    MeshPrimitiveInfo(mesh_entity_type const& meshEntity)
        : M_meshEntity(meshEntity)
    {
        auto verticesUblas = meshEntity.vertices();
        auto G = em_cmatrix_col_type<double>(verticesUblas.data().begin(), nRealDim, mesh_entity_type::numVertices);
        M_bound_min = G.rowwise().minCoeff();
        M_bound_max = G.rowwise().maxCoeff();
        M_bound_min.array() -= 2 * FLT_MIN;
        M_bound_max.array() += 2 * FLT_MIN;
        auto bary = meshEntity.barycenter();
        M_centroid = Eigen::Map<Eigen::Matrix<double, nRealDim, 1>>(bary.data().begin());
    }
    MeshPrimitiveInfo(MeshPrimitiveInfo&&) = default;
    MeshPrimitiveInfo(MeshPrimitiveInfo const&) = default;
    MeshPrimitiveInfo& operator=(MeshPrimitiveInfo&&) = default;
    MeshPrimitiveInfo& operator=(MeshPrimitiveInfo const&) = default;

    mesh_entity_type const& meshEntity() const { return M_meshEntity.get(); }
    vector_realdim_type const& boundMin() const noexcept { return M_bound_min; }
    vector_realdim_type const& boundMax() const noexcept { return M_bound_max; }
    vector_realdim_type const& centroid() const noexcept { return M_centroid; }

private:
    vector_realdim_type M_bound_min;
    vector_realdim_type M_bound_max;
    vector_realdim_type M_centroid;
    std::reference_wrapper<mesh_entity_type const> M_meshEntity;
};


double calculateStepRays(int n, double start, double end) {
    double n_total = std::sqrt(n);
    return (end - start) / (n_total - 1);
}


template <typename RangeType2, typename RangeType>
void distScanToBoundary(RangeType2 const& rangeMesh, RangeType const& range, int number_rays_desired, std::ofstream& file, bool isViewInfo)
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

    std::string filenameB = "all_reults_per_vertex.csv";
    if (remove(filenameB.c_str()) != 0) {
        std::cerr << "Error delete file." << std::endl;
    }
    //std::ofstream myfileB (filenameB);
    std::ofstream myfileB(filenameB, std::ios::app);
    myfileB << "Num Vertex,distanceMinREAL,distanceFastMarching,errFastMarching,distanceMinCPU,errCPU,distanceMinGPU,errGPU,time_RT_CPU,time_RT_GPU\n";

    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    using mesh_entity_type2 = std::remove_const_t<entity_range_t<RangeType2>>;


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
        //if ( (pt0[0]>0.0f) && (pt0[0]<1.0f) && (pt0[1]>0.0f) && (pt0[1]<1.0f) && (pt0[2]>0.0f) && (pt0[2]<1.0f)) { ok = true; }
        //if ( (ray_origin[0]>0.0f) && (ray_origin[0]<1.0f) && (ray_origin[1]>0.0f) && (ray_origin[1]<1.0f) && (ray_origin[2]>0.0f) && (ray_origin[2]<1.0f)) { ok = true; }

        ok = true; // All points

        if (ok) {
            nbValues++;
            //distanceMinREAL= std::min( pt0[0], std::min(1.0-pt0[0], std::min(pt0[1], std::min(1.0-pt0[1], std::min(pt0[2],1.0-pt0[2])))));
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

            //file << "Number of effective rays=" << rays.size() << "\n";

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

            double errFastMarching = abs(distanceFastMarching[k] - distanceMinREAL);

            myfileB << k << ","
                << distanceMinREAL << ","
                << distanceFastMarching[k] << ","
                << errFastMarching << ","
                << distanceMinCPU << ","
                << errCPU << ","
                << distanceMinGPU << ","
                << errGPU << ","
                << t_laps_CPU << ","
                << t_laps_GPU << "\n";

            if (isViewInfo) std::cout << "[INFO] [" << k
                << "]"
                << "<" << ray_origin[0]
                << "," << ray_origin[1]
                << "," << ray_origin[2]
                << ">"
                << " Distance Min REAL=" << distanceMinREAL
                << " FastMarching=" << distanceFastMarching[k]
                << " err=" << errFastMarching
                << " CPU=" << distanceMinCPU
                << " err=" << errCPU
                << " GPU=" << distanceMinGPU
                << " err=" << errGPU
                << " t_laps_CPU=" << t_laps_CPU
                << " t_laps_GPU=" << t_laps_GPU
                << "\n";

            // Memory cleaning
            distance_CPU_mode.clear();
            distance_GPU_mode.clear();
        }
    } // END for k


    t_laps_CPU = std::chrono::duration_cast<std::chrono::microseconds>(t_end_bvh_cpu - t_begin_cpu).count();
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside BVH CPU : " << t_laps_CPU << " us\n";
    file << "TimeBVHCPU=" << t_laps_CPU << "\n";

    t_laps_CPU = t_laps_CPU_Total / nbValues;
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside Ray Tracing CPU : " << t_laps_CPU << " us\n";
    file << "TimeRaytracingCPU=" << t_laps_CPU << "\n";

    t_laps_GPU = std::chrono::duration_cast<std::chrono::microseconds>(t_end_bvh_gpu - t_begin_gpu).count();
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside BVH GPU : " << t_laps_GPU << " us \n";
    file << "TimeBVHGPU=" << t_laps_GPU << "\n";

    t_laps_GPU = t_laps_GPU_Total / nbValues;
    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside Ray Tracing GPU : " << t_laps_GPU << " us\n";
    file << "TimeRaytracingGPU=" << t_laps_GPU << "\n";

    if (isViewInfo) std::cout << "[INFO] Elapsed microseconds inside FastMarching : " << t_laps_FastMarching << " us\n";
    file << "TimeFastMarching=" << t_laps_FastMarching << "\n";

    myfileB.close();
}


BOOST_AUTO_TEST_SUITE(bvh_intersection_gpu_tests)


BOOST_AUTO_TEST_CASE(test_load_mesh3)
{

    bool isViewInfo = true;  //isViewInfo = false;;

    std::string filenameA = "results.csv";
    if (remove(filenameA.c_str()) != 0) {
        std::cerr << "Error delete file." << std::endl;
    }

    std::ofstream myfileA(filenameA);
    //std::ofstream myfileA(filenameA, std::ios::app);

    double hsize = 1.0f / 2.0f;
    int number_rays_desired = 703;


    using namespace Feel;
    using Feel::cout;
    using mesh_type = Mesh<Simplex<3, 1, 3>>; //<Dim,Order,RDim>
    auto mesh = unitCube(hsize);
    if (isViewInfo)
    {
        std::cout << "[INFO] maxNumElement : " << mesh->maxNumElements() << std::endl;
        std::cout << "[INFO] maxNumFace    : " << mesh->maxNumFaces() << std::endl;
        std::cout << "[INFO] maxNumPoints  : " << mesh->maxNumPoints() << std::endl;
        std::cout << "[INFO] maxNumVerices : " << mesh->maxNumVertices() << std::endl;
    }
    //auto rangeFaces = markedfaces(mesh);

    auto rangeFaces = boundaryfaces( mesh ) ;


    auto rangeElements = markedelements(mesh);

    myfileA << "h=" << hsize << "\n";
    myfileA << "number_rays_desired=" << number_rays_desired << "\n";

    auto Vh = Pch<1>(mesh);
    auto exp = exporter( _mesh = mesh, _name = fmt::format( "distance_{}d_o{}", 3, 1 ) );
    exp->addRegions();

    std::chrono::steady_clock::time_point t_begin_FastMarching, t_end_FastMarching;
    t_begin_FastMarching = std::chrono::steady_clock::now();
        //auto distToBoundary = distanceToRange(_space = Vh, _range = rangeFaces);
        auto distToBoundary = distanceToRange( _space=Vh, _range=boundaryfaces( mesh )  );
    t_end_FastMarching = std::chrono::steady_clock::now();
    t_laps_FastMarching = std::chrono::duration_cast<std::chrono::microseconds>(t_end_FastMarching - t_begin_FastMarching).count();
    exp->save();

    //auto smesh = createSubmesh(_mesh=mesh,_range=rangeElements);
    for (auto const& pointPair : mesh->points())
    {
        auto const& point = pointPair.second;
        auto const& coords = point.node();
        allNodeCoordinates.push_back({ coords[0], coords[1], coords[2] });
    }

    /*
    All points
    for (auto const& face : rangeElements)
    {
        for (auto const& point : face.get().points())
        {
            auto const& coords = point->node();
            allCoordinates.push_back({coords[0], coords[1], coords[2]});
        }
    }
    */

    for (int i = 0; i < allNodeCoordinates.size(); ++i)
    {
        distanceFastMarching.push_back(distToBoundary[i]);
    }


    // Test distance if it is ok
    for (int i = 0; i < allNodeCoordinates.size(); ++i)
    {

        double distanceMinREAL = std::min(allNodeCoordinates[i][0],
            std::min(1.0 - allNodeCoordinates[i][0],
                std::min(allNodeCoordinates[i][1],
                    std::min(1.0 - allNodeCoordinates[i][1],
                        std::min(allNodeCoordinates[i][2],
                            1.0 - allNodeCoordinates[i][2])))));
        std::cout << "["
            << i
            << "]"
            << " Point <"
            << allNodeCoordinates[i][0] << ", "
            << allNodeCoordinates[i][1] << ", "
            << allNodeCoordinates[i][2] << "> = "

            << distToBoundary[i]
            << " : "
            << distanceMinREAL
            << std::endl;
    }

    distScanToBoundary(rangeElements, rangeFaces, number_rays_desired, myfileA, isViewInfo);
    myfileA.close();
    std::cout << "\n";
}




BOOST_AUTO_TEST_SUITE_END()




