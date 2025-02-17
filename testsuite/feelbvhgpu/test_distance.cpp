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

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcore/json.hpp>

#include <feel/feelfilters/unitcube.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feells/distancetorange.hpp>

#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"

#include "thrust/device_vector.h"
#include "thrust/execution_policy.h"
#include "thrust/functional.h"
#include "thrust/host_vector.h"
#include "thrust/random.h"
#include "thrust/sort.h"
#include "thrust/transform.h"

#include "thrust/copy.h"
#include "thrust/count.h"
#include "thrust/generate.h"
#include "thrust/sort.h"

#include <hwloc.h>

// #include <rccl.h> // For multi-GPU not ready yet
// #include <roctx.h> //Scan Perf not ready yet

using namespace Feel;

inline AboutData
makeAbout()
{
    AboutData about( "test_distance BVH RT CPU AND GPU",
                     "test_distance",
                     "0.1",
                     "nD(n=3)",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2024 Feel++ Consortium" );

    about.addAuthor( "Noname", "developer", "Noname@cemosis.fr", "" );
    return about;
}

inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )( "number_rays_desired", po::value<int>()->default_value( 703 ), "nbRays" )( "isViewInfo", po::value<bool>()->default_value( true ), "isViewInfo" );
    return opts;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

struct DataDistanceErrTime
{
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

struct DataDistanceErrTimeAll
{
    int id;
    double distanceMinREAL;
    double distanceFastMarching;
    double errFastMarching;
    double distanceMinCPU;
    double errCPU;
    double distanceMinGPU;
    double errGPU;
};

struct DataTimeLapsConfig
{
    int nbRays; // <= uniform distribution
    int nbRaysDesired;
    double hsize;
    long int t_laps_BVH_CPU;
    long int t_laps_BVH_GPU;
    long int t_laps_RT_CPU;
    long int t_laps_RT_GPU;
    long int t_laps_FastMarching;
};

struct StatsResult
{
    double mean;
    double variance;
    double stdDev;
};

__global__ void onKernelNothing( float4* nothing )
{
    // nothing void
}

void runPreheatingGPU( int numDevice )
{

    int nbDevices = 0;
    hipGetDeviceCount( &nbDevices );
    if ( numDevice > nbDevices ) numDevice = 0;
    hipSetDevice( numDevice );

    float4* d_nothing;
    hipMalloc( &d_nothing, 14 * sizeof( float4 ) );
    onKernelNothing<<<1, 1>>>( d_nothing );
    hipFree( d_nothing );
}

void runScanPreheatingGPU()
{
    int nDevices;
    hipGetDeviceCount( &nDevices );
    for ( int i = 0; i < nDevices; ++i )
    {
        hipSetDevice( i );
        float4* d_nothing;
        hipMalloc( &d_nothing, 100 * sizeof( float4 ) );
        onKernelNothing<<<1, 1>>>( d_nothing );
        hipFree( d_nothing );
    }
}

template <typename BvhType, typename RayIntersectionResultType>
std::vector<double> getAllDistanceRayIntersections( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    std::vector<double> distance;
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            // std::cout << " --  Distance: " << rir.distance()<< "\n";
            distance.push_back( rir.distance() );
        }
        bvh->worldComm().barrier();
        // std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    return distance;
}

template <typename BvhType, typename RayIntersectionResultType>
std::vector<int> getId( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    std::vector<int> id;
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            // std::cout << " --  Distance: " << rir.distance()<< "\n";
            id.push_back( rir.get_id() );
        }
        bvh->worldComm().barrier();
        // std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    return id;
}

Eigen::Vector3d sphericalToCartesian( double r, double theta, double alpha )
{
    Eigen::Vector3d position;
    position[0] = r * sin( theta ) * cos( alpha );
    position[1] = r * sin( theta ) * sin( alpha );
    position[2] = r * cos( theta );
    return position;
}

double calculateStepRays( int n, double start, double end )
{
    double n_total = std::sqrt( n );
    return ( end - start ) / ( n_total - 1 );
}

template <typename RangeType2>
void distToBoundaryBVHpu(
    RangeType2 const& range,
    DataTimeLapsConfig& allDataPU,
    std::vector<std::vector<double>>& allNodeCoordinates,
    std::vector<DataDistanceErrTime>& allDataDistanceBVHRT,
    bool isViewInfo )
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
    auto bvhThirdParty = boundingVolumeHierarchy( _range = range, _kind = "third-party", _quality = BVHEnum::Quality::High );
    t_end_bvh_cpu = std::chrono::steady_clock::now();
    // END::Build BVH CPU

    sleep( 1 );

    // BEGIN::Build BVH GPU
    t_begin_gpu = std::chrono::steady_clock::now();
    auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-party" );
    t_end_bvh_gpu = std::chrono::steady_clock::now();
    // END::Build BVH GPU

    sleep( 1 );

    // isViewInfo=true;

    t_laps_CPU = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_cpu - t_begin_cpu ).count();
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH CPU : " << t_laps_CPU << " us\n";

    t_laps_GPU = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_gpu - t_begin_gpu ).count();
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH GPU : " << t_laps_GPU << " us\n";

    const double epsilon = 0.00001f;
    double distanceMinCPU = 0.0f;
    double distanceMinGPU = 0.0f;
    double distanceMinREAL = 0.0f;
    int nbValues = 0;

    for ( int k = 0; k < allNodeCoordinates.size(); ++k )
    {
        // Build ray
        std::vector<bvh_ray_type> rays;
        std::vector<double> distance_Real_mode;
        BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
        Eigen::Vector3d ray_origin = { allNodeCoordinates[k][0], allNodeCoordinates[k][1], allNodeCoordinates[k][2] };

        bool ok = false;
        // if ( (ray_origin[0]>0.0f) && (ray_origin[0]<1.0f) && (ray_origin[1]>0.0f) && (ray_origin[1]<1.0f) && (ray_origin[2]>0.0f) && (ray_origin[2]<1.0f)) { ok = true; }

        ok = true; // All points

        if ( ok )
        {
            nbValues++;
            distanceMinREAL = std::min( ray_origin[0], std::min( 1.0 - ray_origin[0], std::min( ray_origin[1], std::min( 1.0 - ray_origin[1], std::min( ray_origin[2], 1.0 - ray_origin[2] ) ) ) ) );

            double thetaStart = 0.0f;
            double thetaEnd = M_PI;
            double alphaStart = 0.0f;
            double alphaEnd = 2.0f * M_PI;
            double thetaStep = calculateStepRays( number_rays_desired, thetaStart, thetaEnd );
            double alphaStep = calculateStepRays( number_rays_desired, alphaStart, alphaEnd );

            for ( double theta = thetaStart; theta <= thetaEnd; theta += thetaStep )
            {
                for ( double alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep )
                {
                    Eigen::Vector3d ray_direction = sphericalToCartesian( 1.0f, theta, alpha );
                    Eigen::Vector3d ray_origin_Epsilon = ray_origin;
                    rays.push_back( bvh_ray_type( ray_origin_Epsilon, ray_direction ) );
                }
            }

            nbRays = rays.size();

            for ( int i = 0; i < rays.size(); ++i )
            {
                raysDistributed.push_back( rays[i] );
            }

            std::vector<double> dist;
            std::vector<double> distance_CPU_mode;
            std::vector<double> distance_GPU_mode;

            // In normal CPU mode
            // Ray Tracing
            t_begin_raytracing_cpu = std::chrono::steady_clock::now();
            auto multiRayDistributedIntersectionResult = bvhThirdParty->intersect( _ray = raysDistributed );
            for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
            {
                dist = getAllDistanceRayIntersections( bvhThirdParty, rayIntersectionResult );
                distance_CPU_mode.insert( distance_CPU_mode.end(), dist.begin(), dist.end() );
            }
            t_end_raytracing_cpu = std::chrono::steady_clock::now();

            t_begin_raytracing_gpu = std::chrono::steady_clock::now();
            auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect( _ray = raysDistributed );
            for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
            {
                dist = getAllDistanceRayIntersections( bvhHIPParty, rayIntersectionResult );
                distance_GPU_mode.insert( distance_GPU_mode.end(), dist.begin(), dist.end() );
            }
            t_end_raytracing_gpu = std::chrono::steady_clock::now();

            distanceMinCPU = INFINITY;
            for ( int k = 0; k < distance_CPU_mode.size(); ++k )
            {
                distanceMinCPU = fmin( distanceMinCPU, distance_CPU_mode[k] );
            }

            distanceMinGPU = INFINITY;

            for ( int k = 0; k < distance_GPU_mode.size(); ++k )
            {
                distanceMinGPU = fmin( distanceMinGPU, distance_GPU_mode[k] );
            }

            t_laps_CPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_raytracing_cpu - t_begin_raytracing_cpu ).count();
            t_laps_GPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_raytracing_gpu - t_begin_raytracing_gpu ).count();

            t_laps_CPU_Total = t_laps_CPU_Total + t_laps_CPU;
            t_laps_GPU_Total = t_laps_GPU_Total + t_laps_GPU;

            double errCPU = abs( distanceMinCPU - distanceMinREAL );
            double errGPU = abs( distanceMinGPU - distanceMinREAL );

            // double errFastMarching = abs(distanceFastMarching[k] - distanceMinREAL);

            // if (true) std::cout << "[" << k<< "]\n";

            if ( isViewInfo ) std::cout << "[INFO] [" << k
                                        << "]"
                                        << "<" << std::fixed << std::setprecision( 9 ) << ray_origin[0]
                                        << "," << std::fixed << std::setprecision( 9 ) << ray_origin[1]
                                        << "," << std::fixed << std::setprecision( 9 ) << ray_origin[2]
                                        << ">"
                                        << " Distance Min REAL=" << distanceMinREAL
                                        //<< " FastMarching=" << distanceFastMarching[k]
                                        //<< " err=" << errFastMarching
                                        << " CPU=" << std::fixed << std::setprecision( 9 ) << distanceMinCPU
                                        << " err=" << std::fixed << std::setprecision( 9 ) << errCPU
                                        << " GPU=" << std::fixed << std::setprecision( 9 ) << distanceMinGPU
                                        << " err=" << std::fixed << std::setprecision( 9 ) << errGPU
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
                t_laps_GPU };

            allDataDistanceBVHRT.push_back( data );

            // Memory cleaning
            distance_CPU_mode.clear();
            distance_GPU_mode.clear();
        }
    } // END for k

    // Elapse Time BVH - RT - CPU - GPU
    allDataPU.t_laps_BVH_CPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_bvh_cpu - t_begin_cpu ).count();
    // allDataPU.t_laps_RT_CPU = t_laps_CPU_Total / nbValues;
    allDataPU.t_laps_RT_CPU = t_laps_CPU_Total;
    allDataPU.t_laps_BVH_GPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_bvh_gpu - t_begin_gpu ).count();
    allDataPU.t_laps_RT_GPU = t_laps_GPU_Total;
    // allDataPU.t_laps_RT_GPU = t_laps_GPU_Total / nbValues;
    allDataPU.nbRays = nbRays;
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH CPU : " << allDataPU.t_laps_BVH_CPU << " ms\n";
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside Ray Tracing CPU : " << allDataPU.t_laps_RT_CPU << " ms\n";
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH GPU : " << allDataPU.t_laps_BVH_GPU << " ms\n";
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside Ray Tracing GPU : " << allDataPU.t_laps_RT_GPU << " ms\n";
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside FastMarching : " << allDataPU.t_laps_FastMarching << " ms\n";
}

template <typename RangeType2>
void distToBoundaryBVHpuSendAllNode(
    RangeType2 const& range,
    DataTimeLapsConfig& allDataPU,
    std::vector<std::vector<double>>& allNodeCoordinates,
    std::vector<DataDistanceErrTimeAll>& allDataDistanceBVHRTAll,
    bool isViewInfo )
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

    int number_rays_desired = allDataPU.nbRaysDesired;

    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType2>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    // BEGIN::Build BVH CPU
    t_begin_cpu = std::chrono::steady_clock::now();
    tic();
    auto bvhThirdParty = boundingVolumeHierarchy( _range = range, _kind = "third-party", _quality = BVHEnum::Quality::High );
    auto timeBVHCPUDuration = toc( "timeBVHCPUDuration" );
    t_end_bvh_cpu = std::chrono::steady_clock::now();
    // END::Build BVH CPU

    sleep( 1 );

    // BEGIN::Build BVH GPU
    t_begin_gpu = std::chrono::steady_clock::now();
    tic();
    auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-party" );
    auto timeBVHGPUDuration = toc( "timeBVHGPUDuration" );
    t_end_bvh_gpu = std::chrono::steady_clock::now();
    // END::Build BVH GPU

    sleep( 1 );

    // BVH
    t_laps_CPU = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_cpu - t_begin_cpu ).count();
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH CPU : " << t_laps_CPU << " us\n";

    t_laps_GPU = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_gpu - t_begin_gpu ).count();
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH GPU : " << t_laps_GPU << " us\n";

    const double epsilon = 0.00001f;
    int nbValues = 0;

    std::vector<bvh_ray_type> rays;
    std::vector<double> distance_Real_mode;
    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;

    int nbRays = std::sqrt( number_rays_desired );
    nbRays = nbRays * nbRays;

    int kblock = 10000;

    for ( int k = 0; k < allNodeCoordinates.size(); ++k )
    {
        // Build ray
        Eigen::Vector3d ray_origin = { allNodeCoordinates[k][0], allNodeCoordinates[k][1], allNodeCoordinates[k][2] };
        double thetaStart = 0.0f;
        double thetaEnd = M_PI;
        double alphaStart = 0.0f;
        double alphaEnd = 2.0f * M_PI;
        double thetaStep = calculateStepRays( number_rays_desired, thetaStart, thetaEnd );
        double alphaStep = calculateStepRays( number_rays_desired, alphaStart, alphaEnd );

        int j = 0;
        for ( double theta = thetaStart; theta <= thetaEnd; theta += thetaStep )
        {
            for ( double alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep )
            {
                Eigen::Vector3d ray_direction = sphericalToCartesian( 1.0f, theta, alpha );
                Eigen::Vector3d ray_origin_Epsilon = ray_origin;
                rays.push_back( bvh_ray_type( ray_origin_Epsilon, ray_direction ) );
                rays.back().id = j + ( k + 1 ) * kblock;
                j++;
            }
        }
    } // END for k

    for ( int k = 0; k < rays.size(); ++k )
    {
        raysDistributed.push_back( rays[k] );
    }

    std::vector<double> dist;
    std::vector<double> distance_CPU_mode;
    std::vector<double> distance_GPU_mode;

    std::vector<int> id;
    std::vector<int> id_CPU;
    std::vector<int> id_GPU;

    // Ray Tracing BVH CPU
    t_begin_raytracing_cpu = std::chrono::steady_clock::now();
    tic();
    auto multiRayDistributedIntersectionResult = bvhThirdParty->intersect( _ray = raysDistributed );
    auto timeRTCPUDuration = toc( "timeRTCPUDuration" );

    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
    {
        dist = getAllDistanceRayIntersections( bvhThirdParty, rayIntersectionResult );
        distance_CPU_mode.insert( distance_CPU_mode.end(), dist.begin(), dist.end() );
        id = getId( bvhThirdParty, rayIntersectionResult );
        id_CPU.insert( id_CPU.end(), id.begin(), id.end() );
    }

    t_end_raytracing_cpu = std::chrono::steady_clock::now();
    // std::cout << "timeRTCPUDuration " << timeRTCPUDuration << " \n";

    // Ray Tracing BVH GPU
    t_begin_raytracing_gpu = std::chrono::steady_clock::now();
    tic();
    auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect( _ray = raysDistributed, _parallel = false );
    auto timeRTGPUDuration = toc( "timeRTGPUDuration" );

    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
    {
        dist = getAllDistanceRayIntersections( bvhHIPParty, rayIntersectionResult );
        distance_GPU_mode.insert( distance_GPU_mode.end(), dist.begin(), dist.end() );
        id = getId( bvhHIPParty, rayIntersectionResult );
        id_GPU.insert( id_GPU.end(), id.begin(), id.end() );
    }

    t_end_raytracing_gpu = std::chrono::steady_clock::now();
    // std::cout << "timeRTGPUDuration " << timeRTGPUDuration << " \n";

    int k1 = 0;
    int k2 = 0;
    std::vector<double> distanceMinCPU;
    std::vector<double> distanceMinGPU;
    double value = INFINITY;

    if ( true )
    {
        int idBlock = 1;
        int idBlockLast = 1;
        value = distance_CPU_mode[0];
        for ( int i = 0; i < distance_CPU_mode.size(); ++i )
        {
            idBlock = int( id_CPU[i] / kblock );
            if ( idBlock == idBlockLast )
            {
                // std::cout << "CPU NumRay "<<id_CPU[i]<<" blockid="<<idBlock<<" value="<<distance_CPU_mode[i]<<" \n";
                value = fmin( value, distance_CPU_mode[i] );
            }
            if ( idBlock != idBlockLast )
            {
                distanceMinCPU.push_back( value );
                // std::cout << "  CPU NumRay "<<id_CPU[i-1]<<" blockid="<<idBlockLast << " value="<<value<< " \n";
                if ( isViewInfo ) std::cout << "  CPU blockid=" << idBlockLast << " value=" << value << " \n";
                value = INFINITY;
                idBlockLast = idBlock;
            }
        }
        // std::cout << "  CPU NumRay "<<id_CPU[distance_CPU_mode.size()-1]<<" blockid="<<idBlockLast << " value="<<value<< " \n";
        if ( isViewInfo ) std::cout << "  CPU blockid=" << idBlockLast << " value=" << value << " \n";
        distanceMinCPU.push_back( value );

        if ( isViewInfo ) std::cout << "=====================================================================================\n";

        idBlock = 1;
        idBlockLast = 1;
        value = distance_GPU_mode[0];
        for ( int i = 0; i < distance_GPU_mode.size(); ++i )
        {
            idBlock = int( id_GPU[i] / kblock );

            if ( idBlock == idBlockLast )
            {
                // std::cout << "GPU NumRay "<<id_GPU[i]<<" blockid="<<idBlock<<" value="<<distance_GPU_mode[i]<<" \n";
                value = fmin( value, distance_GPU_mode[i] );
            }
            if ( idBlock != idBlockLast )
            {
                distanceMinGPU.push_back( value );
                // std::cout << "  GPU NumRay "<<id_GPU[i-1]<<" blockid="<<idBlockLast << " value="<<value<< " \n";
                if ( isViewInfo ) std::cout << "  GPU blockid=" << idBlockLast << " value=" << value << " \n";
                value = INFINITY;
                idBlockLast = idBlock;
            }
        }
        // std::cout << "  GPU NumRay "<<id_GPU[distance_GPU_mode.size()-1]<<" blockid="<<idBlockLast << " value="<<value<< " \n";
        if ( isViewInfo ) std::cout << "  GPU blockid=" << idBlockLast << " value=" << value << " \n";
        distanceMinGPU.push_back( value );

        std::cout << "Nb Coordinates = " << allNodeCoordinates.size() << "\n";

        for ( int index = 0; index < allNodeCoordinates.size(); ++index )
        {
            double distanceMinREAL;
            distanceMinREAL = std::min( allNodeCoordinates[index][0],
                                        std::min( 1.0 - allNodeCoordinates[index][0],
                                                  std::min( allNodeCoordinates[index][1],
                                                            std::min( 1.0 - allNodeCoordinates[index][1],
                                                                      std::min( allNodeCoordinates[index][2], 1.0 - allNodeCoordinates[index][2] ) ) ) ) );

            double errCPU = abs( distanceMinCPU[index] - distanceMinREAL );
            double errGPU = abs( distanceMinGPU[index] - distanceMinREAL );

            DataDistanceErrTimeAll data = {
                index,
                distanceMinREAL,
                -1,
                -1,
                distanceMinCPU[index],
                errCPU,
                distanceMinGPU[index],
                errGPU };
            allDataDistanceBVHRTAll.push_back( data );
        }
    }

    // Elapse Time BVH - RT - CPU - GPU
    allDataPU.t_laps_BVH_CPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_bvh_cpu - t_begin_cpu ).count();
    allDataPU.t_laps_RT_CPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_raytracing_cpu - t_begin_raytracing_cpu ).count();
    allDataPU.t_laps_BVH_GPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_bvh_gpu - t_begin_gpu ).count();
    allDataPU.t_laps_RT_GPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_raytracing_gpu - t_begin_raytracing_gpu ).count();
    allDataPU.nbRays = nbRays;

    if ( isViewInfo )
    {
        std::cout << "[INFO] Elapsed microseconds\n";
        std::cout << "[INFO] BVH CPU : " << allDataPU.t_laps_BVH_CPU << " ms\n";
        std::cout << "[INFO] RT  CPU : " << allDataPU.t_laps_RT_CPU << " ms\n";
        std::cout << "[INFO] BVH GPU : " << allDataPU.t_laps_BVH_GPU << " ms\n";
        std::cout << "[INFO] RT  GPU : " << allDataPU.t_laps_RT_GPU << " ms\n";
        std::cout << "[INFO] FastMarching : " << allDataPU.t_laps_FastMarching << " ms\n";
        std::cout << "[INFO] ========================================================================\n";
    }
}

// Statistical part mean, std,variance...
StatsResult calculateStats( const std::vector<long int>& values )
{
    StatsResult result;
    int n = values.size();
    double sum = 0.0;
    for ( const auto& value : values )
    {
        sum += value;
    }
    result.mean = sum / n;
    double squaredDiffSum = 0.0;
    for ( const auto& value : values )
    {
        double diff = value - result.mean;
        squaredDiffSum += diff * diff;
    }
    result.variance = squaredDiffSum / n;
    result.stdDev = std::sqrt( result.variance );
    return result;
}

// We combine the two calculations to make only one pass
void calculateCPUGPUStats( const std::vector<DataDistanceErrTime>& allDataDistanceBVHRT, StatsResult& cpuStats, StatsResult& gpuStats )
{
    std::vector<long int> cpuTimes, gpuTimes;

    for ( const auto& data : allDataDistanceBVHRT )
    {
        cpuTimes.push_back( data.t_laps_CPU );
        gpuTimes.push_back( data.t_laps_GPU );
    }
    cpuStats = calculateStats( cpuTimes );
    gpuStats = calculateStats( gpuTimes );
}

BOOST_AUTO_TEST_SUITE( distance_bvh_cpu_gpu_gpu_tests )

BOOST_AUTO_TEST_CASE( all_distance )
{

    // We read the value of "hsize" and "number_rays_desired"
    double hsize = option( _name = "hsize" ).as<double>();
    int number_rays_desired = option( _name = "number_rays_desired" ).as<int>();
    bool isViewInfo = option( _name = "isViewInfo" ).as<bool>();
    // double hsize = 1.0f / 2.0f;
    // int number_rays_desired = 703;

    // hsize = 0.25;
    // hsize = 0.005;

    // number_rays_desired =2000;

    bool isPreheating = true; // isPreheating = false;

    /*
    if (isPreheating) runPreheatingGPU(0);
    if (isPreheating) runPreheatingGPU(1);
    if (isPreheating) runPreheatingGPU(2);
    if (isPreheating) runPreheatingGPU(3);
    */

    if ( isPreheating ) runScanPreheatingGPU();

    using namespace Feel;
    using Feel::cout;

    // 3D object initialization
    using mesh_type = Mesh<Simplex<3, 1, 3>>; //<Dim,Order,RDim>

    // auto mesh = unitCube();

    auto mesh = unitCube( hsize );

    // Small information about the structure
    if ( isViewInfo )
    {
        std::cout << "[INFO] hsize : " << hsize << std::endl;
        std::cout << "[INFO] number_rays_desired : " << number_rays_desired << std::endl;

        std::cout << "[INFO] maxNumElement : " << mesh->maxNumElements() << std::endl;
        std::cout << "[INFO] maxNumFace    : " << mesh->maxNumFaces() << std::endl;
        std::cout << "[INFO] maxNumPoints  : " << mesh->maxNumPoints() << std::endl;
        std::cout << "[INFO] maxNumVerices : " << mesh->maxNumVertices() << std::endl;
    }

    // Selecting what you want to process
    auto rangeFaces = markedfaces( mesh );
    auto submeshFaces = boundaryfaces( mesh );
    auto rangeElements = markedelements( mesh );

    auto Vh = Pch<1>( mesh );
#if 0
    // old version to access the coordinates of the nodes of the points which is not the same for the distancetorange function
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
    allDataPU.nbRaysDesired = number_rays_desired;
    allDataPU.hsize = hsize;
    std::vector<DataDistanceErrTime> allDataDistanceBVHRT;

    // List of node coordinates
    for ( size_type k = 0; k < Vh->nLocalDofWithGhost(); ++k )
    {
        auto const& dofPt = Vh->dof()->dofPoint( k ).template get<0>();
        allNodeCoordinates.push_back( { dofPt[0], dofPt[1], dofPt[2] } );
    }

    int nbNode = allNodeCoordinates.size();

    // Calculates Node points to Surface distances by the method FastMarching
    std::chrono::steady_clock::time_point t_begin_FastMarching, t_end_FastMarching;
    t_begin_FastMarching = std::chrono::steady_clock::now();
    auto distToBoundary = distanceToRange( _space = Vh, _range = submeshFaces );
    t_end_FastMarching = std::chrono::steady_clock::now();
    long int t_laps_FastMarching = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_FastMarching - t_begin_FastMarching ).count();

    allDataPU.t_laps_FastMarching = t_laps_FastMarching;

    bool isTransferAllNodes = true;
    // isTransferAllNodes=false;

    if ( isTransferAllNodes )
    {
        // In this part all data is sent at once from CPU to GPU.
        std::vector<DataDistanceErrTimeAll> allDataDistanceBVHRTAll;

        // We calculate the distances from the edge of the cube and the intersection points.
        // distToBoundaryBVHpuSendAllNode(rangeFaces,allDataPU,allNodeCoordinates,allDataDistanceBVHRTAll,isViewInfo);
        distToBoundaryBVHpuSendAllNode( rangeFaces, allDataPU, allNodeCoordinates, allDataDistanceBVHRTAll, false );

        if ( true )
        {
            for ( int i = 0; i < nbNode; ++i )
            {
                allDataDistanceBVHRTAll[i].distanceFastMarching = distToBoundary[i];
                allDataDistanceBVHRTAll[i].errFastMarching = abs( distToBoundary[i] - allDataDistanceBVHRTAll[i].distanceMinREAL );
            }

            // Debriefing Save all data
            std::string filenameDataDistanceErrTime = "all_results_per_vertex.csv";
            if ( remove( filenameDataDistanceErrTime.c_str() ) != 0 )
            {
                std::cerr << "Error delete file." << std::endl;
            }
            std::ofstream myfileB( filenameDataDistanceErrTime );
            myfileB << "Num Vertex,PosX,PosY,PosZ,distanceMinREAL,distanceFastMarching,errFastMarching,distanceMinCPU,errCPU,distanceMinGPU,errGPU\n";
            for ( int i = 0; i < nbNode; ++i )
            {
                myfileB << allDataDistanceBVHRTAll[i].id << ","
                        << std::fixed << std::setprecision( 9 ) << allNodeCoordinates[i][0] << ","
                        << std::fixed << std::setprecision( 9 ) << allNodeCoordinates[i][1] << ","
                        << std::fixed << std::setprecision( 9 ) << allNodeCoordinates[i][2] << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRTAll[i].distanceMinREAL << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRTAll[i].distanceFastMarching << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRTAll[i].errFastMarching << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRTAll[i].distanceMinCPU << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRTAll[i].errCPU << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRTAll[i].distanceMinGPU << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRTAll[i].errGPU << "\n";
            }
            myfileB.close();

            std::string filenameA = "results.csv";
            if ( remove( filenameA.c_str() ) != 0 )
            {
                std::cerr << "Error delete file." << std::endl;
            }
            std::ofstream myfileA( filenameA );
            myfileA << "hsize=" << allDataPU.hsize << "\n";
            myfileA << "maxNumElement=" << mesh->maxNumElements() << "\n";
            myfileA << "maxNumFace=" << mesh->maxNumFaces() << "\n";
            myfileA << "maxNumPoints=" << mesh->maxNumPoints() << "\n";
            myfileA << "maxNumVerices=" << mesh->maxNumVertices() << "\n";
            myfileA << "nbRaysDesired=" << allDataPU.nbRaysDesired << "\n";
            myfileA << "nbRays=" << allDataPU.nbRays << "\n";

            myfileA << "timeBVHcpu=" << allDataPU.t_laps_BVH_CPU << "\n";
            myfileA << "timeRTcpu=" << allDataPU.t_laps_RT_CPU << "\n";

            myfileA << "timeBVHgpu=" << allDataPU.t_laps_BVH_GPU << "\n";
            myfileA << "timeRTgpu=" << allDataPU.t_laps_RT_GPU << "\n";

            myfileA << "timeFastMarching=" << allDataPU.t_laps_FastMarching << "\n";
            myfileA << "totalTimeBVHRTcpu=" << allDataPU.t_laps_BVH_CPU + allDataPU.t_laps_RT_CPU << "\n";
            myfileA << "totalTimeBVHRTgpu=" << allDataPU.t_laps_BVH_GPU + allDataPU.t_laps_RT_GPU << "\n";
            myfileA.close();

            std::string filenameC = "results2.csv";
            if ( remove( filenameC.c_str() ) != 0 )
            {
                std::cerr << "Error delete file." << std::endl;
            }
            std::ofstream myfileC( filenameC );
            myfileC << "hsize" << ",";
            myfileC << "maxNumElement" << ",";
            myfileC << "maxNumFace" << ",";
            myfileC << "maxNumPoints" << ",";
            myfileC << "maxNumVerices" << ",";
            myfileC << "nbRaysDesired" << ",";
            myfileC << "nbRays" << ",";

            myfileC << "timeBVHcpu" << ",";
            myfileC << "timeRTcpu" << ",";

            myfileC << "timeBVHgpu" << ",";
            myfileC << "timeRTgpu" << ",";

            myfileC << "timeFastMarching" << ",";
            myfileC << "totalTimeBVHRTcpu" << ",";
            myfileC << "totalTimeBVHRTgpu" << "\n";

            myfileC << allDataPU.hsize << ",";
            myfileC << mesh->maxNumElements() << ",";
            myfileC << mesh->maxNumFaces() << ",";
            myfileC << mesh->maxNumPoints() << ",";
            myfileC << mesh->maxNumVertices() << ",";
            myfileC << allDataPU.nbRaysDesired << ",";
            myfileC << allDataPU.nbRays << ",";

            myfileC << allDataPU.t_laps_BVH_CPU << ",";
            myfileC << allDataPU.t_laps_RT_CPU << ",";

            myfileC << allDataPU.t_laps_BVH_GPU << ",";
            myfileC << allDataPU.t_laps_RT_GPU << ",";

            myfileC << allDataPU.t_laps_FastMarching << ",";
            myfileC << allDataPU.t_laps_BVH_CPU + allDataPU.t_laps_RT_CPU << ",";
            myfileC << allDataPU.t_laps_BVH_GPU + allDataPU.t_laps_RT_GPU << "\n";

            myfileC.close();

            // Data backup file distances for paraview
            // Transferring data in the format for the export function
            auto distanceMinCPU = Vh->element();
            auto distanceMinGPU = Vh->element();
            for ( size_t i = 0; i < nbNode; ++i )
            {
                distanceMinCPU[i] = allDataDistanceBVHRTAll[i].distanceMinCPU;
                distanceMinGPU[i] = allDataDistanceBVHRTAll[i].distanceMinGPU;
            }
            // Save the file with all the distance parameters
            auto exp = exporter( _mesh = mesh, _name = fmt::format( "distance_{}d_o{}", 3, 1 ) );
            exp->addRegions();
            exp->add( "distToBoundary", distToBoundary );
            exp->add( "distanceMinCPU", distanceMinCPU );
            exp->add( "distanceMinGPU", distanceMinGPU );
            exp->save();
        }
    }
    else
    {
        // Calculates Node points to Surface distances by the method BVH RT CPU and GPU
        isViewInfo = false;
        distToBoundaryBVHpu( rangeFaces, allDataPU, allNodeCoordinates, allDataDistanceBVHRT, isViewInfo );

        if ( true )
        {
            // We fill the Fast Marching distance data into the allDataDistanceBVHRT data structure
            for ( int i = 0; i < nbNode; ++i )
            {
                allDataDistanceBVHRT[i].distanceFastMarching = distToBoundary[i];
                allDataDistanceBVHRT[i].errFastMarching = abs( distToBoundary[i] - allDataDistanceBVHRT[i].distanceMinREAL );
            }

            // Debriefing Save all data
            std::string filenameDataDistanceErrTime = "all_results_per_vertex.csv";
            if ( remove( filenameDataDistanceErrTime.c_str() ) != 0 )
            {
                std::cerr << "Error delete file." << std::endl;
            }
            std::ofstream myfileB( filenameDataDistanceErrTime );
            myfileB << "Num Vertex,PosX,PosY,PosZ,distanceMinREAL,distanceFastMarching,errFastMarching,distanceMinCPU,errCPU,distanceMinGPU,errGPU,time_RT_CPU,time_RT_GPU\n";
            for ( int i = 0; i < nbNode; ++i )
            {
                myfileB << allDataDistanceBVHRT[i].id << ","
                        << std::fixed << std::setprecision( 9 ) << allNodeCoordinates[i][0] << ","
                        << std::fixed << std::setprecision( 9 ) << allNodeCoordinates[i][1] << ","
                        << std::fixed << std::setprecision( 9 ) << allNodeCoordinates[i][2] << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRT[i].distanceMinREAL << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRT[i].distanceFastMarching << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRT[i].errFastMarching << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRT[i].distanceMinCPU << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRT[i].errCPU << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRT[i].distanceMinGPU << ","
                        << std::fixed << std::setprecision( 9 ) << allDataDistanceBVHRT[i].errGPU << ","
                        << allDataDistanceBVHRT[i].t_laps_CPU << ","
                        << allDataDistanceBVHRT[i].t_laps_GPU << "\n";
            }
            myfileB.close();

            // Statistical part
            StatsResult cpuStats, gpuStats;
            calculateCPUGPUStats( allDataDistanceBVHRT, cpuStats, gpuStats );

            std::string filenameA = "results.csv";
            if ( remove( filenameA.c_str() ) != 0 )
            {
                std::cerr << "Error delete file." << std::endl;
            }
            std::ofstream myfileA( filenameA );
            myfileA << "hsize=" << allDataPU.hsize << "\n";
            myfileA << "maxNumElement=" << mesh->maxNumElements() << "\n";
            myfileA << "maxNumFace=" << mesh->maxNumFaces() << "\n";
            myfileA << "maxNumPoints=" << mesh->maxNumPoints() << "\n";
            myfileA << "maxNumVerices=" << mesh->maxNumVertices() << "\n";
            myfileA << "nbRaysDesired=" << allDataPU.nbRaysDesired << "\n";
            myfileA << "nbRays=" << allDataPU.nbRays << "\n";

            myfileA << "timeBVHcpu=" << allDataPU.t_laps_BVH_CPU << "\n";
            myfileA << "timeMeanRTcpu=" << allDataPU.t_laps_RT_CPU << "\n";
            // myfileA << "timeStandardDeviationRTcpu=" << cpuStats.stdDev  << "\n";
            // myfileA << "timeVarianceRTcpu=" << cpuStats.variance << "\n";

            myfileA << "timeBVHgpu=" << allDataPU.t_laps_BVH_GPU << "\n";
            myfileA << "timeMeanRTgpu=" << allDataPU.t_laps_RT_GPU << "\n";
            // myfileA << "timeStandardDeviationRTgpu=" << gpuStats.stdDev << "\n";
            // myfileA << "timeVarianceRTgpu=" << gpuStats.variance << "\n";

            myfileA << "timeFastMarching=" << allDataPU.t_laps_FastMarching << "\n";
            myfileA << "totalTimeBVHRTcpu=" << allDataPU.t_laps_BVH_CPU + allDataPU.t_laps_RT_CPU << "\n";
            myfileA << "totalTimeBVHRTgpu=" << allDataPU.t_laps_BVH_GPU + allDataPU.t_laps_RT_GPU << "\n";
            myfileA.close();

            std::string filenameC = "results2.csv";
            if ( remove( filenameC.c_str() ) != 0 )
            {
                std::cerr << "Error delete file." << std::endl;
            }
            std::ofstream myfileC( filenameA );
            myfileC << "hsize=" << allDataPU.hsize << "\n";
            myfileC << "maxNumElement=" << mesh->maxNumElements() << "\n";
            myfileC << "maxNumFace=" << mesh->maxNumFaces() << "\n";
            myfileC << "maxNumPoints=" << mesh->maxNumPoints() << "\n";
            myfileC << "maxNumVerices=" << mesh->maxNumVertices() << "\n";
            myfileC << "nbRaysDesired=" << allDataPU.nbRaysDesired << "\n";
            myfileC << "nbRays=" << allDataPU.nbRays << "\n";

            myfileC << "timeBVHcpu,";
            myfileC << "timeMeanRTcpu,";
            // myfileC << "timeStandardDeviationRTcpu=" << cpuStats.stdDev  << "\n";
            // myfileC << "timeVarianceRTcpu=" << cpuStats.variance << "\n";

            myfileC << "timeBVHgpun,";
            myfileC << "timeMeanRTgpu,";
            // myfileA << "timeStandardDeviationRTgpu=" << gpuStats.stdDev << "\n";
            // myfileA << "timeVarianceRTgpu=" << gpuStats.variance << "\n";

            myfileC << "timeFastMarching,";
            myfileC << "totalTimeBVHRTcpu,";
            myfileC << "totalTimeBVHRTgpu\n";

            myfileC << allDataPU.hsize << ",";
            myfileC << mesh->maxNumElements() << ",";
            myfileC << mesh->maxNumFaces() << ",";
            myfileC << mesh->maxNumPoints() << ",";
            myfileC << mesh->maxNumVertices() << ",";
            myfileC << allDataPU.nbRaysDesired << ",";
            myfileC << allDataPU.nbRays << ",";

            myfileC << allDataPU.t_laps_BVH_CPU << ",";
            myfileC << allDataPU.t_laps_RT_CPU << ",";
            // myfileC << "timeStandardDeviationRTcpu=" << cpuStats.stdDev  <<  ",";
            // myfileC << "timeVarianceRTcpu=" << cpuStats.variance << ",";

            myfileC << allDataPU.t_laps_BVH_GPU << ",";
            myfileC << allDataPU.t_laps_RT_GPU << ",";
            // myfileA << "timeStandardDeviationRTgpu=" << gpuStats.stdDev <<  ",";
            // myfileA << "timeVarianceRTgpu=" << gpuStats.variance << ",";

            myfileC << allDataPU.t_laps_FastMarching << ",";
            myfileC << allDataPU.t_laps_BVH_CPU + allDataPU.t_laps_RT_CPU << ",";
            myfileC << allDataPU.t_laps_BVH_GPU + allDataPU.t_laps_RT_GPU << "\n";
            myfileC.close();

            // Data backup file distances for paraview
            // Transferring data in the format for the export function
            auto distanceMinCPU = Vh->element();
            auto distanceMinGPU = Vh->element();
            for ( size_t i = 0; i < nbNode; ++i )
            {
                distanceMinCPU[i] = allDataDistanceBVHRT[i].distanceMinCPU;
                distanceMinGPU[i] = allDataDistanceBVHRT[i].distanceMinGPU;
            }
            // Save the file with all the distance parameters
            auto exp = exporter( _mesh = mesh, _name = fmt::format( "distance_{}d_o{}", 3, 1 ) );
            exp->addRegions();
            exp->add( "distToBoundary", distToBoundary );
            exp->add( "distanceMinCPU", distanceMinCPU );
            exp->add( "distanceMinGPU", distanceMinGPU );
            exp->save();

            // Memory cleaning
            allNodeCoordinates.clear();
            allDataDistanceBVHRT.clear();
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
