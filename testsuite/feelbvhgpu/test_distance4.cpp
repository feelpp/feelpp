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

#include <signal.h>

// #include <rccl.h> // For multi-GPU not ready yet
// #include <roctx.h> //Scan Perf not ready yet

using namespace Feel;

void sigterm_handler( int signum )
{
    MPI_Abort( MPI_COMM_WORLD, 0 );
    exit( 0 );
}

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


struct DataDistanceErrTimeAll
{
    int rank;
    size_t id;
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
    int rank;
    size_t nbRays; // <= uniform distribution
    size_t nbRaysDesired;
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
            distance.push_back( rir.distance() );
        }
        bvh->worldComm().barrier();
        // std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    return distance;
}

template <typename BvhType, typename RayIntersectionResultType>
std::vector<size_t> getId( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    std::vector<size_t> id;
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            // std::cout << " RIR --  id " <<"["<<bvh->worldComm().rank()<<"] = "<< rir.get_id()<< "\n";
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

template <typename BVHRayType>
std::vector<BVHRayType> generateRays( const std::vector<std::vector<double>>& allNodeCoordinates, int number_rays_desired )
{
    std::vector<BVHRayType> rays;
    int nbRays = std::sqrt( number_rays_desired );
    nbRays = nbRays * nbRays;
    const size_t kblock = nbRays + 1;

    for ( size_t k = 0; k < allNodeCoordinates.size(); ++k )
    {
        Eigen::Vector3d ray_origin = { allNodeCoordinates[k][0], allNodeCoordinates[k][1], allNodeCoordinates[k][2] };
        double thetaStart = 0.0;
        double thetaEnd = M_PI;
        double alphaStart = 0.0;
        double alphaEnd = 2.0 * M_PI;
        double thetaStep = calculateStepRays( number_rays_desired, thetaStart, thetaEnd );
        double alphaStep = calculateStepRays( number_rays_desired, alphaStart, alphaEnd );

        size_t j = 0;
        for ( double theta = thetaStart; theta <= thetaEnd; theta += thetaStep )
        {
            for ( double alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep )
            {
                Eigen::Vector3d ray_direction = sphericalToCartesian( 1.0, theta, alpha );
                Eigen::Vector3d ray_origin_Epsilon = ray_origin;
                rays.push_back( BVHRayType( ray_origin_Epsilon, ray_direction ) );
                rays.back().id = j + ( k + 1 ) * kblock;
                j++;
            }
        }
    }

    return rays;
}

double valueFilter( double v )
{
    if ( fabs( v ) < 1.0e-10 ) return 0.0;
    return v;
}

double calculateDist( std::vector<std::vector<double>> allNodeCoordinates, int index )
{
    double distanceMinREAL;
    distanceMinREAL = std::min( allNodeCoordinates[index][0],
                                std::min( 1.0 - allNodeCoordinates[index][0],
                                          std::min( allNodeCoordinates[index][1],
                                                    std::min( 1.0 - allNodeCoordinates[index][1],
                                                              std::min( allNodeCoordinates[index][2], 1.0 - allNodeCoordinates[index][2] ) ) ) ) );

    distanceMinREAL = valueFilter( distanceMinREAL );
    return ( distanceMinREAL );
}

std::vector<double> calculateDistanceMinPU(
    const std::vector<size_t>& id_PU,
    const std::vector<double>& distance_PU_mode,
    size_t kblock,
    const std::vector<std::vector<double>>& allNodeCoordinates,
    bool isViewValues = false )
{
    std::vector<double> distanceMinPU;
    size_t idBlock = 1;
    size_t idBlockLast = 1;
    double value = distance_PU_mode[0];

    for ( size_t i = 0; i < distance_PU_mode.size(); ++i )
    {
        idBlock = static_cast<size_t>( id_PU[i] / kblock );
        if ( idBlock == idBlockLast )
        {
            value = std::min( value, distance_PU_mode[i] );
        }
        if ( idBlock != idBlockLast )
        {
            value = valueFilter( value );
            distanceMinPU.push_back( value );
            if ( isViewValues )
                std::cout << "  Block Point id=" << idBlockLast << " value =" << value << "  <>  " << calculateDist( allNodeCoordinates, idBlockLast - 1 ) << " \n";
            value = std::numeric_limits<double>::infinity();
            idBlockLast = idBlock;
        }
    }

    value = valueFilter( value );
    distanceMinPU.push_back( value );
    if ( isViewValues )
        std::cout << "  Block Point id=" << idBlockLast << " value =" << value << "  <>  " << calculateDist( allNodeCoordinates, idBlockLast - 1 ) << " \n";

    return distanceMinPU;
}

void barrierAlpha()
{
    mpi::environment env;
    mpi::communicator world;
    for ( int r = 0; r < world.size(); ++r )
    {
        world.barrier();
        if ( r == world.rank() )
        {
            std::cout << "RRRRRRRRRRRRRRRRRRRRR ALPhA Rank [" << r << "] BARRIER RRRRRRRRRRRRRRRRRRRRR\n";
        }
    }
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
    int numRank = 0;
    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType2>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    //******************************************************************************************************************/
    // BUILD RAYS
    const double epsilon = 0.00001f;
    int nbValues = 0;

    std::vector<bvh_ray_type> rays;
    std::vector<double> distance_Real_mode;
    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;

    int nbRays = std::sqrt( number_rays_desired );
    nbRays = nbRays * nbRays;

    // size_t kblock = 10000;
    size_t kblock = nbRays + 1;

    if ( isViewInfo ) std::cout << "[INFO] Number of ray  = " << nbRays << "\n";
    if ( isViewInfo ) std::cout << "[INFO] Number of allNodeCoordinates   = " << allNodeCoordinates.size() << "\n";
    if ( isViewInfo ) std::cout << "[INFO] Max Size of allNodeCoordinates =" << allNodeCoordinates.max_size() << "\n";

    if ( isViewInfo ) std::cout << "[INFO] Generate Rays"
                                << "\n";
    size_t estimatedRayCount = allNodeCoordinates.size() * nbRays;
    if ( isViewInfo ) std::cout << "[INFO] Estimation Nb Rays =" << estimatedRayCount << "\n";
    rays = generateRays<bvh_ray_type>( allNodeCoordinates, number_rays_desired );
    for ( size_t k = 0; k < rays.size(); ++k )
    {
        raysDistributed.push_back( rays[k] );
    }
    if ( isViewInfo ) std::cout << "[INFO] Generate " << rays.size() << " Rays Done"
                                << "\n";
    //******************************************************************************************************************/

    barrierAlpha();

    //******************************************************************************************************************/
    //==================================================================================================================/
    //******************************************************************************************************************/

#if 1
    //******************************************************************************************************************/
    // BUILD BVH CPU
    // BEGIN::Build BVH CPU
    t_begin_cpu = std::chrono::steady_clock::now();
    tic();
    auto bvhThirdParty = boundingVolumeHierarchy( _range = range, _kind = "third-party", _quality = BVHEnum::Quality::High );
    auto timeBVHCPUDuration = toc( "timeBVHCPUDuration" );
    t_end_bvh_cpu = std::chrono::steady_clock::now();
    // END::Build BVH CPU
    numRank = bvhThirdParty->worldComm().rank();
    t_laps_CPU = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_cpu - t_begin_cpu ).count();
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH CPU : " << t_laps_CPU << " us"
                                << " rank=[" << numRank << "]\n";
    sleep( 1 );
    //******************************************************************************************************************/
#endif

#if 1
    //******************************************************************************************************************/
    // BUILD BVH GPU
    // BEGIN::Build BVH GPU
    t_begin_gpu = std::chrono::steady_clock::now();
    tic();
    auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-party" );
    // auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-multi-gpu-party" );
    auto timeBVHGPUDuration = toc( "timeBVHGPUDuration" );
    t_end_bvh_gpu = std::chrono::steady_clock::now();
    // END::Build BVH GPU
    numRank = bvhHIPParty->worldComm().rank();
    t_laps_GPU = std::chrono::duration_cast<std::chrono::microseconds>( t_end_bvh_gpu - t_begin_gpu ).count();
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH GPU : " << t_laps_GPU << " us"
                                << " rank=[" << numRank << "]\n";
    sleep( 1 );
    //******************************************************************************************************************/
#endif

    //******************************************************************************************************************/
    //==================================================================================================================/
    //******************************************************************************************************************/

    std::vector<double> dist;
    std::vector<double> distance_CPU_mode;
    std::vector<double> distance_GPU_mode;

    std::vector<size_t> id;
    std::vector<size_t> id_CPU;
    std::vector<size_t> id_GPU;

#if 1
    //******************************************************************************************************************/
    // BUILD CPU RAY TRACKING
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
    if ( isViewInfo ) std::cout << "HHHHHHHHHHHHHHHHHHH timeRTCPUDuration " << timeRTCPUDuration << " \n";
        //******************************************************************************************************************/
#endif

#if 1
    //******************************************************************************************************************/
    // BUILD GPU RAY TRACKING
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
    if ( isViewInfo ) std::cout << "HHHHHHHHHHHHHHHHHHH timeRTGPUDuration " << timeRTGPUDuration << " \n";
        //******************************************************************************************************************/
#endif

        //******************************************************************************************************************/
        //==================================================================================================================/
        //******************************************************************************************************************/

#if 1
    //******************************************************************************************************************/
    // CALCULATE DISTANCE MIN TO SURFACE CUBE CPU AND GPU
    if ( isViewInfo ) std::cout << "\n\n";
    if ( isViewInfo ) std::cout << "[INFO] Calul MinDist CPU\n";
    std::vector<double> distanceMinCPU;
    distanceMinCPU = calculateDistanceMinPU( id_CPU, distance_CPU_mode, kblock, allNodeCoordinates, true );
    if ( isViewInfo ) std::cout << "[INFO] Calul MinDist GPU\n";
    std::vector<double> distanceMinGPU;
    distanceMinGPU = calculateDistanceMinPU( id_GPU, distance_GPU_mode, kblock, allNodeCoordinates, true );
    //******************************************************************************************************************/
#endif

    //******************************************************************************************************************/
    //==================================================================================================================/
    //******************************************************************************************************************/

    //******************************************************************************************************************/
    // DEBRIEFING

    for ( size_t index = 0; index < allNodeCoordinates.size(); ++index )
    {
        double distanceMinREAL = calculateDist( allNodeCoordinates, index );
        double errCPU = abs( distanceMinCPU[index] - distanceMinREAL );
        double errGPU = abs( distanceMinGPU[index] - distanceMinREAL );

        DataDistanceErrTimeAll data = {
            numRank,
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

    allDataPU.rank = numRank;
    allDataPU.t_laps_BVH_CPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_bvh_cpu - t_begin_cpu ).count();
    allDataPU.t_laps_RT_CPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_raytracing_cpu - t_begin_raytracing_cpu ).count();
    allDataPU.t_laps_BVH_GPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_bvh_gpu - t_begin_gpu ).count();
    allDataPU.t_laps_RT_GPU = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_raytracing_gpu - t_begin_raytracing_gpu ).count();
    allDataPU.nbRays = nbRays;
    //******************************************************************************************************************/

    printf( "FINISHED\n" );
}

void saveAllData(
    const int maxNumElements,
    const int maxNumFaces,
    const int maxNumPoints,
    const int maxNumVertices,
    const std::vector<std::vector<double>>& allNodeCoordinates,
    const std::vector<DataDistanceErrTimeAll>& allDataDistanceBVHRTAll,
    const DataTimeLapsConfig& allDataPU )
{

    const std::vector<std::pair<std::string, std::function<void( std::ofstream& )>>> files = {
        { "all_results_per_vertex.csv", [&]( std::ofstream& file )
          {
              file << "Num Rank,Num Vertex,PosX,PosY,PosZ,distanceMinREAL,distanceFastMarching,errFastMarching,distanceMinCPU,errCPU,distanceMinGPU,errGPU\n";
              for ( size_t i = 0; i < allNodeCoordinates.size(); ++i )
              {
                  file << allDataDistanceBVHRTAll[i].rank << ","
                       << allDataDistanceBVHRTAll[i].id << ","
                       << std::fixed << std::setprecision( 9 )
                       << allNodeCoordinates[i][0] << ","
                       << allNodeCoordinates[i][1] << ","
                       << allNodeCoordinates[i][2] << ","
                       << allDataDistanceBVHRTAll[i].distanceMinREAL << ","
                       << allDataDistanceBVHRTAll[i].distanceFastMarching << ","
                       << allDataDistanceBVHRTAll[i].errFastMarching << ","
                       << allDataDistanceBVHRTAll[i].distanceMinCPU << ","
                       << allDataDistanceBVHRTAll[i].errCPU << ","
                       << allDataDistanceBVHRTAll[i].distanceMinGPU << ","
                       << allDataDistanceBVHRTAll[i].errGPU << "\n";
              }
          } },
        { "results.csv", [&]( std::ofstream& file )
          {
              file << "rank=" << allDataPU.rank << "\n"
                   << "hsize=" << allDataPU.hsize << "\n"
                   << "maxNumElement=" << maxNumElements << "\n"
                   << "maxNumFace=" << maxNumFaces << "\n"
                   << "maxNumPoints=" << maxNumPoints << "\n"
                   << "maxNumVerices=" << maxNumVertices << "\n"
                   << "nbRaysDesired=" << allDataPU.nbRaysDesired << "\n"
                   << "nbRays=" << allDataPU.nbRays << "\n"
                   << "timeBVHcpu=" << allDataPU.t_laps_BVH_CPU << "\n"
                   << "timeRTcpu=" << allDataPU.t_laps_RT_CPU << "\n"
                   << "timeBVHgpu=" << allDataPU.t_laps_BVH_GPU << "\n"
                   << "timeRTgpu=" << allDataPU.t_laps_RT_GPU << "\n"
                   << "timeFastMarching=" << allDataPU.t_laps_FastMarching << "\n"
                   << "totalTimeBVHRTcpu=" << allDataPU.t_laps_BVH_CPU + allDataPU.t_laps_RT_CPU << "\n"
                   << "totalTimeBVHRTgpu=" << allDataPU.t_laps_BVH_GPU + allDataPU.t_laps_RT_GPU << "\n";
          } },
        { "results2.csv", [&]( std::ofstream& file )
          {
              file << "rank,hsize,maxNumElement,maxNumFace,maxNumPoints,maxNumVerices,nbRaysDesired,nbRays,"
                   << "timeBVHcpu,timeRTcpu,timeBVHgpu,timeRTgpu,timeFastMarching,totalTimeBVHRTcpu,totalTimeBVHRTgpu\n"
                   << allDataPU.rank << ","
                   << allDataPU.hsize << ","
                   << maxNumElements << ","
                   << maxNumFaces << ","
                   << maxNumPoints << ","
                   << maxNumVertices << ","
                   << allDataPU.nbRaysDesired << ","
                   << allDataPU.nbRays << ","
                   << allDataPU.t_laps_BVH_CPU << ","
                   << allDataPU.t_laps_RT_CPU << ","
                   << allDataPU.t_laps_BVH_GPU << ","
                   << allDataPU.t_laps_RT_GPU << ","
                   << allDataPU.t_laps_FastMarching << ","
                   << allDataPU.t_laps_BVH_CPU + allDataPU.t_laps_RT_CPU << ","
                   << allDataPU.t_laps_BVH_GPU + allDataPU.t_laps_RT_GPU << "\n";
          } } };

    for ( const auto& [filename, writeFunc] : files )
    {
        std::ofstream file( filename );
        if ( !file.is_open() )
        {
            std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
            continue;
        }
        writeFunc( file );
        file.close();
    }
}

BOOST_AUTO_TEST_SUITE( distance_bvh_cpu_gpu_gpu_tests )

BOOST_AUTO_TEST_CASE( all_distance )
{
    mpi::environment env;
    mpi::communicator world;

    int numRank = world.rank();

    // signal(SIGTERM, sigterm_handler);

    // We read the value of "hsize" and "number_rays_desired"
    double hsize = option( _name = "hsize" ).as<double>();
    int number_rays_desired = option( _name = "number_rays_desired" ).as<int>();
    bool isViewInfo = option( _name = "isViewInfo" ).as<bool>();

    bool isPreheating = true; // isPreheating = false;
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

    std::vector<std::vector<double>> allNodeCoordinates;
    DataTimeLapsConfig allDataPU;
    allDataPU.nbRaysDesired = number_rays_desired;
    allDataPU.hsize = hsize;
    //std::vector<DataDistanceErrTime> allDataDistanceBVHRT;

    // List of node coordinates
    for ( size_type k = 0; k < Vh->nLocalDofWithGhost(); ++k )
    {
        auto const& dofPt = Vh->dof()->dofPoint( k ).template get<0>();
        allNodeCoordinates.push_back( { dofPt[0], dofPt[1], dofPt[2] } );
    }

    int nbNode = allNodeCoordinates.size();

    //******************************************************************************************************************/
    // Calculates Node points to Surface distances by the method FastMarching
    std::chrono::steady_clock::time_point t_begin_FastMarching, t_end_FastMarching;
    t_begin_FastMarching = std::chrono::steady_clock::now();
    auto distToBoundary = distanceToRange( _space = Vh, _range = submeshFaces );
    t_end_FastMarching = std::chrono::steady_clock::now();
    long int t_laps_FastMarching = std::chrono::duration_cast<std::chrono::milliseconds>( t_end_FastMarching - t_begin_FastMarching ).count();
    //******************************************************************************************************************/

    //******************************************************************************************************************/
    //==================================================================================================================/
    //******************************************************************************************************************/

    bool isOn = true; // isOn = false;

    if ( isOn )
    {
        // In this part all data is sent at once from CPU to GPU.
        std::vector<DataDistanceErrTimeAll> allDataDistanceBVHRTAll;
        // We calculate the distances from the edge of the cube and the intersection points.
        distToBoundaryBVHpuSendAllNode( rangeFaces, allDataPU, allNodeCoordinates, allDataDistanceBVHRTAll, true );

        // Add FastMarching Part
        allDataPU.t_laps_FastMarching = t_laps_FastMarching;
        for ( int i = 0; i < nbNode; ++i )
        {
            allDataDistanceBVHRTAll[i].distanceFastMarching = distToBoundary[i];
            allDataDistanceBVHRTAll[i].errFastMarching = abs( distToBoundary[i] - allDataDistanceBVHRTAll[i].distanceMinREAL );
        }

        //******************************************************************************************************************/
        // Save All Data
        saveAllData( mesh->maxNumElements(), mesh->maxNumFaces(), mesh->maxNumPoints(), mesh->maxNumVertices(),
                      allNodeCoordinates, allDataDistanceBVHRTAll, allDataPU );

        //******************************************************************************************************************/
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
        //******************************************************************************************************************/

        //******************************************************************************************************************/
        //==================================================================================================================/
        //******************************************************************************************************************/

        if ( isViewInfo )
        {
            std::cout << "\n";
            std::cout << "[INFO] Elapsed microseconds\n";
            std::cout << "[INFO] BVH CPU : " << allDataPU.t_laps_BVH_CPU << " ms\n";
            std::cout << "[INFO] RT  CPU : " << allDataPU.t_laps_RT_CPU << " ms\n";
            std::cout << "[INFO] BVH GPU : " << allDataPU.t_laps_BVH_GPU << " ms\n";
            std::cout << "[INFO] RT  GPU : " << allDataPU.t_laps_RT_GPU << " ms\n";
            std::cout << "[INFO] FastMarching : " << allDataPU.t_laps_FastMarching << " ms\n";
            std::cout << "\n";
        }
    }


}

BOOST_AUTO_TEST_SUITE_END()
