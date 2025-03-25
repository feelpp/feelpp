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


#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
//#include <boost/math/statistics/shapiro_wilk.hpp>

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

    about.addAuthor( "Patrick Lemoine", "developer", "Noname@cemosis.fr", "" );
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
    double px;
    double py;
    double pz;
    double distanceMinREAL;
    double distanceFastMarching;
    double errFastMarching;
    double distanceMinCPU;
    double errCPU;
    double distanceMinGPU;
    double errGPU;

    // Add this part for serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        ar& rank;
        ar& id;
        ar& px;
        ar& py;
        ar& pz;
        ar& distanceMinREAL;
        ar& distanceFastMarching;
        ar& errFastMarching;
        ar& distanceMinCPU;
        ar& errCPU;
        ar& distanceMinGPU;
        ar& errGPU;
    }
};

struct DataTimeLapsConfig
{
    int rank;
    size_t nbRays;
    size_t nbRaysDesired;
    double hsize;
    int maxNumElements;
    int maxNumFaces;
    int maxNumPoints;
    int maxNumVertices;
    long int t_laps_BVH_CPU;
    long int t_laps_BVH_GPU;
    long int t_laps_RT_CPU;
    long int t_laps_RT_GPU;
    long int t_laps_FastMarching;

    // Add this part for serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        ar& rank;
        ar& nbRays;
        ar& nbRaysDesired;
        ar& hsize;
        ar& maxNumElements;
        ar& maxNumFaces;
        ar& maxNumPoints;
        ar& maxNumVertices;
        ar& t_laps_BVH_CPU;
        ar& t_laps_BVH_GPU;
        ar& t_laps_RT_CPU;
        ar& t_laps_RT_GPU;
        ar& t_laps_FastMarching;
    }
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

struct BoundingBoxMesh {
    Point min;
    Point max;
};


BoundingBoxMesh calculateBoundingBoxMesh(const Mesh<Simplex<3, 1, 3>>& mesh) {
    // This will be useful later ...
    BoundingBoxMesh bbox;
    bbox.min = Point(std::numeric_limits<double>::max());
    bbox.max = Point(std::numeric_limits<double>::lowest());
    if (mesh.maxNumPoints() > 0) {
        for (auto pointIndex = 0; pointIndex < mesh.maxNumPoints(); ++pointIndex) {
            auto point = mesh.point(pointIndex);
            bbox.min[0] = std::min(bbox.min[0], point[0]);
            bbox.min[1] = std::min(bbox.min[1], point[1]);
            bbox.min[2] = std::min(bbox.min[2], point[2]);
            bbox.max[0] = std::max(bbox.max[0], point[0]);
            bbox.max[1] = std::max(bbox.max[1], point[1]);
            bbox.max[2] = std::max(bbox.max[2], point[2]);
        }
    } else {
        std::cout << "No points in the mesh." << std::endl;
    }

    return bbox;
}

template <typename MeshEntityType>
struct MeshPrimitiveInfo
{
    // For more information and to establish a ray tracing strategy according to the box cpu meshs.
    using mesh_entity_type = std::decay_t<typename MeshEntityType::type>;
    static constexpr uint16_type nDim = mesh_entity_type::nDim;
    static constexpr uint16_type nRealDim = mesh_entity_type::nRealDim;
    using vector_realdim_type = Eigen::Matrix<double, nRealDim, 1>;

    MeshPrimitiveInfo(MeshEntityType const& meshEntity)
        : M_meshEntity(meshEntity.get()) 
    {
        auto verticesUblas = M_meshEntity.vertices();
        auto G = Feel::em_cmatrix_col_type<double>(verticesUblas.data().begin(), nRealDim, mesh_entity_type::numVertices);
        M_bound_min = G.rowwise().minCoeff();
        M_bound_max = G.rowwise().maxCoeff();
        M_bound_min.array() -= 2 * Feel::type_traits<double>::epsilon();
        M_bound_max.array() += 2 * Feel::type_traits<double>::epsilon();
        auto bary = M_meshEntity.barycenter();
        M_centroid = Eigen::Map<vector_realdim_type>(bary.data().begin(), nRealDim);
    }

    mesh_entity_type const& meshEntity() const { return M_meshEntity; }
    vector_realdim_type const& boundMin() const noexcept { return M_bound_min; }
    vector_realdim_type const& boundMax() const noexcept { return M_bound_max; }
    vector_realdim_type const& centroid() const noexcept { return M_centroid; }

private:
    vector_realdim_type M_bound_min;
    vector_realdim_type M_bound_max;
    vector_realdim_type M_centroid;
    mesh_entity_type const& M_meshEntity;
};


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

template <typename BVHRayType>
std::vector<BVHRayType> generateCameraRays(
    const Eigen::Vector3d& cameraPosition,
    const Eigen::Vector3d& cameraTarget,
    const Eigen::Vector3d& cameraUp,
    int imageWidth,
    int imageHeight,
    double fieldOfViewDegrees )
{
    std::vector<BVHRayType> rays;
    rays.reserve( imageWidth * imageHeight );
    double fieldOfViewRadians = fieldOfViewDegrees * M_PI / 180.0;

    Eigen::Vector3d forward = ( cameraTarget - cameraPosition ).normalized();
    Eigen::Vector3d right = forward.cross( cameraUp ).normalized();
    Eigen::Vector3d up = right.cross( forward );

    double aspectRatio = static_cast<double>( imageWidth ) / imageHeight;
    double halfFovTan = std::tan( fieldOfViewRadians / 2.0 );
    Eigen::Vector3d horizontal = right * ( 2.0 * halfFovTan * aspectRatio );
    Eigen::Vector3d vertical = up * ( 2.0 * halfFovTan );
    Eigen::Vector3d viewportUpperLeft = forward - horizontal / 2.0 + vertical / 2.0;

    int rayId = 0;
    for ( int y = 0; y < imageHeight; ++y )
    {
        for ( int x = 0; x < imageWidth; ++x )
        {
            double ndcX = ( 2.0 * ( x + 0.5 ) / imageWidth - 1.0 ) * aspectRatio;
            double ndcY = 1.0 - 2.0 * ( y + 0.5 ) / imageHeight;

            Eigen::Vector3d rayDirection = viewportUpperLeft + horizontal * ( ( x + 0.5 ) / imageWidth ) - vertical * ( ( y + 0.5 ) / imageHeight );
            rayDirection.normalize();
            rays.emplace_back( cameraPosition, rayDirection );
            rays.back().id = rayId++;
        }
    }

    return rays;
}

void savePPM( const std::string& filename, unsigned char* data, int width,
              int height )
{
    std::ofstream file( filename, std::ios::binary );
    file << "P6\n"
         << width << " " << height << "\n255\n";
    file.write( reinterpret_cast<char*>( data ), width * height * 3 );
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

void barrierAlpha( int numFlag )
{
    mpi::environment env;
    mpi::communicator world;
    for ( int r = 0; r < world.size(); ++r )
    {
        world.barrier();
        if ( r == world.rank() )
        {
            std::cout << "RRRRRRRRRRRRRRRRRRRRR ALPhA Rank [" << r << "] BARRIER " << world.size() << " numFlag [" << numFlag << "] RRRRRRRRRRRRRRRRRRRRR" << std::endl;
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

    barrierAlpha( 1 );

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
    //auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-multi-gpu-party" );
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
    bool isViewInfoLevel2 = false; // isViewInfoLevel2 = true;
    std::vector<double> distanceMinCPU;
    distanceMinCPU = calculateDistanceMinPU( id_CPU, distance_CPU_mode, kblock, allNodeCoordinates, isViewInfoLevel2 );
    if ( isViewInfo ) std::cout << "[INFO] Calul MinDist GPU\n";
    std::vector<double> distanceMinGPU;
    distanceMinGPU = calculateDistanceMinPU( id_GPU, distance_GPU_mode, kblock, allNodeCoordinates, isViewInfoLevel2 );
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
            allNodeCoordinates[index][0],
            allNodeCoordinates[index][1],
            allNodeCoordinates[index][2],
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

// ===============
// Debriefing part
// ===============

void saveAllData(
    const std::string nameFile,
    const std::vector<DataDistanceErrTimeAll>& allDataDistanceBVHRTAll,
    const DataTimeLapsConfig& allDataPU )

{
    const std::vector<std::pair<std::string, std::function<void( std::ofstream& )>>> files = {
        { nameFile + "_all_data.csv", [&]( std::ofstream& file )
          {
              file << "Num Rank,Num Vertex,PosX,PosY,PosZ,distanceMinREAL,distanceFastMarching,errFastMarching,distanceMinCPU,errCPU,distanceMinGPU,errGPU\n";
              for ( size_t i = 0; i < allDataDistanceBVHRTAll.size(); ++i )
              {
                  file << allDataDistanceBVHRTAll[i].rank << ","
                       << allDataDistanceBVHRTAll[i].id << ","
                       << std::fixed << std::setprecision( 9 )
                       << allDataDistanceBVHRTAll[i].px << ","
                       << allDataDistanceBVHRTAll[i].py << ","
                       << allDataDistanceBVHRTAll[i].pz << ","
                       << allDataDistanceBVHRTAll[i].distanceMinREAL << ","
                       << allDataDistanceBVHRTAll[i].distanceFastMarching << ","
                       << allDataDistanceBVHRTAll[i].errFastMarching << ","
                       << allDataDistanceBVHRTAll[i].distanceMinCPU << ","
                       << allDataDistanceBVHRTAll[i].errCPU << ","
                       << allDataDistanceBVHRTAll[i].distanceMinGPU << ","
                       << allDataDistanceBVHRTAll[i].errGPU << "\n";
              }
          } },
        { nameFile + "_in_column.csv", [&]( std::ofstream& file )
          {
              file << "rank=" << allDataPU.rank << "\n"
                   << "hsize=" << allDataPU.hsize << "\n"
                   << "maxNumElement=" << allDataPU.maxNumElements << "\n"
                   << "maxNumFace=" << allDataPU.maxNumFaces << "\n"
                   << "maxNumPoints=" << allDataPU.maxNumPoints << "\n"
                   << "maxNumVerices=" << allDataPU.maxNumVertices << "\n"
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
        { nameFile + "_in_line.csv", [&]( std::ofstream& file )
          {
              file << "rank,hsize,maxNumElement,maxNumFace,maxNumPoints,maxNumVerices,nbRaysDesired,nbRays,"
                   << "timeBVHcpu,timeRTcpu,timeBVHgpu,timeRTgpu,timeFastMarching,totalTimeBVHRTcpu,totalTimeBVHRTgpu\n"
                   << allDataPU.rank << ","
                   << allDataPU.hsize << ","
                   << allDataPU.maxNumElements << ","
                   << allDataPU.maxNumFaces << ","
                   << allDataPU.maxNumPoints << ","
                   << allDataPU.maxNumVertices << ","
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

    // Add other files if needed for debriefing ...

    for ( const auto& [filename, writeFunc] : files )
    {
        std::ofstream file( filename );
        if ( !file.is_open() )
        {
            std::cerr << "Error opening file " << filename << std::endl;
            continue;
        }
        writeFunc( file );
        file.close();
    }
}

void saveAllDataDebriefing(
    const std::string nameFile,
    const std::vector<DataDistanceErrTimeAll>& allDataDistanceBVHRTAll,
    const std::vector<DataTimeLapsConfig>& allDataPU )

{
    const std::vector<std::pair<std::string, std::function<void( std::ofstream& )>>> files = {
        { nameFile + "_all_data.csv", [&]( std::ofstream& file )
          {
              file << "Num Rank,Num Vertex,PosX,PosY,PosZ,distanceMinREAL,distanceFastMarching,errFastMarching,distanceMinCPU,errCPU,distanceMinGPU,errGPU\n";
              for ( size_t i = 0; i < allDataDistanceBVHRTAll.size(); ++i )
              {
                  file << allDataDistanceBVHRTAll[i].rank << ","
                       << allDataDistanceBVHRTAll[i].id << ","
                       << std::fixed << std::setprecision( 9 )
                       << allDataDistanceBVHRTAll[i].px << ","
                       << allDataDistanceBVHRTAll[i].py << ","
                       << allDataDistanceBVHRTAll[i].pz << ","
                       << allDataDistanceBVHRTAll[i].distanceMinREAL << ","
                       << allDataDistanceBVHRTAll[i].distanceFastMarching << ","
                       << allDataDistanceBVHRTAll[i].errFastMarching << ","
                       << allDataDistanceBVHRTAll[i].distanceMinCPU << ","
                       << allDataDistanceBVHRTAll[i].errCPU << ","
                       << allDataDistanceBVHRTAll[i].distanceMinGPU << ","
                       << allDataDistanceBVHRTAll[i].errGPU << "\n";
              }
          } },
        { nameFile + "_in_line.csv", [&]( std::ofstream& file )
          {
              file << "rank,hsize,maxNumElement,maxNumFace,maxNumPoints,maxNumVerices,nbRaysDesired,nbRays,"
                   << "timeBVHcpu,timeRTcpu,timeBVHgpu,timeRTgpu,timeFastMarching,totalTimeBVHRTcpu,totalTimeBVHRTgpu\n";
              for ( size_t i = 0; i < allDataPU.size(); ++i )
              {
                  file << allDataPU[i].rank << ","
                       << allDataPU[i].hsize << ","
                       << allDataPU[i].maxNumElements << ","
                       << allDataPU[i].maxNumFaces << ","
                       << allDataPU[i].maxNumPoints << ","
                       << allDataPU[i].maxNumVertices << ","
                       << allDataPU[i].nbRaysDesired << ","
                       << allDataPU[i].nbRays << ","
                       << allDataPU[i].t_laps_BVH_CPU << ","
                       << allDataPU[i].t_laps_RT_CPU << ","
                       << allDataPU[i].t_laps_BVH_GPU << ","
                       << allDataPU[i].t_laps_RT_GPU << ","
                       << allDataPU[i].t_laps_FastMarching << ","
                       << allDataPU[i].t_laps_BVH_CPU + allDataPU[i].t_laps_RT_CPU << ","
                       << allDataPU[i].t_laps_BVH_GPU + allDataPU[i].t_laps_RT_GPU << "\n";
              }
          } } };

    for ( const auto& [filename, writeFunc] : files )
    {
        std::ofstream file( filename );
        if ( !file.is_open() )
        {
            std::cerr << "Error opening file " << filename << std::endl;
            continue;
        }
        writeFunc( file );
        file.close();
    }
}

void saveAllDataDebriefingJSON(
    const std::string nameFile,
    const std::vector<DataDistanceErrTimeAll>& allDataDistanceBVHRTAll,
    const std::vector<DataTimeLapsConfig>& allDataPU )
{
    using json = nlohmann::json;

    json allData;

    json allDataDistanceJson = json::array();
    for ( const auto& data : allDataDistanceBVHRTAll )
    {
        allDataDistanceJson.push_back( { { "rank", data.rank },
                                         { "id", data.id },
                                         { "px", data.px },
                                         { "py", data.py },
                                         { "pz", data.pz },
                                         { "distanceMinREAL", data.distanceMinREAL },
                                         { "distanceFastMarching", data.distanceFastMarching },
                                         { "errFastMarching", data.errFastMarching },
                                         { "distanceMinCPU", data.distanceMinCPU },
                                         { "errCPU", data.errCPU },
                                         { "distanceMinGPU", data.distanceMinGPU },
                                         { "errGPU", data.errGPU } } );
    }
    allData["allDataDistanceBVHRTAll"] = allDataDistanceJson;

    json allDataPUJson = json::array();
    for ( const auto& data : allDataPU )
    {
        allDataPUJson.push_back( { { "rank", data.rank },
                                   { "hsize", data.hsize },
                                   { "maxNumElements", data.maxNumElements },
                                   { "maxNumFaces", data.maxNumFaces },
                                   { "maxNumPoints", data.maxNumPoints },
                                   { "maxNumVertices", data.maxNumVertices },
                                   { "nbRaysDesired", data.nbRaysDesired },
                                   { "nbRays", data.nbRays },
                                   { "t_laps_BVH_CPU", data.t_laps_BVH_CPU },
                                   { "t_laps_RT_CPU", data.t_laps_RT_CPU },
                                   { "t_laps_BVH_GPU", data.t_laps_BVH_GPU },
                                   { "t_laps_RT_GPU", data.t_laps_RT_GPU },
                                   { "t_laps_FastMarching", data.t_laps_FastMarching },
                                   { "totalTimeBVHRTcpu", data.t_laps_BVH_CPU + data.t_laps_RT_CPU },
                                   { "totalTimeBVHRTgpu", data.t_laps_BVH_GPU + data.t_laps_RT_GPU } } );
    }
    allData["allDataPU"] = allDataPUJson;

    std::ofstream file( nameFile + "_debriefing_data.json" );
    if ( !file.is_open() )
    {
        std::cerr << "Error opening file " << nameFile + "_debriefing_data.json" << std::endl;
        return;
    }
    file << std::setw( 4 ) << allData << std::endl;
    file.close();
}

// =================
// Meshs fusion part
// =================

template <typename MeshType>
std::shared_ptr<MeshType> concatenate( const std::shared_ptr<MeshType>& mesh1, const std::shared_ptr<MeshType>& mesh2 )
{
    auto result_mesh = std::make_shared<MeshType>( "concatenated_mesh", mesh1->worldCommPtr() );
    for ( auto const& elt : elements( mesh1 ) )
    {
        result_mesh->addElement( typename MeshType::element_type( elt ) );
    }
    for ( auto const& elt : elements( mesh2 ) )
    {
        result_mesh->addElement( typename MeshType::element_type( elt ) );
    }
    for ( auto const& pt : points( mesh1 ) )
    {
        result_mesh->addPoint( pt );
    }
    for ( auto const& pt : points( mesh2 ) )
    {
        result_mesh->addPoint( pt );
    }
    result_mesh->updateForUse();
    return result_mesh;
}

template <typename MeshType>
std::shared_ptr<MeshType> gatherMeshes( const std::shared_ptr<MeshType>& local_mesh )
{
    boost::mpi::communicator world;
    int rank = world.rank();
    int size = world.size();

    // Serialize the local mesh
    std::vector<char> local_mesh_data;
    {
        std::ostringstream oss;
        boost::archive::binary_oarchive oa( oss );
        oa << local_mesh;
        std::string str = oss.str();
        local_mesh_data.assign( str.begin(), str.end() );
    }

    // Collect mesh sizes
    std::vector<int> sizes( size );
    int local_size = static_cast<int>( local_mesh_data.size() );
    boost::mpi::gather( world, local_size, sizes, 0 );

    std::vector<char> received_data;
    if ( rank == 0 )
    {
        int total_size = std::accumulate( sizes.begin(), sizes.end(), 0 );
        received_data.resize( total_size );
    }

    boost::mpi::gatherv( world, local_mesh_data.data(), local_mesh_data.size(), received_data.data(), sizes, 0 );

    std::shared_ptr<MeshType> global_mesh;
    if ( rank == 0 )
    {
        std::vector<std::shared_ptr<MeshType>> meshes;
        int current_position = 0;
        for ( int i = 0; i < size; ++i )
        {
            std::vector<char> current_mesh_data( received_data.begin() + current_position,
                                                 received_data.begin() + current_position + sizes[i] );
            std::istringstream iss( std::string( current_mesh_data.begin(), current_mesh_data.end() ) );
            boost::archive::binary_iarchive ia( iss );
            std::shared_ptr<MeshType> mesh_part;
            ia >> mesh_part;
            meshes.push_back( mesh_part );
            current_position += sizes[i];
        }

        // We concatenate all meshes
        global_mesh = meshes[0];
        for ( size_t i = 1; i < meshes.size(); ++i )
        {
            global_mesh = concatenate( global_mesh, meshes[i] );
        }
    }

    return global_mesh;
}

// ==================
// Build picture part
// ==================

template <typename RangeType2>
void builtPicture(
    RangeType2 const& range )
{
    bool isViewInfo = true;
    int numRank = 0;
    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType2>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    std::chrono::steady_clock::time_point t_begin_BuildCameraRays, t_end_BuildCameraRays;
    std::chrono::steady_clock::time_point t_begin_BVH, t_end_BVH;
    std::chrono::steady_clock::time_point t_begin_RT, t_end_RT;
    std::chrono::steady_clock::time_point t_begin_BuildPicture, t_end_BuildPicture;
    std::chrono::steady_clock::time_point t_begin_AllProcess, t_end_AllProcess;

    long int t_laps_BuildCameraRays = 0;
    long int t_laps_BVH = 0;
    long int t_laps_RT = 0;
    long int t_laps_BuildPicture = 0;
    long int t_laps_AllProcess = 0;

    t_begin_AllProcess = std::chrono::steady_clock::now();
    //******************************************************************************************************************/
    // Build Rays Camera
    // const double epsilon = 0.00001f;
    // int nbValues = 0;

    std::vector<bvh_ray_type> rays;
    std::vector<double> distance_Real_mode;
    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;

    t_begin_BuildCameraRays = std::chrono::steady_clock::now();
    Eigen::Vector3d cameraPosition( 2.0, 2.0, 2.0 );
    Eigen::Vector3d cameraLookAt( 0.5, 0.5, 0.5 );
    Eigen::Vector3d cameraUp( 0.0, 1.0, 0.0 ); // Assuming Y is up

    int imageWidth = 800;
    int imageHeight = 800;
    double fieldOfViewDegrees = 60.0;
    rays = generateCameraRays<bvh_ray_type>( cameraPosition, cameraLookAt, cameraUp, imageWidth, imageHeight, fieldOfViewDegrees );
    for ( size_t k = 0; k < rays.size(); ++k )
    {
        raysDistributed.push_back( rays[k] );
    }
    t_end_BuildCameraRays = std::chrono::steady_clock::now();
    t_laps_BuildCameraRays = std::chrono::duration_cast<std::chrono::microseconds>( t_end_BuildCameraRays - t_begin_BuildCameraRays ).count();

    if ( isViewInfo ) std::cout << "[INFO] Generate " << rays.size() << " Rays Done"
                                << "\n";
    barrierAlpha( 0 );

    //******************************************************************************************************************/
    // Build BVH GPU
    t_begin_BVH = std::chrono::steady_clock::now();
    auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-party" );
    // auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-multi-gpu-party" );
    // BVH::CPU auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "third-party", _quality = BVHEnum::Quality::High );
    t_end_BVH = std::chrono::steady_clock::now();
    t_laps_BVH = std::chrono::duration_cast<std::chrono::microseconds>( t_end_BVH - t_begin_BVH ).count();
    numRank = bvhHIPParty->worldComm().rank();
    sleep( 1 );
    //******************************************************************************************************************/

    //******************************************************************************************************************/
    // Build GPU Rays Tracing
    std::vector<double> dist;
    std::vector<double> distance_GPU_mode;
    std::vector<size_t> id;
    std::vector<size_t> id_GPU;
    t_begin_RT = std::chrono::steady_clock::now();
    auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect( _ray = raysDistributed, _parallel = false );
    t_end_RT = std::chrono::steady_clock::now();
    t_laps_RT = std::chrono::duration_cast<std::chrono::microseconds>( t_end_RT - t_begin_RT ).count();

    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
    {
        dist = getAllDistanceRayIntersections( bvhHIPParty, rayIntersectionResult );
        distance_GPU_mode.insert( distance_GPU_mode.end(), dist.begin(), dist.end() );
        id = getId( bvhHIPParty, rayIntersectionResult );
        id_GPU.insert( id_GPU.end(), id.begin(), id.end() );
    }
    if ( isViewInfo ) std::cout << "[INFO] id_GPU " << id_GPU.size() << "\n";
    //******************************************************************************************************************/

    //******************************************************************************************************************/
    // Build Picture
    t_begin_BuildPicture = std::chrono::steady_clock::now();
    {
        // Normalize distances to the range [0, 1]
        double max_distance = 0.0;
        for ( const auto& d : distance_GPU_mode )
        {
            max_distance = std::max( max_distance, d );
        }

        double min_distance = max_distance;
        for ( const auto& d : distance_GPU_mode )
        {
            min_distance = std::min( min_distance, d );
        }

        // Create image data (RGB)
        unsigned char* image_data = new unsigned char[imageWidth * imageHeight * 3];
        for ( int i = 0; i < imageWidth * imageHeight; ++i )
        {
            unsigned char r, g, b;
            image_data[i * 3 + 0] = 0;
            image_data[i * 3 + 1] = 0;
            image_data[i * 3 + 2] = 0;
        }

        for ( int k = 0; k < distance_GPU_mode.size(); ++k )
        {
            unsigned char r, g, b;
            double distance = ( distance_GPU_mode[k] - min_distance ) / ( max_distance - min_distance );
            g = static_cast<unsigned char>( distance * 255.0 );
            r = static_cast<unsigned char>( ( 1.0 - distance ) * 255.0 );
            b = 0;
            int i = id_GPU[k];
            image_data[i * 3 + 0] = r;
            image_data[i * 3 + 1] = g;
            image_data[i * 3 + 2] = b;
        }

        savePPM( "distance_image.ppm", image_data, imageWidth, imageHeight );
        delete[] image_data;
    }
    t_end_BuildPicture = std::chrono::steady_clock::now();
    t_laps_BuildPicture = std::chrono::duration_cast<std::chrono::microseconds>( t_end_BuildPicture - t_begin_BuildPicture ).count();
    t_end_AllProcess = std::chrono::steady_clock::now();
    t_laps_AllProcess = std::chrono::duration_cast<std::chrono::microseconds>( t_end_AllProcess - t_begin_AllProcess ).count();

    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside Build Camera Rays : " << t_laps_BuildCameraRays << " us"
                                << " rank=[" << numRank << "]\n";
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside BVH : " << t_laps_BVH << " us"
                                << " rank=[" << numRank << "]\n";
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside RT : " << t_laps_RT << " us"
                                << " rank=[" << numRank << "]\n";
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside Build Picture : " << t_laps_BuildPicture << " us"
                                << " rank=[" << numRank << "]\n";
    if ( isViewInfo ) std::cout << "[INFO] Elapsed microseconds inside All Process : " << t_laps_AllProcess << " us"
                                << " rank=[" << numRank << "]\n";
}

// =========================
// Statistical analysis part
// =========================

double calculateMean( const std::vector<double>& vec )
{
    return std::accumulate( vec.begin(), vec.end(), 0.0 ) / vec.size();
}

double calculateStandardDeviation( const std::vector<double>& vec, double mean )
{
    double sumSquares = std::accumulate( vec.begin(), vec.end(), 0.0,
                                         [mean]( double acc, double val )
                                         { return acc + std::pow( val - mean, 2 ); } );
    return std::sqrt( sumSquares / ( vec.size() - 1 ) );
}

double calculateStandardDeviation( const std::vector<double>& values )
{
    double mean = std::accumulate( values.begin(), values.end(), 0.0 ) / values.size();
    double variance = std::accumulate( values.begin(), values.end(), 0.0,
                                       [mean]( double acc, double val )
                                       { return acc + std::pow( val - mean, 2 ); } );
    return std::sqrt( variance / ( values.size() - 1 ) );
}

std::pair<double, double> calculateConfidenceInterval( const std::vector<double>& data, double mean, double stdDev, double confidenceLevel = 0.95 )
{
    if ( data.empty() ) return std::make_pair( std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() );

    boost::math::students_t_distribution<> t_dist( data.size() - 1 );
    double t_value = boost::math::quantile( boost::math::complement( t_dist, ( 1.0 - confidenceLevel ) / 2.0 ) );
    double marginOfError = t_value * ( stdDev / std::sqrt( data.size() ) );

    return std::make_pair( mean - marginOfError, mean + marginOfError );
}

double calculateMeanRelativeError( const std::vector<double>& exact, const std::vector<double>& calculated )
{
    double errorSum = 0.0;
    int validCount = 0;
    for ( size_t i = 0; i < exact.size(); ++i )
    {
        if ( exact[i] != 0.0 )
        {
            errorSum += std::abs( exact[i] - calculated[i] ) / exact[i];
            validCount++;
        }
    }
    return validCount > 0 ? errorSum / validCount : 0.0;
}

std::pair<double, double> shapiro_wilk( std::vector<double>& x )
{
    int n = x.size();
    std::sort( x.begin(), x.end() );

    std::vector<double> a( n );
    double m = 0.0, s = 0.0;

    for ( int i = 0; i < n; i++ )
    {
        m += x[i];
    }
    m /= n;

    for ( int i = 0; i < n; i++ )
    {
        s += ( x[i] - m ) * ( x[i] - m );
    }

    if ( n == 3 )
    {
        a[0] = 0.7071;
        a[1] = 0;
    }
    else
    {
        double an = 0.5641896 * ( n - 3.0 ) / ( n + 1.0 );
        a[n - 1] = -an;
        a[0] = -a[n - 1];
        double alpha = ( a[0] - an ) * ( a[0] - an );
        double nn2 = n * n;
        for ( int i = 1; i < ( n - 1 ) / 2; i++ )
        {
            a[i] = a[0] - ( alpha / nn2 ) * ( nn2 - 1 - 2 * ( n - 1 - i ) * ( n - i ) );
            a[n - 1 - i] = -a[i];
        }
    }

    double w = 0.0;
    for ( int i = 0; i < n; i++ )
    {
        w += a[i] * x[i];
    }
    w = w * w / s;

    double mu = 0.0038915 * std::log( n ) * std::log( n ) * std::log( n ) - 0.083751 * std::log( n ) * std::log( n ) + 0.31082 * std::log( n ) - 1.5861;
    double sigma = std::exp( 0.0030302 * std::log( n ) * std::log( n ) - 0.082676 * std::log( n ) - 0.4803 );

    double z = ( std::log( 1 - w ) - mu ) / sigma;
    double p = 1 - 0.5 * ( 1 + std::erf( z / std::sqrt( 2 ) ) );

    return std::make_pair( w, p );
}

void chiSquaredTest( const std::vector<double>& exactDistances, const std::vector<double>& calculatedDistances )
{
    const int numCategories = 5;
    std::vector<int> exactCounts( numCategories, 0 );
    std::vector<int> calculatedCounts( numCategories, 0 );

    double minDistance = std::min( *std::min_element( exactDistances.begin(), exactDistances.end() ), *std::min_element( calculatedDistances.begin(), calculatedDistances.end() ) );
    double maxDistance = std::max( *std::max_element( exactDistances.begin(), exactDistances.end() ), *std::max_element( calculatedDistances.begin(), calculatedDistances.end() ) );

    double intervalSize = ( maxDistance - minDistance ) / numCategories;

    for ( double distance : exactDistances )
    {
        int category = std::min( static_cast<int>( ( distance - minDistance ) / intervalSize ), numCategories - 1 );
        exactCounts[category]++;
    }

    for ( double distance : calculatedDistances )
    {
        int category = std::min( static_cast<int>( ( distance - minDistance ) / intervalSize ), numCategories - 1 );
        calculatedCounts[category]++;
    }

    double chi2 = 0.0;
    for ( int i = 0; i < numCategories; ++i )
    {
        double expectedCount = ( exactCounts[i] + calculatedCounts[i] ) / 2.0;
        if ( expectedCount > 0 )
        {
            chi2 += ( ( exactCounts[i] - expectedCount ) * ( exactCounts[i] - expectedCount ) ) / expectedCount;
            chi2 += ( ( calculatedCounts[i] - expectedCount ) * ( calculatedCounts[i] - expectedCount ) ) / expectedCount;
        }
    }

    int degreesOfFreedom = numCategories - 1;

    boost::math::chi_squared_distribution<> chi2_dist( degreesOfFreedom );
    double p_value = 1 - boost::math::cdf( chi2_dist, chi2 );

    std::cout << "Chi-squared test:" << std::endl;
    std::cout << "  Chi-squared statistic: " << chi2 << std::endl;
    std::cout << "  Degrees of freedom: " << degreesOfFreedom << std::endl;
    std::cout << "  p-value: " << p_value << std::endl;

    if ( p_value < 0.05 )
    {
        std::cout << "  The distributions are likely different." << std::endl;
    }
    else
    {
        std::cout << "  There is not enough evidence to reject the hypothesis that the distributions are equal." << std::endl;
    }
}

void analyzeMethod( const std::vector<double>& exactDistances,
                    const std::vector<double>& calculatedDistances,
                    const std::string& methodName )
{
    std::vector<double> differences;
    for ( size_t i = 0; i < exactDistances.size(); ++i )
    {
        differences.push_back( calculatedDistances[i] - exactDistances[i] );
    }

    double meanDifference = calculateMean( differences );
    double stdDevDifference = calculateStandardDeviation( differences, meanDifference );
    double meanRelativeError = calculateMeanRelativeError( exactDistances, calculatedDistances );

    std::cout << "[INFO STAT]: Analysis for " << methodName << ":\n";
    std::cout << "  Mean of differences: " << meanDifference << std::endl;
    std::cout << "  Standard deviation of differences: " << stdDevDifference << std::endl;
    std::cout << "  Mean relative error: " << meanRelativeError << std::endl;

    // Confidence Interval
    auto confidenceInterval = calculateConfidenceInterval( differences, meanDifference, stdDevDifference );
    std::cout << "  Confidence Interval (95%): [" << confidenceInterval.first << ", " << confidenceInterval.second << "]" << std::endl;

    // Normality test (Shapiro-Wilk)
    auto result_shapiro = shapiro_wilk( differences );
    double w_shapiro = result_shapiro.first;
    double p_value_shapiro = result_shapiro.second;
    if ( p_value_shapiro < 0.05 )
    {
        std::cout << "  The data probably do not follow a normal distribution." << std::endl;
    }
    else
    {
        std::cout << "  There is not enough evidence to reject the normality of the data." << std::endl;
    }

    // Student's t-test
    double t_stat = meanDifference / ( stdDevDifference / std::sqrt( differences.size() ) );
    boost::math::students_t_distribution<> t_dist( differences.size() - 1 );
    double p_value = 2 * ( 1 - boost::math::cdf( t_dist, std::abs( t_stat ) ) );
    std::cout << "  Student's t-test: p-value = " << p_value << std::endl;

    // Chi-squared test
    chiSquaredTest( exactDistances, calculatedDistances );
}

void compareMethods( const std::vector<double>& exactDistances,
                     const std::vector<double>& calculatedDistances1,
                     const std::vector<double>& calculatedDistances2 )
{

    std::cout << std::endl;

    if ( exactDistances.size() != calculatedDistances1.size() ||
         exactDistances.size() != calculatedDistances2.size() )
    {
        std::cerr << "[INFO STAT]: Error: Distance vectors do not have the same size." << std::endl;
        return;
    }

    std::cout << "[INFO STAT]: Comparison of the two methods:\n\n";

    analyzeMethod( exactDistances, calculatedDistances1, "Method 1" );
    std::cout << std::endl;
    analyzeMethod( exactDistances, calculatedDistances2, "Method 2" );

    // Direct comparison of errors from both methods
    std::vector<double> errorDifferences;
    for ( size_t i = 0; i < exactDistances.size(); ++i )
    {
        double error1 = std::abs( calculatedDistances1[i] - exactDistances[i] );
        double error2 = std::abs( calculatedDistances2[i] - exactDistances[i] );
        errorDifferences.push_back( error1 - error2 );
    }

    double meanErrorDifference = calculateMean( errorDifferences );
    double stdDevErrorDifference = calculateStandardDeviation( errorDifferences, meanErrorDifference );

    std::cout << "\n[INFO STAT]: Direct comparison of errors (Method 1 - Method 2):\n";
    std::cout << "  Mean of error differences: " << meanErrorDifference << std::endl;
    std::cout << "  Standard deviation of error differences: " << stdDevErrorDifference << std::endl;

    // Paired t-test
    double t_stat = meanErrorDifference / ( stdDevErrorDifference / std::sqrt( errorDifferences.size() ) );
    boost::math::students_t_distribution<> t_dist( errorDifferences.size() - 1 );
    double p_value = 2 * ( 1 - boost::math::cdf( t_dist, std::abs( t_stat ) ) );
    std::cout << "  Paired t-test: p-value = " << p_value << std::endl;

    std::cout << "\n[INFO STAT]: Conclusion: ";
    if ( p_value < 0.05 )
    {
        if ( meanErrorDifference < 0 )
        {
            std::cout << "Method 2 is statistically better than Method 1." << std::endl;
        }
        else
        {
            std::cout << "Method 1 is statistically better than Method 2." << std::endl;
        }
    }
    else
    {
        std::cout << "There is no statistically significant difference between the two methods." << std::endl;
    }
}

// Monte Carlo Part

std::pair<double, double> calculateConfidenceInterval( const std::vector<double>& data, double confidenceLevel = 0.95 )
{
    if ( data.empty() ) return std::make_pair( std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() );

    double mean = calculateMean( data );
    double stdDev = calculateStandardDeviation( data );
    boost::math::students_t_distribution<> t_dist( data.size() - 1 );
    double t_value = boost::math::quantile( boost::math::complement( t_dist, ( 1.0 - confidenceLevel ) / 2.0 ) );
    double marginOfError = t_value * ( stdDev / std::sqrt( data.size() ) );

    return std::make_pair( mean - marginOfError, mean + marginOfError );
}

double calculateConvergenceRate( const std::vector<int>& numSamples, const std::vector<double>& errors )
{
    std::vector<double> logSamples, logErrors;
    for ( size_t i = 0; i < numSamples.size(); ++i )
    {
        logSamples.push_back( std::log( numSamples[i] ) );
        logErrors.push_back( std::log( errors[i] ) );
    }

    double n = logSamples.size();
    double sumX = std::accumulate( logSamples.begin(), logSamples.end(), 0.0 );
    double sumY = std::accumulate( logErrors.begin(), logErrors.end(), 0.0 );
    double sumXY = std::inner_product( logSamples.begin(), logSamples.end(), logErrors.begin(), 0.0 );
    double sumX2 = std::inner_product( logSamples.begin(), logSamples.end(), logSamples.begin(), 0.0 );

    double slope = ( n * sumXY - sumX * sumY ) / ( n * sumX2 - sumX * sumX );
    return -slope; // <=== The rate of convergence is the opposite of the slope
}

template <typename RangeType2>
void analyzeMonteCarloMethodConvergenceRate(
    RangeType2 const& range,
    std::vector<std::vector<double>>& allNodeCoordinates,
    DataTimeLapsConfig& allDataPU,
    const std::vector<int>& numRaysList )
{
    std::vector<double> convergenceErrorsCPU;
    std::vector<double> convergenceErrorsGPU;

    for ( const auto& nbRay : numRaysList )
    {
        std::vector<DataDistanceErrTimeAll> allDataDistanceBVHRTAlltmp;
        allDataPU.nbRaysDesired = nbRay;
        distToBoundaryBVHpuSendAllNode( range, allDataPU, allNodeCoordinates, allDataDistanceBVHRTAlltmp, true );

        std::vector<double> distancesMinREAL, distancesMinCPU, distancesMinGPU;
        for ( const auto& data : allDataDistanceBVHRTAlltmp )
        {
            distancesMinREAL.push_back( data.distanceMinREAL );
            distancesMinCPU.push_back( data.distanceMinCPU );
            distancesMinGPU.push_back( data.distanceMinGPU );
        }

        convergenceErrorsCPU.push_back( calculateMeanRelativeError( distancesMinREAL, distancesMinCPU ) );
        convergenceErrorsGPU.push_back( calculateMeanRelativeError( distancesMinREAL, distancesMinGPU ) );
    }

    double convergenceRateCPU = calculateConvergenceRate( numRaysList, convergenceErrorsCPU );
    double convergenceRateGPU = calculateConvergenceRate( numRaysList, convergenceErrorsGPU );

    std::cout << "Estimated convergence rate CPU: " << convergenceRateCPU << std::endl;
    std::cout << "Estimated convergence rate GPU: " << convergenceRateGPU << std::endl;
    std::cout << "Theoretical convergence rate for Monte Carlo: 0.5" << std::endl;

    auto ciCPU = calculateConfidenceInterval( convergenceErrorsCPU );
    auto ciGPU = calculateConfidenceInterval( convergenceErrorsGPU );

    std::cout << "95% Confidence Interval for CPU errors: [" << ciCPU.first << ", " << ciCPU.second << "]" << std::endl;
    std::cout << "95% Confidence Interval for GPU errors: [" << ciGPU.first << ", " << ciGPU.second << "]" << std::endl;
}

// ==========================
// Read TicToc Time file part
// ==========================

struct Data
{
    double count;
    double total;
    double max;
    double min;
    double mean;
    double stddev;
};

std::map<std::string, Data> loadData( const std::string& file )
{
    std::map<std::string, Data> data;
    std::ifstream fileStream( file );
    if ( !fileStream )
    {
        std::cerr << "Error opening file." << std::endl;
        return data;
    }

    std::string line;
    // Skip the first two lines (header)
    std::getline( fileStream, line );
    std::getline( fileStream, line );

    while ( std::getline( fileStream, line ) )
    {
        std::istringstream iss( line );
        std::string name;
        Data dataItem;
        std::getline( iss, name, '|' );
        std::getline( iss, line, '|' );
        std::getline( iss, line, '|' );
        iss >> dataItem.count;
        iss.ignore();
        iss >> dataItem.total;
        iss.ignore();
        iss >> dataItem.max;
        iss.ignore();
        iss >> dataItem.min;
        iss.ignore();
        iss >> dataItem.mean;
        iss.ignore();
        iss >> dataItem.stddev;
        data[name] = dataItem;
    }
    fileStream.close();
    return data;
}

void query( const std::map<std::string, Data>& data, const std::string& name )
{
    if ( data.find( name ) != data.end() )
    {
        const Data& dataItem = data.at( name );
        std::cout << "Name: " << name << std::endl;
        std::cout << "Count: " << dataItem.count << std::endl;
        std::cout << "Total: " << dataItem.total << std::endl;
        std::cout << "Max: " << dataItem.max << std::endl;
        std::cout << "Min: " << dataItem.min << std::endl;
        std::cout << "Mean: " << dataItem.mean << std::endl;
        std::cout << "StdDev: " << dataItem.stddev << std::endl;
    }
    else
    {
        std::cout << "Name not found." << std::endl;
    }
}






BOOST_AUTO_TEST_SUITE( distance_bvh_cpu_gpu_gpu_tests )

BOOST_AUTO_TEST_CASE( all_distance )
{
    mpi::environment env;
    mpi::communicator world;

    int numRank = world.rank();
    int nbCPUs = world.size();

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


    // Just to CTRL if it works for instance... stored for MPI CPU boxes
    // après on fera du lancé de rayons dessus pour optiniser et résuire les temps
    /*
    if ( isViewInfo )
    {
        Eigen::Matrix<double, 3, 1> global_min = Eigen::Matrix<double, 3, 1>::Constant(std::numeric_limits<double>::max());
        Eigen::Matrix<double, 3, 1> global_max = Eigen::Matrix<double, 3, 1>::Constant(std::numeric_limits<double>::lowest());

        for (auto const& element : elements(mesh))
        {
            MeshPrimitiveInfo<std::decay_t<decltype(element)>> element_info(element);
            global_min = global_min.cwiseMin(element_info.boundMin());
            global_max = global_max.cwiseMax(element_info.boundMax());
        }

        std::cout << "[INFO] Bounding Box Min: " << global_min.transpose() << std::endl;
        std::cout << "[INFO] Bounding Box Max: " << global_max.transpose() << std::endl;
    }
    */


    // Selecting what you want to process
    auto rangeFaces = markedfaces( mesh );
    auto submeshFaces = boundaryfaces( mesh );
    auto rangeElements = markedelements( mesh );

    auto Vh = Pch<1>( mesh );

    std::vector<std::vector<double>> allNodeCoordinates;
    DataTimeLapsConfig allDataPU;
    allDataPU.nbRaysDesired = number_rays_desired;
    allDataPU.hsize = hsize;
    allDataPU.maxNumElements = mesh->maxNumElements();
    allDataPU.maxNumFaces = mesh->maxNumFaces();
    allDataPU.maxNumPoints = mesh->maxNumPoints();
    allDataPU.maxNumVertices = mesh->maxNumVertices();

    // std::vector<DataDistanceErrTime> allDataDistanceBVHRT;

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
    bool isBuildPictureOn = true;
    isBuildPictureOn = false;
    bool isStatisticalAnalysisOn = true;
    //isStatisticalAnalysisOn = false; // add limit inf
    bool isSaveTicTocTime = true;    // isSaveTicTocTime  = false;
    bool isStatisticalAnalysisConvergenceRateMonteCarloOn = true;
    isStatisticalAnalysisConvergenceRateMonteCarloOn = false;

    // In this part all data is sent at once from CPU to GPU.
    std::vector<DataDistanceErrTimeAll> allDataDistanceBVHRTAll;

    if ( isOn )
    {
        // 3D rendering
        if ( isBuildPictureOn ) builtPicture( rangeFaces );
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
        // Statistical analysis
        if ( (isStatisticalAnalysisOn)  && (number_rays_desired>16) )
        {
            std::vector<double> distancesMinREAL;
            std::vector<double> distancesMinCPU;
            std::vector<double> distancesMinGPU;

            for ( const auto& data : allDataDistanceBVHRTAll )
            {
                // std::cout<<"data.distanceMinREAL="<<data.distanceMinREAL<<" data.distanceMinCPU="<<data.distanceMinCPU<<" data.distanceMinGPU="<<data.distanceMinGPU<<std::endl;
                distancesMinREAL.push_back( data.distanceMinREAL );
                distancesMinCPU.push_back( data.distanceMinCPU );
                distancesMinGPU.push_back( data.distanceMinGPU );
            }

            // analyzeMethod(distancesMinREAL, distancesMinCPU, "BVH CPU comparison");
            // analyzeMethod(distancesMinREAL, distancesMinGPU, "BVH GPU comparison");
            compareMethods( distancesMinREAL, distancesMinCPU, distancesMinGPU );
        }
        //******************************************************************************************************************/

        //******************************************************************************************************************/
        // Statistical analysis Convergence Rate Monte Carlo
        if (isStatisticalAnalysisConvergenceRateMonteCarloOn)
        {
            std::vector<int> numRaysList = { 100, 1000, 10000, 100000 };
            analyzeMonteCarloMethodConvergenceRate( rangeFaces, allNodeCoordinates, allDataPU, numRaysList );
        }
        //******************************************************************************************************************/

        //******************************************************************************************************************/
        // Save All tic toc time information
        if ( isSaveTicTocTime )
        {
            std::ofstream os( "tictoc.md" );
            Environment::saveTimersMD( os );
        }
        //******************************************************************************************************************/

        //******************************************************************************************************************/
        /*
            Todo: define the parameters that I will recover to put it in the debriefing.
            std::map<std::string, Data> data = loadData("tictoc.md");
                std::string name;
                std::cout << "Enter the query name: ";
                std::getline(std::cin, name);
                query(data, name);
        */
        //******************************************************************************************************************/

        //******************************************************************************************************************/
        // Save All Data for rank n
        saveAllData( "rank_" + std::to_string( numRank ) + "_results", allDataDistanceBVHRTAll, allDataPU );

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
            std::cout << "[INFO] Elapsed microseconds for Rank : " << numRank << "\n";
            std::cout << "[INFO] BVH CPU : " << allDataPU.t_laps_BVH_CPU << " ms\n";
            std::cout << "[INFO] RT  CPU : " << allDataPU.t_laps_RT_CPU << " ms\n";
            std::cout << "[INFO] BVH GPU : " << allDataPU.t_laps_BVH_GPU << " ms\n";
            std::cout << "[INFO] RT  GPU : " << allDataPU.t_laps_RT_GPU << " ms\n";
            std::cout << "[INFO] FastMarching : " << allDataPU.t_laps_FastMarching << " ms\n";
            std::cout << "\n";
        }
    }

    //******************************************************************************************************************/
    //==================================================================================================================/
    //******************************************************************************************************************/

    //******************************************************************************************************************/
    // Gathering results from all MPI CPUs. Then synthesis of the results...  Will see if it works properly ;-)

    barrierAlpha( 2 );

    // We collect the sizes of the local vectors.
    std::vector<int> sizes( world.size() );
    int local_size = allDataDistanceBVHRTAll.size();
    mpi::gather( world, local_size, sizes, 0 );

    // We prepare the vector to receive all the data on rank 0.
    std::vector<DataDistanceErrTimeAll> gatheredData;
    std::vector<DataTimeLapsConfig> gatheredDataTimeLaps;

    if ( world.rank() == 0 )
    {
        // we calculate the displacement.
        std::vector<int> displacements( world.size(), 0 );
        for ( int i = 1; i < world.size(); ++i )
        {
            displacements[i] = displacements[i - 1] + sizes[i - 1];
        }

        // we resize the reception vector.
        int total_size = std::accumulate( sizes.begin(), sizes.end(), 0 );
        gatheredData.resize( total_size );

        // we gather data from all ranks.
        mpi::gatherv( world, allDataDistanceBVHRTAll.data(), allDataDistanceBVHRTAll.size(),
                      gatheredData.data(), sizes, displacements, 0 );

        // we gather data from allDataPU from all ranks.
        gatheredDataTimeLaps.resize( world.size() );
        mpi::gather( world, allDataPU, gatheredDataTimeLaps, 0 );

        //... Save all data
        saveAllDataDebriefing( "debriefing_results", gatheredData, gatheredDataTimeLaps );
        saveAllDataDebriefingJSON( "debriefing_results", gatheredData, gatheredDataTimeLaps );
        // saveAllDataDebriefing( "debriefing_results",gatheredData);
    }
    else
    {
        // Other ranks simply send their data.
        mpi::gatherv( world, allDataDistanceBVHRTAll.data(), allDataDistanceBVHRTAll.size(), 0 );
        mpi::gather( world, allDataPU, 0 );
    }

    //******************************************************************************************************************/
    // Gathering Meshs results... Will see if it works properly ;-)

#if 0
    std::shared_ptr<decltype(mesh)::element_type> global_mesh;

    if (nbCPUs > 1) {
        global_mesh = gatherMeshes(mesh);
    } else {
        global_mesh = mesh;
    }

    // The goal is to check if the mesh is correct.
    if (numRank == 0) {
        std::cout << "Well done mesh assembled on master process.\n";
        std::cout << "Number of elements in the global mesh : " << global_mesh->numElements() << std::endl;
        auto e = exporter(_mesh=global_mesh, _name="my_global_mesh");
    }
#endif
}

BOOST_AUTO_TEST_SUITE_END()
