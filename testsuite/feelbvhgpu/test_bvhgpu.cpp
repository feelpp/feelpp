#define BOOST_TEST_MODULE bvhgpu_tests

// NOTA : Objective: Ray tracing from inside a cube using BVH Ray Tracing with a CPU and a GPU method. Compare the performances of the two methods.

#include <ranges>
#include <fmt/chrono.h>
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



#if defined(FEELPP_HAS_HIP)
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#endif

using namespace Feel;

// Section to use if you want to make meshes, otherwise it doesn't work.
inline AboutData
makeAbout()
{
    AboutData about( "test_bvhgpu",
                     "test_bvhgpu",
                     "0.2",
                     "nD(n=2,3) test bvh",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2024 Feel++ Consortium" );

    about.addAuthor( "Noname", "developer", "Noname@cemosis.fr", "" );
    return about;
}

inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()( "mesh2D.filename", po::value<std::string>(), "mesh2D.filename" )( "mesh3D.filename", po::value<std::string>(), "mesh3D.filename" );
    return opts;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

template <typename BvhType, typename RayIntersectionResultType>
void printRayIntersectionResults( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs, std::vector<typename std::remove_pointer_t<BvhType>::vector_realdim_type /*Eigen::Vector3d*/> const& pointIntersection )
{
    if ( bvh->worldComm().isMasterRank() )
        BOOST_TEST_MESSAGE( "Number of intersection: " << rirs.size() );
    BOOST_CHECK_MESSAGE( pointIntersection.size() == rirs.size(), fmt::format( "Number of intersection between ray and BVH tree is not correct : {} vs {}", pointIntersection.size(), rirs.size() ) );

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
                if constexpr ( std::decay_t<decltype( *bvh )>::nRealDim == 3 )
                    BOOST_TEST_MESSAGE( " --  Coordinates: " << rir.coordinates()[0] << "," << rir.coordinates()[1] << "," << rir.coordinates()[2] );
                BOOST_CHECK_SMALL( ( rir.coordinates() - pointIntersection[counter] ).norm(), 1e-8 );
            }

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            BOOST_TEST_MESSAGE( " --  Mesh entity id: " << prim.meshEntity().id() );
            BOOST_TEST_MESSAGE( " --  Mesh entity barycenter: " << prim.meshEntity().barycenter() );
            BOOST_TEST_MESSAGE( " ---------------------------------" );
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );
        counter++;
    }
}

template <typename BvhType, typename RayIntersectionResultType>
void printRayIntersectionResults2( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs, std::vector<typename std::remove_pointer_t<BvhType>::vector_realdim_type /*Eigen::Vector3d*/> const& pointIntersection )
{
    if ( bvh->worldComm().isMasterRank() )
        LOG( INFO ) << "Number of intersection: " << rirs.size() << "\n";

    int counter = 0;
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            LOG( INFO ) << " --  ProcessId: " << rir.processId() << "\n";
            LOG( INFO ) << " --  PrimitiveId: " << rir.primitiveId() << "\n";
            LOG( INFO ) << " --  Distance: " << rir.distance() << "\n";

            if ( rir.hasCoordinates() )
            {
                if constexpr ( std::decay_t<decltype( *bvh )>::nRealDim == 3 )
                    LOG( INFO ) << " --  Coordinates: " << rir.coordinates()[0] << "," << rir.coordinates()[1] << "," << rir.coordinates()[2] << "\n";
                ;
            }

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            LOG( INFO ) << " --  Mesh entity id: " << prim.meshEntity().id() << "\n";
            LOG( INFO ) << " --  Mesh entity barycenter: " << prim.meshEntity().barycenter() << "\n";
            LOG( INFO ) << " ---------------------------------"
                        << "\n";
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );
        counter++;
    }
}

template <typename BvhType, typename RayIntersectionResultType>
void printRayIntersectionResults( BvhType const& bvh, std::vector<RayIntersectionResultType> const& rirs )
{
    if ( bvh->worldComm().isMasterRank() )
        LOG( INFO ) << "Number of intersection: " << rirs.size() << "\n";
    for ( auto const& rir : rirs )
    {
        if ( rir.processId() == bvh->worldComm().rank() )
        {
            LOG( INFO ) << " --  ProcessId: " << rir.processId() << "\n";
            LOG( INFO ) << " --  PrimitiveId: " << rir.primitiveId() << "\n";
            LOG( INFO ) << " --  Distance: " << rir.distance() << "\n";

            if ( rir.hasCoordinates() )
            {
                if constexpr ( std::decay_t<decltype( *bvh )>::nRealDim == 3 )
                    LOG( INFO ) << " --  Coordinates: " << rir.coordinates()[0] << "," << rir.coordinates()[1] << "," << rir.coordinates()[2] << "\n";
                ;
            }

            auto const& prim = bvh->primitiveInfo( rir.primitiveId() );
            LOG( INFO ) << " --  Mesh entity id: " << prim.meshEntity().id() << "\n";
            LOG( INFO ) << " --  Mesh entity barycenter: " << prim.meshEntity().barycenter() << "\n";
            // LOG(INFO) << " --  Mesh entity barycenter: " << prim.meshEntity()
            LOG( INFO ) << " ---------------------------------"
                        << "\n";
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );
    }
}

template <typename RangeType>
void test3DWithHybrid( RangeType const& range )
{
    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    Eigen::Vector3d origin1 = { 1000.0, 0.0, 0.0 };
    Eigen::Vector3d direction_perp_1 = { 1., 0., 0. };

    Eigen::Vector3d origin2 = { -10.0, -0.25, -0.25 };
    Eigen::Vector3d direction_perp_2 = { 1., 0., 0. };

    Eigen::Vector3d origin3 = { -5.0, -0.20, -0.20 };
    Eigen::Vector3d direction_perp_3 = { 1., 0., 0. };

    std::vector<bvh_ray_type> rays;

    rays.push_back( bvh_ray_type( origin1, direction_perp_1 ) );
    // rays.push_back( bvh_ray_type(origin2,direction_perp_2) );
    // rays.push_back( bvh_ray_type(origin3,direction_perp_3) );

    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;

    for ( int k = 0; k < rays.size(); ++k )
    {
        raysDistributed.push_back( rays[k] );
    }

    LOG( INFO ) << "\n";

    // In normal CPU mode
    LOG( INFO ) << "In normal CPU mode\n";
    auto bvhThirdPartyLow = boundingVolumeHierarchy( _range = range, _kind = "third-party" );
    auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect( _ray = raysDistributed );
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
    {
        LOG( INFO ) << "multiRayDistributedIntersectionResult parts"
                    << "\n";
        BOOST_CHECK_MESSAGE( !rayIntersectionResult.empty(), fmt::format( "Intersection between ray and BVH tree has been found" ) );
        LOG( INFO ) << "Intersection between ray and BVH tree has been found"
                    << "\n";
        printRayIntersectionResults( bvhThirdPartyLow, rayIntersectionResult );
    }

    // In GPU mode with AMD HIP
    LOG( INFO ) << "In GPU mode with AMD HIP\n";
    auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-party" );
    auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect( _ray = raysDistributed );
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
    {
        LOG( INFO ) << "Hip parts"
                    << "\n";
        BOOST_CHECK_MESSAGE( !rayIntersectionResult.empty(), fmt::format( "Intersection between ray and BVH tree has been found" ) );
        LOG( INFO ) << "Intersection between ray and BVH tree has been found"
                    << "\n";
        printRayIntersectionResults( bvhHIPParty, rayIntersectionResult );
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
            // LOG(INFO) << " --  Distance: " << rir.distance()<< "\n";
            distance.push_back( rir.distance() );
        }
        bvh->worldComm().barrier();
        std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );
    }
    return distance;
}

Eigen::Vector3d sphericalToCartesian( double r, double theta, double alpha )
{
    Eigen::Vector3d position;
    position[0] = r * sin( theta ) * cos( alpha );
    position[1] = r * sin( theta ) * sin( alpha );
    position[2] = r * cos( theta );
    return position;
}

template <typename ExecSpace, typename RangeType>
void test3DInsideObjectWithHybrid( RangeType const& range )
{
    std::chrono::steady_clock::time_point t_begin_cpu, t_begin_gpu;
    std::chrono::steady_clock::time_point t_end_cpu, t_end_gpu;

    std::chrono::steady_clock::time_point t_begin_raytracing_cpu, t_begin_raytracing_gpu;
    std::chrono::steady_clock::time_point t_end_raytracing_cpu, t_end_raytracing_gpu;
    std::chrono::steady_clock::time_point t_end_bvh_cpu, t_end_bvh_gpu;

    long int t_laps;

    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    // Save all informations
    std::string filename = "results.txt";
    std::ofstream myfile( filename );

    for ( int kkk : std::views::iota(1, 5) )
    {

        Eigen::Vector3d ray_origin = { 0.0f, 0.0f, 0.0f };
        double thetaStart = 0.0f;
        double thetaEnd = M_PI;
        double alphaStart = 0.0f;
        double alphaEnd = 2.0f * M_PI;
        double thetaStep = M_PI / ( 18.0f * 0.125f * float( kkk ) );
        double alphaStep = M_PI / ( 18.0f * 0.125f * float( kkk ) );
        bool isViewInfo = false;

        std::vector<bvh_ray_type> rays;
        for ( double theta = thetaStart; theta <= thetaEnd; theta += thetaStep )
        {
            for ( double alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep )
            {
                Eigen::Vector3d ray_direction = sphericalToCartesian( 1.0f, theta, alpha );

                if ( isViewInfo )
                {
                    LOG( INFO ) << "Origin: <" << ray_origin[0] << ", " << ray_origin[1] << ", " << ray_origin[2] << ">"
                                << " ";
                    LOG( INFO ) << "Direction: <" << ray_direction[0] << ", " << ray_direction[1] << ", " << ray_direction[2] << ">" << std::endl;
                }
                rays.push_back( bvh_ray_type( ray_origin, ray_direction ) );
            }
        }

        BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
        for ( int k = 0; k < rays.size(); ++k )
        {
            raysDistributed.push_back( rays[k] );
        }

        std::vector<double> dist;

        if constexpr ( std::is_same_v<ExecSpace, Kokkos::Serial> )
        {
            // In normal CPU mode
            LOG( INFO ) << "In normal CPU mode\n";
            Kokkos::Timer timer;

            auto bvhThirdPartyLow = boundingVolumeHierarchy( _range = range, _kind = "third-party" );
            double time_bvh = timer.seconds();

            auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect( _ray = raysDistributed );
            double time_raytracing = timer.seconds();

            std::vector<double> distance_CPU_mode;
            for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
            {
                dist = getAllDistanceRayIntersections( bvhThirdPartyLow, rayIntersectionResult );
                // printRayIntersectionResults(bvhThirdPartyLow,rayIntersectionResult);
                distance_CPU_mode.insert( distance_CPU_mode.end(), dist.begin(), dist.end() );
            }
            double time_end = timer.seconds();
            LOG( INFO ) << fmt::format("[cpu] bvh : {}s, rt: {}s, total: {}s, distance GPU={}", time_bvh, time_raytracing, time_end, distance_CPU_mode.size() );
        }
#if defined(FEELPP_HAS_HIP)
        if constexpr ( std::is_same_v<ExecSpace, Kokkos::HIP> )
        {
            // In GPU mode with AMD HIP
            LOG( INFO ) << "In GPU mode with AMD HIP\n";

            Kokkos::Timer timer;
            auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-party" );
            double time_bvh = timer.seconds();timer.reset();

            auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect( _ray = raysDistributed );
            double time_raytracing = timer.seconds();timer.reset();

            std::vector<double> distance_GPU_mode;
            for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
            {
                dist = getAllDistanceRayIntersections( bvhHIPParty, rayIntersectionResult );
                // printRayIntersectionResults(bvhHIPParty,rayIntersectionResult);
                distance_GPU_mode.insert( distance_GPU_mode.end(), dist.begin(), dist.end() );
            }
            double time_end = timer.seconds();

            LOG( INFO ) << fmt::format("[gpu] bvh : {}s, rt: {}s, total: {}s, distance GPU={}", time_bvh, time_raytracing, time_end, distance_GPU_mode.size() );
        }
#endif
    }
#if 0
        // Distance comparison
        double deltaError = 0.00001f;
        double sumErrors = 0.0f;
        bool isError = false;

        if ( distance_GPU_mode.size() != distance_CPU_mode.size() )
        {
            isError = true;
            LOG( INFO ) << "Error size vector distance CPU vs GPU\n";
        }

        for ( int k; k < distance_GPU_mode.size(); ++k )
        {
            double e = fabs( distance_GPU_mode[k] - distance_CPU_mode[k] );
            sumErrors = sumErrors + e;
            if ( e > deltaError )
            {
                LOG( ERROR ) << k << " Distance CPU = " << distance_CPU_mode[k] << " Dist GPU = " << distance_GPU_mode[k] << " e=" << e << "\n";
                isError = true;
            }
        }
        if ( !isError )
        {
            LOG( INFO ) << "WELL DONE :-) No error (same distance). \n";
        }

#endif         
}

template <typename RangeType>
void test3D_AutoDecisionCPUorGPU( RangeType const& range )
{
    std::chrono::steady_clock::time_point t_begin;
    std::chrono::steady_clock::time_point t_end;

    std::chrono::steady_clock::time_point t_begin_raytracing;
    std::chrono::steady_clock::time_point t_end_raytracing;
    std::chrono::steady_clock::time_point t_end_bvh;

    long int t_laps;

    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHRay<mesh_entity_type::nRealDim>;

    bool isModeGPU = true;//isThereAnyGPUhere( false );

    int kkk = 1;

    Eigen::Vector3d ray_origin = { 0.0f, 0.0f, 0.0f };
    double thetaStart = 0.0f;
    double thetaEnd = M_PI;
    double alphaStart = 0.0f;
    double alphaEnd = 2.0f * M_PI;
    double thetaStep = M_PI / ( 18.0f * 0.125f * float( kkk ) );
    double alphaStep = M_PI / ( 18.0f * 0.125f * float( kkk ) );
    bool isViewInfo = false;

    std::vector<bvh_ray_type> rays;
    for ( double theta = thetaStart; theta <= thetaEnd; theta += thetaStep )
    {
        for ( double alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep )
        {
            Eigen::Vector3d ray_direction = sphericalToCartesian( 1.0f, theta, alpha );

            if ( isViewInfo )
            {
                LOG( INFO ) << "Origin: <" << ray_origin[0] << ", " << ray_origin[1] << ", " << ray_origin[2] << ">"
                            << " ";
                LOG( INFO ) << "Direction: <" << ray_direction[0] << ", " << ray_direction[1] << ", " << ray_direction[2] << ">" << std::endl;
            }
            rays.push_back( bvh_ray_type( ray_origin, ray_direction ) );
        }
    }

    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
    for ( int k = 0; k < rays.size(); ++k )
    {
        raysDistributed.push_back( rays[k] );
    }

    std::vector<double> dist;

    t_begin = std::chrono::steady_clock::now();

    if ( !isModeGPU )
    {
        LOG( INFO ) << "Run in normal CPU mode\n";
        auto bvhThirdPartyLow = boundingVolumeHierarchy( _range = range, _kind = "third-party" );
        auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect( _ray = raysDistributed );

        std::vector<double> distance_CPU_mode;
        for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            dist = getAllDistanceRayIntersections( bvhThirdPartyLow, rayIntersectionResult );
            // printRayIntersectionResults(bvhThirdPartyLow,rayIntersectionResult);
            distance_CPU_mode.insert( distance_CPU_mode.end(), dist.begin(), dist.end() );
        }
    }
#if defined(FEELPP_HAS_HIP)
    if ( isModeGPU )
    {
        // In GPU mode with AMD HIP
        LOG( INFO ) << "Run in GPU mode with AMD HIP\n";
        auto bvhHIPParty = boundingVolumeHierarchy( _range = range, _kind = "hip-party" );
        auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect( _ray = raysDistributed );

        std::vector<double> distance_GPU_mode;
        for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
        {
            dist = getAllDistanceRayIntersections( bvhHIPParty, rayIntersectionResult );
            // printRayIntersectionResults(bvhHIPParty,rayIntersectionResult);
            distance_GPU_mode.insert( distance_GPU_mode.end(), dist.begin(), dist.end() );
        }
    }
#endif    
    t_end = std::chrono::steady_clock::now();

    t_laps = std::chrono::duration_cast<std::chrono::microseconds>( t_end - t_begin ).count();
    LOG( INFO ) << "Elapsed microseconds inside BVH Ray Tracing GPU : " << t_laps << " us\n";
}

BOOST_AUTO_TEST_SUITE( bvh_intersection_gpu_tests )

BOOST_AUTO_TEST_CASE( test_load_mesh3 )
{

    // TEST Thrust
    thrust::device_vector<float> d_x( 10, 1.0f );
    using namespace Feel;
    using Feel::cout;
    // typedef Mesh<Simplex<2> > mesh_type;

    using mesh_type = Mesh<Simplex<3, 1, 3>>; //<Dim,Order,RDim>

    auto mesh = loadMesh( _mesh = new mesh_type );

    LOG( INFO ) << "[INFO] maxNumElement : " << mesh->maxNumElements() << std::endl;
    LOG( INFO ) << "[INFO] maxNumFace    : " << mesh->maxNumFaces() << std::endl;
    LOG( INFO ) << "[INFO] maxNumPoints  : " << mesh->maxNumPoints() << std::endl;
    LOG( INFO ) << "[INFO] maxNumVerices : " << mesh->maxNumVertices() << std::endl;

    auto rangeFaces = markedfaces( mesh ); //,{"CavityBottom","CavitySides","CavityTop","Up3","Down3","Back3","Left3","Rigth3"});
    auto submesh = createSubmesh( _mesh = mesh, _range = rangeFaces );

    auto Xhd0 = Pdh<0>( mesh );
    auto measures = Xhd0->element();
    measures.on( _range = elements( mesh ), _expr = vf::meas() );
    double measMin = measures.min();
    double measMax = measures.max();
    size_type nbdyfaces = nelements( boundaryfaces( mesh ) );

    LOG( INFO ) << "mesh entities" << std::endl;
    LOG( INFO ) << "   number of elements : " << mesh->numGlobalElements() << std::endl;
    LOG( INFO ) << "   number of faces : " << mesh->numGlobalFaces() << std::endl;
    LOG( INFO ) << "   number of boundary faces : " << nbdyfaces << std::endl;
    LOG( INFO ) << "   number of points : " << mesh->numGlobalPoints() << std::endl;
    LOG( INFO ) << "   number of vertices : " << mesh->numGlobalVertices() << std::endl;
    LOG( INFO ) << "mesh sizes" << std::endl;
    LOG( INFO ) << "   h max : " << mesh->hMax() << std::endl;
    LOG( INFO ) << "   h min : " << mesh->hMin() << std::endl;
    LOG( INFO ) << "   h avg : " << mesh->hAverage() << std::endl;
    LOG( INFO ) << "   measure : " << mesh->measure() << "\t" << measMin << " : " << measMax << std::endl;
    LOG( INFO ) << "Number of Partitions : " << mesh->numberOfPartitions() << std::endl;
    LOG( INFO ) << "nbdyfaces : " << nbdyfaces << "\n";
    LOG( INFO ) << "\n";

    
    LOG( INFO ) << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    LOG( INFO ) << "Execute on CPU" << std::endl;
    test3DInsideObjectWithHybrid<Kokkos::Serial>( rangeFaces );
    LOG( INFO ) << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";

#if defined(FEELPP_HAS_HIP)
    LOG( INFO ) << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    LOG( INFO ) << "Execute on GPU" << std::endl;
    test3DInsideObjectWithHybrid<Kokkos::HIP>( rangeFaces );
    LOG( INFO ) << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
#endif

    LOG( INFO ) << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
