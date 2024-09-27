#define BOOST_TEST_MODULE bvhgpu_tests


/*
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/bvh.hpp>

#include <feel/feeldiscr/pdh.hpp>
*/

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



//#define COMPILE_WITH_HIP


#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"
//#include "hipblas.h"
//#include "hipsolver.h"
//#include "hipblas-export.h"


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


#include <feel/feeltask/taskpu.hpp>

using namespace Feel;

#define BLOCK_SIZE 128


#include "bvhHybrid.hpp"




__global__ void kernel(float* x,float* y,int n){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    //int idx = hipblockIdx.x * hipblockDim.x + hipthreadIdx.x;
    if (idx < n) y[idx] += 1;
}


void __device__ vector_add_device(const double vecA, const double vecB, double &vecC)
{
    vecC = vecA + vecB;
}

void __global__ vector_add(const double *vecA, const double *vecB, double *vecC, const int nb)
{
    const int i = hipBlockDim_x*hipBlockIdx_x+hipThreadIdx_x;

    if (i < nb)
        vector_add_device(vecA[i], vecB[i], vecC[i]); 
}


struct saxpy_functor
{
    const float a;
    saxpy_functor(float _a) : a(_a) {}
    __host__ __device__
    float operator()(const float& x, const float& y) const { return a * x + y; }
};

void check_solution(double coeff,double* a_in,double* b_in,double* c_in,int nb)
{
  printf("[INFO]: ");
	int errors = 0;
  	for (int i = 0; i < nb; i++) {
	    if (c_in[i] != coeff*(a_in[i] + b_in[i])) { errors++; }
	}
  	if (errors!=0) { printf("FAILED: %d errors\n",errors); } else { printf ("WELL DONE PASSED! :-)\n"); }
}

void write_vector(std::string ch,double* v,int nb)
{
  std::cout<<"[INFO]: "<<ch<<"> ";
	for (int i = 0; i < nb; i++) { std::cout<<int(v[i]); }
  std::cout<<std::endl;
}



//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================



//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
 
 

template <typename RangeType>
void test3D( RangeType const& range )
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

    //rays.push_back( bvh_ray_type(origin1,direction_perp_1) );
    rays.push_back( bvh_ray_type(origin2,direction_perp_2) );
    rays.push_back( bvh_ray_type(origin3,direction_perp_3) );

    BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;
    std::vector<std::size_t> raysDistributedIndices;
    for (int k=0;k<rays.size();++k)
    {
        raysDistributed.push_back( rays[k] );
        raysDistributedIndices.push_back( k );
    }

    std::vector<std::vector<Eigen::Vector3d>> pointIntersections;
    pointIntersections.push_back( { Eigen::Vector3d( 0, 0, 0 ) } ); //SEE
    pointIntersections.push_back( { Eigen::Vector3d( 0, 0, 0 ) } );
    pointIntersections.push_back( { Eigen::Vector3d( 0, 0, 0 ) } );


    std::cout<<"=====================================================================================+============\n";
    std::cout<<"++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++\n";
    std::cout<<"[INFO]: Bounding Colume Hierarchy\n";

    auto bvhInHouse = boundingVolumeHierarchy(_range=range,_kind="in-house");
    //auto bvhThirdParty = boundingVolumeHierarchy(_range=range,_kind="third-party");
    auto bvhThirdPartyLow = boundingVolumeHierarchy(_range=range,_kind="third-party",_quality=BVHEnum::Quality::High);

    //auto bvhGpuParty = boundingVolumeHierarchy(_range=range,_kind="gpu-party",_quality=BVHEnum::Quality::High);
    std::size_t counter = 0;
    
    for ( auto const& ray : rays )
        {
            std::cout << "rays parts"<< "\n";
            auto rayIntersectionResult1 = bvhInHouse->intersect(_ray=ray) ;
            BOOST_CHECK_MESSAGE(!rayIntersectionResult1.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhInHouse,rayIntersectionResult1 );

            //printRayIntersectionResults2(bvhInHouse,rayIntersectionResult1, pointIntersections[counter++] );
        }


    std::cout<<"*************************************************************************************************\n";

    auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            std::cout << "multiRayDistributedIntersectionResult parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhThirdPartyLow,rayIntersectionResult);
        }
    std::cout<<"=================================================================================================\n";

/*
    std::cout<<"*************************************************************************************************\n";

    auto multiRayDistributedIntersectionGpuResult = bvhGpuParty->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionGpuResult )
        {
            std::cout << "GPU parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhGpuParty,rayIntersectionResult);
        }

    std::cout<<"=================================================================================================\n";
*/


}




template <typename RangeType>
void test3DWithHybrid( RangeType const& range )
{

    using mesh_entity_type = std::remove_const_t<entity_range_t<RangeType>>;
    using bvh_ray_type = BVHHybrid::BVHRay<mesh_entity_type::nRealDim>;

    Eigen::Vector3d origin1={1000.0,0.0,0.0};
    Eigen::Vector3d direction_perp_1={1.,0.,0.};

    Eigen::Vector3d origin2={-10.0,-0.25,-0.25};
    Eigen::Vector3d direction_perp_2={1.,0.,0.};

    Eigen::Vector3d origin3={-5.0,-0.20,-0.20};
    Eigen::Vector3d direction_perp_3={1.,0.,0.};


    std::vector<bvh_ray_type> rays;

    rays.push_back( bvh_ray_type(origin1,direction_perp_1) );
    rays.push_back( bvh_ray_type(origin2,direction_perp_2) );
    rays.push_back( bvh_ray_type(origin3,direction_perp_3) );

    BVHHybrid::BVHRaysDistributed<mesh_entity_type::nRealDim> raysDistributed;

    for (int k=0;k<rays.size();++k)
    {
        raysDistributed.push_back( rays[k] );
    }

    std::cout<<"\n";
    std::cout<<"\n";
    std::cout<<"[INFO]: Bounding Colume Hierarchy GPU\n";
    auto bvhHIPParty       = BVHHybrid::boundingVolumeHierarchy(_range=range,_kind="hip-party");

    auto multiRayDistributedIntersectionHipResult = bvhHIPParty->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionHipResult )
        {
            std::cout << "Hip parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhHIPParty,rayIntersectionResult);
        }



    auto bvhThirdPartyLow = BVHHybrid::boundingVolumeHierarchy(_range=range,_kind="third-party");

    auto multiRayDistributedIntersectionResult = bvhThirdPartyLow->intersect(_ray=raysDistributed);
    for ( auto const& rayIntersectionResult : multiRayDistributedIntersectionResult )
        {
            std::cout << "multiRayDistributedIntersectionResult parts"<< "\n";
            BOOST_CHECK_MESSAGE(!rayIntersectionResult.empty(), fmt::format("Intersection between ray and BVH tree has been found"));
            std::cout << "Intersection between ray and BVH tree has been found"<< "\n";
            printRayIntersectionResults(bvhThirdPartyLow,rayIntersectionResult);
        }



}




BOOST_AUTO_TEST_SUITE( bvh_intersection_gpu_tests )

/*
BOOST_AUTO_TEST_CASE( test_gpu_with_specx )
{
    using namespace Feel;
    std::cout<<std::endl;
    int numDevices=0;
    HIP_CHECK(hipGetDeviceCount(&numDevices));
    std::cout<<"[INFO]: Get numDevice = "<<numDevices<<"\n";
    int deviceID=0;
    HIP_CHECK(hipGetDevice(&deviceID));
    std::cout<<"[INFO]: Get deviceID activated = "<<deviceID<<"\n";

    int nbThreads=3;
    int nbElements=10;
    double c0=2.0; 
    double c1=1.5; 
    int size_array = sizeof(double) * nbElements;

    double *h_vecA = (double *)malloc(size_array);
    double *h_vecB = (double *)malloc(size_array);
    double *h_vecC = (double *)malloc(size_array);

    for (int i = 0; i < nbElements; i++)
    {
        h_vecA[i] = 1;
        h_vecB[i] = 2;
        h_vecC[i] = 0;
    }

    write_vector("Vector A",h_vecA,nbElements);
    write_vector("Vector B",h_vecB,nbElements);


    auto FC1=[size_array,nbElements](const double* vh_vecA,const double* vh_vecB,double* vh_vecC0) {  
          double *d_vecA,*d_vecB,*d_vecC;
          hipMalloc((void **)&d_vecA, size_array);
          hipMalloc((void **)&d_vecB, size_array);
          hipMalloc((void **)&d_vecC, size_array);
          hipMemcpy(d_vecA, vh_vecA, size_array, hipMemcpyHostToDevice);
          hipMemcpy(d_vecB, vh_vecB, size_array, hipMemcpyHostToDevice);
          int grid_size = (nbElements + BLOCK_SIZE - 1) / BLOCK_SIZE;
          hipLaunchKernelGGL(vector_add, grid_size,BLOCK_SIZE, 0, 0, d_vecA, d_vecB, d_vecC, nbElements);
          hipMemcpy(vh_vecC0, d_vecC, size_array, hipMemcpyDeviceToHost);
          hipFree(d_vecA); hipFree(d_vecB); hipFree(d_vecC);
        return true;
    };
    
    int numMode=3;
    //NOTA: normal configuration
        if (numMode==1) {
            FC1(h_vecA,h_vecB,h_vecC);
        }

    //NOTA: we use normal specx without gpu
        if (numMode==2) {
            SpRuntime runtime(nbThreads);
            runtime.task(SpRead(h_vecA),SpRead(h_vecB),SpWrite(h_vecC),FC1);
            runtime.waitAllTasks();
            runtime.stopAllThreads();
            runtime.generateDot("Test_hip_AMD.dot", true);
            runtime.generateTrace("Test_hip_AMD.svg");  
        }


    //NOTA: we use my task class
        if (numMode==3) {
            int numTypeThread=3;  //0:No thread  1:Std::thread 2:std:async 3:Specx  
            Task::Task TsK( nbThreads, numTypeThread );
            TsK.setSave( false );
            TsK.setInfo( false );
            TsK.setFileName( "./Test_hip_AMD" );
            TsK.add( _param(h_vecA,h_vecB,h_vecC), _tasks = FC1 );
            TsK.run();
            TsK.close();
        }



    write_vector("Vector C = (A + B)",h_vecC,nbElements);
    check_solution(1.0,h_vecA,h_vecB,h_vecC,nbElements);
    free(h_vecA); free(h_vecB); free(h_vecC);
    std::cout<<"\n";
}
*/

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
    test3D( rangeFaces );
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    //test3D( elements(submesh) );
    std::cout<<"+*-*+*-*+*-*+*+*-*+*-*+*-*+*-*+*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n";
    test3DWithHybrid( rangeFaces );
    std::cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n";
    std::cout<<"\n";
}

/*
BOOST_AUTO_TEST_CASE( test_thrust )
{
    // Taille des vecteurs 
    const int N = 100;
    // Allocation et initialisation des vecteurs sur le device
    thrust::device_vector<float> d_x(N, 1.0f); 
    thrust::device_vector<float> d_y(N, 2.0f);
    // Scalaire a 
    float a = 2.0f;

    for (int i = 0; i < d_x.size(); i++) { std::cout << d_x[i] << " "; }
    std::cout << std::endl;

    for (int i = 0; i < d_y.size(); i++) { std::cout << d_y[i] << " "; }
    std::cout << std::endl;

    // Exécution de l'opération SAXPY
    thrust::transform(thrust::device,
        d_x.begin(), 
        d_x.end(), 
        d_y.begin(), 
        d_y.begin(), 
        saxpy_functor(a));

    // Affichage des résultats
    for (int i = 0; i < d_y.size(); i++) { std::cout << d_y[i] << " "; }
    std::cout << std::endl;

    std::cout << "[INFO]: WELL DONE :-) Test_thrust!"<<"\n";
}
*/



BOOST_AUTO_TEST_SUITE_END()




