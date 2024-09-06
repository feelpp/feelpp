#define BOOST_TEST_MODULE bvhgpu2_tests


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
//#include <feel/feeltask/taskpu.hpp>


//#include <hip/hip_runtime.h>
//#include <hip/hip_runtime_api.h>

#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"
//#include "hipblas.h"
//#include "hipsolver.h"
//#include "hipblas-export.h"


#include <thrust/device_vector.h> 
#include <thrust/transform.h> 
#include <thrust/functional.h> 
#include <thrust/execution_policy.h>
#include <thrust/random.h>


#include <feel/feeltask/taskpu.hpp>

#include "lbvh.cuh"
#include <random>

using namespace Feel;

#define BLOCK_SIZE 128



//#include "Gu.hpp"
//#include "GuHip.hpp"


//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================


//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


struct saxpy_functor
{
    const float a;
    saxpy_functor(float _a) : a(_a) {}
    __host__ __device__
    float operator()(const float& x, const float& y) const { return a * x + y; }
};

struct aabb_getter
{
    __device__
    lbvh::aabb<float> operator()(const float4 f) const noexcept
    {
        lbvh::aabb<float> retval;
        retval.upper = f;
        retval.lower = f;
        return retval;
    }
};

struct distance_calculator
{
    __device__
    float operator()(const float4 point, const float4 object) const noexcept
    {
        return (point.x - object.x) * (point.x - object.x) +
               (point.y - object.y) * (point.y - object.y) +
               (point.z - object.z) * (point.z - object.z);
    }
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================





//=========================================================================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=========================================================================================================================





BOOST_AUTO_TEST_SUITE( bvh_intersection_gpu_thrust_tests )


BOOST_AUTO_TEST_CASE( test_gpu_thrust_with_specx2 )
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



BOOST_AUTO_TEST_CASE( test_gpu_thrust_saxpy )
{
    const int N = 10;
    thrust::device_vector<float> d_x(N, 1.0f); 
    thrust::device_vector<float> d_y(N, 2.0f);
    float a = 2.0f;

    for (int i = 0; i < d_x.size(); i++) { std::cout << d_x[i] << " "; }
    std::cout << std::endl;

    for (int i = 0; i < d_y.size(); i++) { std::cout << d_y[i] << " "; }
    std::cout << std::endl;

    thrust::transform(thrust::device,
        d_x.begin(), 
        d_x.end(), 
        d_y.begin(), 
        d_y.begin(), 
        saxpy_functor(a));

    for (int i = 0; i < d_y.size(); i++) { std::cout << d_y[i] << " "; }
    std::cout << std::endl;
}


BOOST_AUTO_TEST_CASE( test_gpu_thrust_lbvh )
{
    constexpr std::size_t N=10;
    std::vector<float4> ps(N);
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<float> uni(0.0, 1.0);

    for(auto& p : ps)
    {
        p.x = uni(mt);
        p.y = uni(mt);
        p.z = uni(mt);
    }

    lbvh::bvh<float, float4, aabb_getter> bvh(ps.begin(), ps.end(), true);

    const auto bvh_dev = bvh.get_device_repr();

    std::cout << "testing query_device:overlap ...\n";
    thrust::for_each(thrust::device,
        thrust::make_counting_iterator<std::size_t>(0),
        thrust::make_counting_iterator<std::size_t>(N),
        [bvh_dev] __device__ (std::size_t idx) {
            unsigned int buffer[10];
            const auto self = bvh_dev.objects[idx];
            const float  dr = 0.1f;
            for(std::size_t i=1; i<10; ++i)
            {
                for(unsigned int j=0; j<10; ++j)
                {
                    buffer[j] = 0xFFFFFFFF;
                }
                const float r = dr * i;
                lbvh::aabb<float> query_box;
                query_box.lower = make_float4(self.x-r, self.y-r, self.z-r, 0);
                query_box.upper = make_float4(self.x+r, self.y+r, self.z+r, 0);
                const auto num_found = lbvh::query_device(
                        bvh_dev, lbvh::overlaps(query_box), buffer, 10);

                for(unsigned int j=0; j<10; ++j)
                {
                    const auto jdx    = buffer[j];
                    if(j >= num_found)
                    {
                        assert(jdx == 0xFFFFFFFF);
                        continue;
                    }
                    else
                    {
                        assert(jdx != 0xFFFFFFFF);
                        assert(jdx < bvh_dev.num_objects);
                    }
                    const auto other  = bvh_dev.objects[jdx];
                    assert(fabsf(self.x - other.x) < r); // check coordinates
                    assert(fabsf(self.y - other.y) < r); // are in the range
                    assert(fabsf(self.z - other.z) < r); // of query box
                }
            }
            return ;
        });

    std::cout << "testing query_device:nearest_neighbor ...\n";
    thrust::for_each(thrust::device,
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(N),
        [bvh_dev] __device__ (const unsigned int idx) {
            const auto self = bvh_dev.objects[idx];
            const auto nest = lbvh::query_device(bvh_dev, lbvh::nearest(self),
                                                 distance_calculator());
            assert(nest.first != 0xFFFFFFFF);
            const auto other   = bvh_dev.objects[nest.first];
            // of course, the nearest object is itself.
            assert(nest.second == 0.0f);
            assert(self.x == other.x);
            assert(self.y == other.y);
            assert(self.z == other.z);
            return ;
       });

    thrust::device_vector<float4> random_points(N);
    thrust::transform(
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(N),
        random_points.begin(), [] __device__(const unsigned int idx) {
            thrust::default_random_engine rand;
            thrust::uniform_real_distribution<float> uni(0.0f, 1.0f);
            rand.discard(idx);
            const float x = uni(rand);
            const float y = uni(rand);
            const float z = uni(rand);
            return make_float4(x, y, z, 0);
        });

    thrust::for_each(random_points.begin(), random_points.end(),
        [bvh_dev] __device__ (const float4 pos) {
            const auto calc = distance_calculator();
            const auto nest = lbvh::query_device(bvh_dev, lbvh::nearest(pos), calc);
            assert(nest.first != 0xFFFFFFFF);

            for(unsigned int i=0; i<bvh_dev.num_objects; ++i)
            {
                const auto dist = calc(bvh_dev.objects[i], pos);
                //printf("%f\n",dist);
                if(i == nest.first)
                {
                    assert(dist == nest.second);
                }
                else
                {
                    assert(dist >= nest.second);
                }
            }
            return ;
        });

    std::cout << "testing query_host:overlap ...\n";
    {
        for(std::size_t i=0; i<10; ++i)
        {
            const auto self = bvh.objects_host()[i];
            const float dr = 0.1f;
            for(unsigned int cnt=1; cnt<10; ++cnt)
            {
                const float r = dr * cnt;
                lbvh::aabb<float> query_box;
                query_box.lower = make_float4(self.x-r, self.y-r, self.z-r, 0);
                query_box.upper = make_float4(self.x+r, self.y+r, self.z+r, 0);

                std::vector<std::size_t> buffer;
                const auto num_found = lbvh::query_host(bvh,
                        lbvh::overlaps(query_box), std::back_inserter(buffer));

                for(unsigned int jdx : buffer)
                {
                    assert(jdx < bvh.objects_host().size());

                    const auto other  = bvh.objects_host()[jdx];
                    assert(fabsf(self.x - other.x) < r); // check coordinates
                    assert(fabsf(self.y - other.y) < r); // are in the range
                    assert(fabsf(self.z - other.z) < r); // of query box
                }
            }
        }
    }
}



BOOST_AUTO_TEST_SUITE_END()






