
//#define COMPILE_WITH_CUDA
//#define COMPILE_WITH_HIP
//#define COMPILE_WITH_CXX_20

//#define USE_MPI
//#define USE_OpenMP
//#define USE_GPU_HIP  // <**** temporary
//#define UseCUDA
//#define UseHIP



#include <specx/Data/SpDataAccessMode.hpp>
#include <specx/Legacy/SpRuntime.hpp>
#include <specx/Task/SpPriority.hpp>
#include <specx/Task/SpProbability.hpp>
#include <specx/Utils/SpArrayView.hpp>
#include <specx/Utils/SpTimer.hpp>
#include <specx/Utils/small_vector.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpHeapBuffer.hpp>
#include <specx/Utils/SpUtils.hpp>
#include <specx/Utils/SpConsumerThread.hpp>
#include <specx/Legacy/SpRuntime.hpp>


#include <napp/na.hpp>



#ifdef USE_MPI
    #include <mpi.h>
#endif

#ifdef USE_OpenMP
	#include <omp.h>
#endif

#ifdef COMPILE_WITH_HIP
    #include "hip/hip_runtime.h"
    #include "hip/hip_runtime_api.h"
#endif

#ifdef COMPILE_WITH_CUDA
    #include "cuda_runtime.h"
    #include "cuda.h"
#endif



// Some macro functions for the AMD HIP GPU
#define HIP_KERNEL_NAME(...) __VA_ARGS__

#define HIP_CHECK(command) {               \
  hipError_t status = command;             \
  if (status!=hipSuccess) {                \
    std::cerr <<"Error: HIP reports "<< hipGetErrorString(status)<< std::endl; \
    std::abort(); } }


#ifdef NDEBUG
    #define HIP_ASSERT(x) x
#else
    #define HIP_ASSERT(x) (assert((x)==hipSuccess))
#endif


#define CUDA_CHECK(command) {               \
  cudaError_t status = command;             \
  if (status!=hipSuccess) {                \
    std::cerr <<"Error: CUDA reports "<< cudaGetErrorString(status)<< std::endl; \
    std::abort(); } }


#ifdef NDEBUG
    #define CUDA_ASSERT(x) x
#else
    #define CUDA_ASSERT(x) (assert((x)==cudaSuccess))
#endif


#define _param(...) _parameters=Sbtask::parameters(__VA_ARGS__)


//Links internal
#include "hardware.hpp"
#include "taskgpu.hpp"
#include "taskcpu.hpp"

//================================================================================================================================
// THE END.
//================================================================================================================================

