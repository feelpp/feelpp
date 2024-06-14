#include <assert.h>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <iostream>


//Links for dev

#include <thread>
#include <vector>
#include <array>
#include <typeinfo>
#include <iostream>
#include <mutex>
#include <sched.h>
#include <pthread.h>

#include <algorithm> 
#include <string>
#include <utility>
#include <functional>
#include <future>
#include <cassert>
#include <chrono>
#include <type_traits>
#include <list>
#include <ranges>


//Links Specx

/*
#include "SpDataAccessMode.hpp"
#include "Utils/SpUtils.hpp"
#include "Task/SpTask.hpp"
#include "Legacy/SpRuntime.hpp"
#include "Utils/SpTimer.hpp"
#include "Utils/small_vector.hpp"
#include "Utils/SpConsumerThread.hpp"
#include "SpComputeEngine.hpp"
#include "Speculation/SpSpeculativeModel.hpp"
*/

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


//Links Eigen
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>




//#define COMPILE_WITH_CUDA
//#define COMPILE_WITH_HIP
//#define COMPILE_WITH_CXX_20
//#define USE_GPU_HIP  // <**** temporary

//Links mpi
//#define USE_MPI
#ifdef USE_MPI
    #include <mpi.h>
#endif

//Links omp
//#define USE_OpenMP
#ifdef USE_OpenMP
	#include <omp.h>
#endif



//#define COMPILE_WITH_CUDA
//#define COMPILE_WITH_HIP
//#define USE_GPU_HIP  // <**** temporary

//#define UseCUDA
//#define UseHIP


//Links HIP
#ifdef COMPILE_WITH_HIP
    #include "hip/hip_runtime.h"
    #include "hip/hip_runtime_api.h"
#endif

//Links CUDA
//#include "cuda_runtime.h"
//#include "cuda.h"





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


#define _param(...) _parameters=SUBTASK::parameters(__VA_ARGS__)



//Links internal
#include "hardware.hpp"
#include "taskgpu.hpp"
#include "taskcpu.hpp"

//================================================================================================================================
// THE END.
//================================================================================================================================

