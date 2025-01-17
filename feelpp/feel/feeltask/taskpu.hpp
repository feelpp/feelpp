/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

  @file
  @author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  @Date: 2024-06-04
  @copyright 2019 Feel++ Consortium
*/

// Enable compiler warnings for unused results and unknown attributes
#pragma GCC diagnostic warning "-Wunused-result"
#pragma clang diagnostic ignored "-Wunused-result"

#pragma GCC diagnostic warning "-Wunknown-attributes"
#pragma clang diagnostic ignored "-Wunknown-attributes"

// Uncomment to enable specific compilation options
//#define COMPILE_WITH_CUDA       // Enable CUDA support
//#define COMPILE_WITH_HIP        // Enable HIP support (AMD)
//#define COMPILE_WITH_CXX_20     // Enable C++20 features

//#define USE_MPI                // Enable MPI support for parallel computing
//#define USE_OpenMP             // Enable OpenMP support for multi-threading
//#define USE_GPU_HIP            // Temporary definition for HIP GPU support
//#define UseCUDA                // Define to use CUDA APIs
//#define UseHIP                 // Define to use HIP APIs
//#define USE_jthread            // Define to use jthread (C++20)

// Include necessary headers for various functionalities
#include <specx/Data/SpDataAccessMode.hpp>
#include <specx/Legacy/SpRuntime.hpp>
#include <specx/Task/SpPriority.hpp>
#include <specx/Task/SpProbability.hpp>
#include <specx/Utils/SpArrayView.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpConsumerThread.hpp>
#include <specx/Utils/SpHeapBuffer.hpp>
#include <specx/Utils/SpTimer.hpp>
#include <specx/Utils/SpUtils.hpp>
#include <specx/Utils/small_vector.hpp>

#include <napp/na.hpp> // Include NAPP library for numerical applications

#ifdef USE_MPI
#include <mpi.h>      // Include MPI header if MPI support is enabled
#endif

#ifdef USE_OpenMP
#include <omp.h>      // Include OpenMP header if OpenMP support is enabled
#endif

#ifdef COMPILE_WITH_HIP
#include "hip/hip_runtime.h"          // Include HIP runtime if HIP support is enabled
#include "hip/hip_runtime_api.h"      // Include HIP runtime API definitions
#endif

#ifdef COMPILE_WITH_CUDA
#include "cuda.h"                      // Include CUDA header if CUDA support is enabled
#include "cuda_runtime.h"              // Include CUDA runtime definitions
#endif

// Macro functions for AMD HIP error handling and kernel naming 
#define HIP_KERNEL_NAME( ... ) __VA_ARGS__

#define HIP_CHECK( command )                                                                \
    {                                                                                       \
        hipError_t status = command;                                                        \
        if ( status != hipSuccess )                                                         \
        {                                                                                   \
            std::cerr << "Error: HIP reports " << hipGetErrorString( status ) << std::endl; \
            std::abort();                                                                   \
        }                                                                                   \
    }

// Conditional assertion based on debug mode for HIP operations 
#ifdef NDEBUG
#define HIP_ASSERT( x ) x                  // No assertion in release mode
#else
#define HIP_ASSERT( x ) ( assert( ( x ) == hipSuccess ) ) // Assert in debug mode
#endif

// Macro functions for NVIDIA CUDA error handling 
#define CUDA_CHECK( command )                                                                 \
    {                                                                                         \
        cudaError_t status = command;                                                         \
        if ( status != cudaSuccess )                                                          \
        {                                                                                     \
            std::cerr << "Error: CUDA reports " << cudaGetErrorString( status ) << std::endl; \
            std::abort();                                                                     \
        }                                                                                     \
    }

// Conditional assertion based on debug mode for CUDA operations 
#ifdef NDEBUG
#define CUDA_ASSERT( x ) x                  // No assertion in release mode
#else
#define CUDA_ASSERT( x ) ( assert( ( x ) == cudaSuccess ) ) // Assert in debug mode
#endif

// Macro to set parameters for tasks 
#define _param( ... ) _parameters = Sbtask::parameters( __VA_ARGS__ )

// Internal links for hardware and task management 
#include "hardware.hpp"   // Hardware-specific configurations and utilities 
#include "taskcpu.hpp"    // CPU task management functionalities 
#include "taskgpu.hpp"    // GPU task management functionalities 

// links internal scheduling management
//...


