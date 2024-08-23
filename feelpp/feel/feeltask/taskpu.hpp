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

#pragma GCC diagnostic warning "-Wunused-result"
#pragma clang diagnostic ignored "-Wunused-result"

#pragma GCC diagnostic warning "-Wunknown-attributes"
#pragma clang diagnostic ignored "-Wunknown-attributes"

//#define COMPILE_WITH_CUDA
//#define COMPILE_WITH_HIP
//#define COMPILE_WITH_CXX_20

//#define USE_MPI
//#define USE_OpenMP
//#define USE_GPU_HIP  // <**** temporary
//#define UseCUDA
//#define UseHIP
//#define USE_jthread

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
#include "cuda.h"
#include "cuda_runtime.h"
#endif

// Some macro functions for the AMD HIP
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

#ifdef NDEBUG
#define HIP_ASSERT( x ) x
#else
#define HIP_ASSERT( x ) ( assert( ( x ) == hipSuccess ) )
#endif

// Some macro functions for the NVida CUDA
#define CUDA_CHECK( command )                                                                 \
    {                                                                                         \
        cudaError_t status = command;                                                         \
        if ( status != hipSuccess )                                                           \
        {                                                                                     \
            std::cerr << "Error: CUDA reports " << cudaGetErrorString( status ) << std::endl; \
            std::abort();                                                                     \
        }                                                                                     \
    }

#ifdef NDEBUG
#define CUDA_ASSERT( x ) x
#else
#define CUDA_ASSERT( x ) ( assert( ( x ) == cudaSuccess ) )
#endif

#define _param( ... ) _parameters = Sbtask::parameters( __VA_ARGS__ )

// Links internal
#include "hardware.hpp"
#include "taskcpu.hpp"
#include "taskgpu.hpp"

// links internal scheduling management
//...
