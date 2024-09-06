###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 30 Aug 2018
#
#  Copyright (C) 2018 Feel++ Consortium
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
project(Feelpp)





SET( CMAKE_MODULE_PATH  ${FEELPP_HOME_DIR}/feelpp/cmake/modules ${FEELPP_HOME_DIR}/feelpp/cmake/machines )

include(feelpp.dirs)

# Use of standard headers for GiNaC compiler (for supermuc this option should be OFF)
option(USE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER "Enable use of standard include headers in ginac-excompiler" ON)

OPTION(FEELPP_ENABLE_HOMEBREW "Enable Feel++/Homebrew support" OFF)
if (APPLE)
  set(FEELPP_ENABLE_HOMEBREW ON)
endif()

# possibly set the Feel++ configuration
include(feelpp.machines.config)

# RelWithDebInfo  is the default value for CMAKE_BUILD_TYPE
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

if(APPLE)
  set(FEELPP_OS "MacOSX")
elseif(CYGWIN)
  set(FEELPP_OS "Windows")
else(APPLE)
  set(FEELPP_OS "${CMAKE_SYSTEM_NAME}")
endif(APPLE)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-1024 -Wno-inline" CACHE STRING "Default flags" )
  SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fPIC " CACHE STRING "Debug flags" FORCE)
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O1" CACHE STRING "Release With Debug information flags" FORCE)
  SET(CMAKE_CXX_FLAGS_ASAN "-g3 -O0  -fsanitize=address,leak -fno-omit-frame-pointer -fno-optimize-sibling-calls " CACHE STRING "Asan flags" FORCE)
  SET(CMAKE_CXX_FLAGS_DEBUGFULL "-g -O0 -fno-omit-frame-pointer -fno-optimize-sibling-calls " CACHE STRING "Debug Full flags" FORCE)
  SET(CMAKE_CXX_FLAGS_COVERAGE "-g -O0 --coverage " CACHE STRING "Coverage reporting for gcc/clang" FORCE)
endif()


#  set(CTEST_BUILD_FLAGS -j${N2})
#  set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N2})

# set test timeout to 300s
set(CTEST_TIMEOUT "20")

#else()
#  OPTION(FEELPP_ENABLE_MPI_MODE "Instantiation mode" OFF )
#endif()


if (TRUE)
    MESSAGE("{COMPILE_WITH_HIP STEP 0}")
    if(NOT DEFINED HIP_PATH)
        if(NOT DEFINED ENV{HIP_PATH})
            set(HIP_PATH "/opt/rocm-5.7.2/hip" CACHE PATH "Path to which HIP has been installed")
        else()
            set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
        endif()
    endif()

    if(NOT DEFINED HCC_PATH)
      if(DEFINED ENV{HCC_PATH})
          set(HCC_PATH $ENV{HCC_PATH} CACHE PATH "Path to which HCC has been installed")
      else()
          set(HCC_PATH "${ROCM_PATH}/hcc" CACHE PATH "Path to which HCC has been installed")
      endif()
      set(HCC_HOME "${HCC_PATH}")
    endif()

    if(NOT DEFINED HIP_CLANG_PATH)
        if(NOT DEFINED ENV{HIP_CLANG_PATH})
            set(HIP_CLANG_PATH "${ROCM_PATH}/llvm/bin" CACHE PATH "Path to which HIP compatible clang binaries have been installed")
        else()
            set(HIP_CLANG_PATH $ENV{HIP_CLANG_PATH} CACHE PATH "Path to which HIP compatible clang binaries have been installed")
        endif()
    endif()


    set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})

    list(APPEND CMAKE_PREFIX_PATH
    "${HIP_PATH}/lib/cmake"
    "${HIP_PATH}/../lib/cmake" # hopefully catches all extra HIP dependencies
    )


    find_package(BLAS)
    find_package(LAPACK)
    find_package(hip)
    find_package(HIP QUIET)

    if(HIP_FOUND)
        message(STATUS "Found HIP: " ${HIP_VERSION})
    else()
        message(FATAL_ERROR "Could not find HIP. Ensure that HIP is either installed in /opt/rocm/hip or the variable HIP_PATH is set to point to the right location.")
    endif()


    message(STATUS "HIP Path: " ${HIP_PATH})

   

    if($ENV{VERBOSE})
            MESSAGE(STATUS "GEMM -- BLAS_LIBRARIES   : ${BLAS_LIBRARIES}")
            MESSAGE(STATUS "GEMM -- BLAS_LINKER_FLAGS   : ${BLAS_LINKER_FLAGS}")
            MESSAGE(STATUS "GEMM -- LAPACK_LIBRARIES   : ${LAPACK_LIBRARIES}")
            MESSAGE(STATUS "GEMM -- LAPACK_LINKER_FLAGS   : ${LAPACK_LINKER_FLAGS}")
    endif()

    MESSAGE("")
    MESSAGE("{BEGIN::COMPILE_WITH_HIP}")

    include_directories(${HIP_PATH}/include)
    include_directories(/opt/rocm-5.7.2/hip/include)
    include_directories(/opt/rocm-5.7.2/hipblas/include)
    include_directories(/opt/rocm-5.7.2/hipblas/lib)

    list(APPEND CMAKE_PREFIX_PATH /opt/rocm-5.7.2 /opt/rocm-5.7.2)

    # Set the GPU to compile for
    #set(GPU_TARGETS "gfx900:gfx906:gfx1032;gfx1035" CACHE STRING "GPU targets to compile for")

    set(GPU_RUNTIME "HIP" CACHE STRING "Switches between HIP and CUDA")
    set(GPU_RUNTIMES "HIP" "CUDA")
    set_property(CACHE GPU_RUNTIME PROPERTY STRINGS ${GPU_RUNTIMES})

    if(NOT "${GPU_RUNTIME}" IN_LIST GPU_RUNTIMES)
        set(ERROR_MESSAGE "GPU_RUNTIME is set to \"${GPU_RUNTIME}\".\nGPU_RUNTIME must be either HIP or CUDA.")
        message(FATAL_ERROR ${ERROR_MESSAGE})
    endif()

    #enable_language(${GPU_RUNTIME})

    set(CMAKE_${GPU_RUNTIME}_STANDARD 17)
    set(CMAKE_${GPU_RUNTIME}_EXTENSIONS OFF)
    set(CMAKE_${GPU_RUNTIME}_STANDARD_REQUIRED ON)

    set(ROCM_ROOT "/opt/rocm-5.7.2" CACHE PATH "Root directory of the ROCm installation")
    list(APPEND CMAKE_PREFIX_PATH "${ROCM_ROOT}")



    if(FALSE)
      find_package(rocprim REQUIRED)
      if (rocprim_FOUND)
        message("{*************************** rocprim_LIBRARIES *************************************}") 
      else()
        message(FATAL_ERROR "xxxxxxxxxxxxxxxx Could not find rocprim xxxxxxxxxxxxxxxx")     
      endif()
    endif()

    if(FALSE)
      find_package(rocthrust REQUIRED)
      if (rocthrust_FOUND)
        message("{*************************** rocthrust_LIBRARIES *************************************}") 
      else()
        message(FATAL_ERROR "xxxxxxxxxxxxxxxx Could not find rocThrust xxxxxxxxxxxxxxxx")     
      endif()
    endif()


    if(TRUE)
      find_package(rocprim CONFIG PATHS "/nvme0/lemoinep/soft/rocPRIM")
      if (rocprim_FOUND)
        message("{*************************** rocprim_LIBRARIES *************************************}") 
        set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/nvme0/lemoinep/soft/rocPRIM/include")
      else()
        message(FATAL_ERROR "xxxxxxxxxxxxxxxx Could not find rocprim xxxxxxxxxxxxxxxx")     
      endif()
    endif()

    if(TRUE)
      #find_package(rocthrust REQUIRED)
      #find_package(rocthrust REQUIRED CONFIG PATHS "/opt/rocm-5.7.2/rocthrust") in future
      find_package(rocthrust REQUIRED CONFIG PATHS "/nvme0/lemoinep/soft/rocThrust")
      if (rocthrust_FOUND)
        message("{*************************** rocthrust_LIBRARIES *************************************}") 
        #list(APPEND include_dirs "/nvme0/lemoinep/soft/rocThrust/include")
        #list(APPEND CMAKE_PREFIX_PATH /nvme0/lemoinep/soft/rocThrust/include)
        #include_directories(/nvme0/lemoinep/soft/rocThrust/include)
        #include_directories(/opt/rocm-5.7.2/rocthrust/include)
        set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/nvme0/lemoinep/soft/rocThrust/include")
      else()
        message(FATAL_ERROR "xxxxxxxxxxxxxxxx Could not find rocThrust xxxxxxxxxxxxxxxx")     
      endif()
    endif()

    add_compile_definitions( __HIP_PLATFORM_AMD__)
    #list(APPEND GPU_HIP_LIBRARIES hip::device hip::host hipblas)
    list(APPEND GPU_HIP_LIBRARIES hip::device hip::host hipblas roc::rocthrust)
    list(APPEND GPU_ROCTHRUST_LIBRARIES roc::rocthrust)

    if(GPU_RUNTIME STREQUAL "HIP")
      list(APPEND include_dirs "${ROCM_ROOT}/include")
      list(APPEND include_dirs "/nvme0/lemoinep/soft/rocPRIM/include")
      list(APPEND include_dirs "/nvme0/lemoinep/soft/rocThrust/include")
      MESSAGE(${ROCM_ROOT}/include)
    endif()

    MESSAGE("{END::COMPILE_WITH_HIP}")
    MESSAGE("")
endif()






if ( FEELPP_BENCHMARK_FLAGS )
  set(CMAKE_BUILD_TYPE Release )
  set(GCC_PARAM_INLINE_UNIT_GROWTH 150)
  set(GCC_PARAM_MAX_INLINE_INSNS_SINGLE 500)
  set(GCC_PARAM_LARGE_FUNCTION_GROWTH 600)
  #  SET(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wshadow -Woverloaded-virtual -std=c++0x -O3 -DNDEBUG --param max-inline-recursive-depth=256 --param max-gcse-memory=8000 --param max-inline-insns-single=${GCC_PARAM} --param inline-unit-growth=${GCC_PARAM} --param large-unit-insns=${GCC_PARAM} --param large-function-growth=${GCC_PARAM} --param large-function-insns=${GCC_PARAM} " CACHE STRING "Benchmarks Release flags" FORCE)
  #  SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -O3 -DNDEBUG --param max-inline-recursive-depth=256 --param max-gcse-memory=8000 --param max-inline-insns-single=${GCC_PARAM} --param inline-unit-growth=${GCC_PARAM} --param large-unit-insns=${GCC_PARAM} --param large-function-growth=${GCC_PARAM} --param large-function-insns=${GCC_PARAM} " CACHE STRING "Benchmarks Release flags" FORCE)
  SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -O3 -DNDEBUG --param max-inline-insns-single=${GCC_PARAM_MAX_INLINE_INSNS_SINGLE} --param inline-unit-growth=${GCC_PARAM_INLINE_UNIT_GROWTH} --param large-function-growth=${GCC_PARAM_LARGE_FUNCTION_GROWTH} " CACHE STRING "Benchmarks Release flags" FORCE)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  option(FEELPP_ENABLE_SSE2 "Enable/Disable SSE2 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    message(STATUS "[feelpp] Enabling SSE2 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE3 "Enable/Disable SSE3 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
    message(STATUS "[feelpp] Enabling SSE3 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSSE3 "Enable/Disable SSSE3 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSSE3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mssse3")
    message(STATUS "[feelpp] Enabling SSSE3 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE4_1 "Enable/Disable SSE4.1 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE4_1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
    message(STATUS "[feelpp] Enabling SSE4.1 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE4_2 "Enable/Disable SSE4.2 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE4_2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
    message(STATUS "[feelpp] Enabling SSE4.2 in tests/examples")
  endif()

  option(FEELPP_ENABLE_ALTIVEC "Enable/Disable AltiVec in tests/examples" OFF)
  if(FEELPP_ENABLE_ALTIVEC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -maltivec -mabi=altivec")
    message(STATUS "[feelpp] Enabling AltiVec in tests/examples")
  endif()

  option(FEELPP_ENABLE_NEON "Enable/Disable Neon in tests/examples" OFF)
  if(FEELPP_ENABLE_NEON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfloat-abi=softfp -mfpu=neon -mcpu=cortex-a8")
    message(STATUS "[feelpp] Enabling NEON in tests/examples")
  endif()
endif(CMAKE_COMPILER_IS_GNUCXX)
#INCLUDE(PackageArchGlobalMacros)
# INCLUDE(FeelGlobalMacros)
#INCLUDE(AdvancedSet)
#INCLUDE(AdvancedOption)

include(feelpp.extra.warnings)
include(feelpp.extra.astyle)
add_definitions(${FEELPP_FLAGS})

if("${CMAKE_CXX_FLAGS} ${FEELPP_FLAGS}" MATCHES "[^ ]")
  message(STATUS "[feelpp] [feel++] Global flags: ${CMAKE_CXX_FLAGS} ${FEELPP_FLAGS}")
endif()

message(STATUS "[feelpp] [feel++] Debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "[feelpp] [feel++] Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")

INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckLibraryExists)
INCLUDE(ParseArguments)

# MACRO(CAR var)
#   SET(${var} ${ARGV1})
# ENDMACRO(CAR)

# MACRO(CDR var junk)
#   SET(${var} ${ARGN})
# ENDMACRO(CDR)
INCLUDE(CheckTypeSize)
CHECK_TYPE_SIZE(int SIZE_INT )
CHECK_TYPE_SIZE(uint SIZE_UINT )
CHECK_TYPE_SIZE(size_t SIZE_SIZE_T )
CHECK_TYPE_SIZE(long SIZE_LONG )
CHECK_TYPE_SIZE(float SIZE_FLOAT )
CHECK_TYPE_SIZE(double SIZE_DOUBLE )
CHECK_TYPE_SIZE("long double" SIZE_LONG_DOUBLE)
MESSAGE(STATUS "[feelpp] SIZE_INT=${SIZE_INT}")
MESSAGE(STATUS "[feelpp] SIZE_UINT=${SIZE_UINT}")
MESSAGE(STATUS "[feelpp] SIZE_SIZE_T=${SIZE_SIZE_T}")
MESSAGE(STATUS "[feelpp] SIZE_LONG=${SIZE_LONG}")
MESSAGE(STATUS "[feelpp] SIZE_FLOAT=${SIZE_FLOAT}")
MESSAGE(STATUS "[feelpp] SIZE_DOUBLE=${SIZE_DOUBLE}")
MESSAGE(STATUS "[feelpp] SIZE_LONG_DOUBLE=${SIZE_LONG_DOUBLE}")

include( feelpp.macros )

include( feelpp.dependencies )

# Find feel++ dependencies 
add_custom_target(contrib)

add_subdirectory(contrib)
set(FEELPP_ENABLED_OPTIONS ${FEELPP_ENABLED_OPTIONS} ${FEELPP_ENABLED_OPTIONS_CONTRIB} PARENT_SCOPE )
feelpp_set_options( Feelpp::feelpp_contrib "Feel++" )

if(NOT FEELPP_ENABLE_BUILD_STATIC)
	set(CMAKE_C++_CREATE_SHARED_LIBRARY 1)
endif()

option(FEELPP_ENABLE_JSON_DIAGNOSTICS "enable extended json diagnostics" ON)

# feel++ lib
ADD_SUBDIRECTORY ( feel )

# get feel++ environment of cmake variable (FEELPP_HAS_...)
feelpp_get_environment()

#install cmake files
add_subdirectory(cmake)

OPTION(FEELPP_ENABLE_QUICKSTART "enable Feel++ quickstart" ON)
OPTION(FEELPP_ENABLE_PYFEELPP_LIBFEELPP "enable PyFeel++" ON)

add_subdirectory(tools)
include(feelpp.module.quickstart)
if ( FEELPP_HAS_QUICKSTART )
  add_custom_target(quickstart)
  add_subdirectory(quickstart)
endif()



# The configuration File needs to be done at the very end otherwise we won't have the correct values for headers and libraries
# The add_subdirectory is way too soon to get those values
# The install command of this file is located in cmake/module/CMakeLists (because of feelpp-install-lib target)
set(FEELPP_CONFIG_LIB_FILE ${FEELPP_BINARY_DIR}/feelpp/cmake/modules/Feel++Config.cmake )
configure_file(${FEELPP_SOURCE_DIR}/feelpp/cmake/modules/Feel++Config.cmake.in ${FEELPP_CONFIG_LIB_FILE} @ONLY)
configure_file(${FEELPP_SOURCE_DIR}/feelpp/cmake/modules/feelpp.install.config.cmake.in ${FEELPP_BINARY_DIR}/feelpp/cmake/modules/feelpp.install.config.cmake @ONLY)


