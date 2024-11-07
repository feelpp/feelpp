if(FEELPP_ENABLE_ROCM)
    message(STATUS "[feelpp] ROCm support enabled")
    find_package(HIP REQUIRED)
    if ( HIP_FOUND )
        message(STATUS "HIP found")
        enable_language(HIP)
        set(CMAKE_HIP_STANDARD ${CPPSTD})
        message(STATUS "HIP standard: ${CMAKE_HIP_STANDARD}")
        set(CMAKE_HIP_STANDARD_REQUIRED ON)
    endif()
    
    find_package(rocblas REQUIRED)
    find_package(rocthrust REQUIRED)
    # Add HIP-specific settings
    set(FEELPP_ENABLE_GPU "rocm"   PARENT_SCOPE)
endif()

if(FEELPP_ENABLE_CUDA)
    find_package(CUDA REQUIRED)
    # Add CUDA-specific settings
    set(FEELPP_ENABLE_GPU "cuda"   PARENT_SCOPE )
endif()

add_library(feelpp_gpu INTERFACE)
add_library(Feelpp::feelpp_gpu ALIAS feelpp_gpu)
if (hip_FOUND)
  target_link_libraries(feelpp_gpu INTERFACE hip::device hip::host hipblas roc::rocthrust)
  target_compile_definitions(feelpp_gpu INTERFACE FEELPP_HAS_HIP FEELPP_HAS_ROCM)
endif()
if(CUDA_FOUND)
    target_link_libraries(Feelpp::feelpp_gpu INTERFACE CUDA::CUDA)
    target_compile_definitions(Feelpp::feelpp_gpu INTERFACE FEELPP_HAS_CUDA)
else()
    message(STATUS "CUDA not found.")
endif()


