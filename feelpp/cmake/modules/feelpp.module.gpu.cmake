# Detect GPU support
option(FEELPP_ENABLE_HIP "Enable HIP support" OFF)
option(FEELPP_ENABLE_ROCM "Enable ROCm support" OFF)
option(FEELPP_ENABLE_CUDA "Enable CUDA support" OFF)





if(FEELPP_ENABLE_ROCM)
    message(STATUS "[feelpp] ROCm support enabled")
    find_package(HIP REQUIRED)
    if ( HIP_FOUND )
        message(STATUS "HIP found")
        enable_language(HIP)
        add_definitions(-DUSE_HIP)
        set(CMAKE_HIP_STANDARD ${CPPSTD})
        set(CMAKE_HIP_STANDARD_REQUIRED ON)
    endif()
    
    find_package(rocblas REQUIRED)
    find_package(rocthrust REQUIRED)
    # Add HIP-specific settings
    set(FEELPP_ENABLE_GPU "rocm" PARENT_SCOPE FORCE)
endif()

if(FEELPP_ENABLE_CUDA)
    find_package(CUDA REQUIRED)
    add_definitions(-DUSE_CUDA)
    # Add CUDA-specific settings
    set(FEELPP_ENABLE_GPU "cuda" PARENT_SCOPE FORCE)
endif()
set(FEELPP_ENABLE_GPU "cpu" PARENT_SCOPE)
