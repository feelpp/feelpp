
project(shading_masks_gpu)

include(CheckLanguage)
check_language(HIP)
if(CMAKE_HIP_COMPILER)
   enable_language(HIP)
endif()

find_package(hip)


# A VOIR SI ...
#set(SOLARSHADING_TARGET_NAME feelpp_solarshading)
#add_library(${SOLARSHADING_TARGET_NAME} INTERFACE)
#target_include_directories(${SOLARSHADING_TARGET_NAME} INTERFACE
#  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#  $<INSTALL_INTERFACE:include>
#)
#target_link_libraries(${SOLARSHADING_TARGET_NAME} INTERFACE Feelpp::feelpp)
#feelpp_add_application(shading_masks_gpu SRCS shading_masks_gpu.cpp
#  DEFS FEELPP_TOP_DIM=3
#  LINK_LIBRARIES ${SOLARSHADING_TARGET_NAME} ${GPU_HIP_LIBRARIES} Feelpp::feelpp  hip::device hip::host hipblas)
#  set_source_files_properties(shading_masks_gpu.cpp 
#  PROPERTIES LANGUAGE ${GPU_RUNTIME})


#3D OK
if(TRUE)
feelpp_add_application(shading_masks_gpu_3D SRCS shading_masks_gpu.cpp
  DEFS FEELPP_TOP_DIM=3 
  LINK_LIBRARIES ${GPU_HIP_LIBRARIES} Feelpp::feelpp  hip::device hip::host hipblas)
  set_source_files_properties(shading_masks_gpu.cpp 
  PROPERTIES LANGUAGE ${GPU_RUNTIME})
endif()


if(TRUE)
feelpp_add_application(shading_masks_gpu_2D SRCS shading_masks_gpu.cpp
  DEFS FEELPP_TOP_DIM=2 
  LINK_LIBRARIES ${GPU_HIP_LIBRARIES} Feelpp::feelpp  hip::device hip::host hipblas)
  set_source_files_properties(shading_masks_gpu.cpp 
  PROPERTIES LANGUAGE ${GPU_RUNTIME})
endif()




if(TRUE)
    feelpp_add_testcase(cases)
endif()


