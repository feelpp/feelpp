
add_library(
  feelpp_toolbox_heatbase
  SHARED
  heatboundaryconditions.cpp
)
target_link_libraries(feelpp_toolbox_heatbase PUBLIC feelpp_modelcore )
add_library(Feelpp::feelpp_toolbox_heatbase ALIAS feelpp_toolbox_heatbase )  # to match exported target
set_target_properties(feelpp_toolbox_heatbase PROPERTIES OUTPUT_NAME "feelpp_toolbox_heatbase")
target_include_directories(feelpp_toolbox_heatbase PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
set_property(TARGET feelpp_toolbox_heatbase PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
  add_precompiled_header( feelpp_toolbox_heatbase )
endif()
#include(GNUInstallDirs)
set_target_properties(feelpp_toolbox_heatbase PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_heatbase EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )








# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  message(STATUS "[feelpp] HEAT TRANSFER ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    genLibHeat(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_HEAT  ${FEELPP_MODELS_LIBRARIES_HEAT} ${HEAT_LIB_NAME}  )
  endforeach()
endforeach()

# isoparametric
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_T_ORDER})
  message(STATUS "[feelpp] HEAT TRANSFER ${_DIM}D PNGN N=2..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
    genLibHeat(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      GEO_ORDER ${_T_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_HEAT  ${FEELPP_MODELS_LIBRARIES_HEAT} ${HEAT_LIB_NAME} )
  endforeach()
endforeach()

add_library(feelpp_toolbox_heat_lib SHARED libheat.cpp )
add_library(Feelpp::feelpp_toolbox_heat_lib ALIAS feelpp_toolbox_heat_lib )  # to match exported target
set_target_properties(feelpp_toolbox_heat_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_heat_lib")
target_link_libraries(feelpp_toolbox_heat_lib ${FEELPP_MODELS_LIBRARIES_HEAT})
set_property(TARGET feelpp_toolbox_heat_lib PROPERTY MACOSX_RPATH ON)
set_target_properties(feelpp_toolbox_heat_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_heat_lib EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/heat COMPONENT Devel)
