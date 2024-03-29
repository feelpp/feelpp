#############################################################################

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/feelmodelscoreconstconfig.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/feelmodelscoreconstconfig.hpp )
INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/feelmodelscoreconstconfig.hpp DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelcore COMPONENT Devel)

#############################################################################
#############################################################################


set( FEELPP_TOOLBOXES_CORE_SRC  options.cpp log.cpp timertool.cpp
  modelbase.cpp modelalgebraic.cpp modelnumerical.cpp modelphysics.cpp modelgenericpde.cpp modelmeshes.cpp modelmeshadaptation.cpp
  modelalgebraicfactory.cpp markermanagement.cpp modelmeasures.cpp genericboundaryconditions.cpp )

if( FEELPP_MODELS_HAS_MESHALE )
  # harmonic extension src
  if( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
    list(APPEND FEELPP_TOOLBOXES_CORE_SRC ../modelmesh/harmonicextension.cpp )
  endif()
  # winslow src
  if( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
    list(APPEND FEELPP_TOOLBOXES_CORE_SRC  ../modelmesh/winslow_inst2d.cpp ../modelmesh/winslow_inst3d.cpp )
  endif()

  # lib src
  list(APPEND FEELPP_TOOLBOXES_CORE_SRC  ../modelmesh/ale2d.cpp ../modelmesh/ale3d.cpp )
endif()

add_library(feelpp_modelcore
  SHARED
  ${FEELPP_TOOLBOXES_CORE_SRC}
  )
target_link_libraries(feelpp_modelcore Feelpp::feelpp)
add_library(Feelpp::feelpp_modelcore ALIAS feelpp_modelcore )  # to match exported target
set_target_properties(feelpp_modelcore PROPERTIES OUTPUT_NAME "feelpp_modelcore")
target_include_directories(feelpp_modelcore PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
set_property(TARGET feelpp_modelcore PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
  add_precompiled_header( feelpp_modelcore )
endif()
include(GNUInstallDirs)
install(TARGETS feelpp_modelcore EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
set_target_properties(feelpp_modelcore PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelcore COMPONENT Devel)
