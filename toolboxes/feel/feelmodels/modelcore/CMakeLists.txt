#############################################################################

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/feelmodelscoreconstconfig.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/feelmodelscoreconstconfig.hpp )
INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/feelmodelscoreconstconfig.hpp DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelcore COMPONENT Devel)

#############################################################################
#############################################################################

add_library(feelpp_modelcore
  SHARED
  options.cpp log.cpp timertool.cpp
  modelbase.cpp modelalgebraic.cpp modelnumerical.cpp
  markermanagement.cpp modelmeasures.cpp
  )
target_link_libraries(feelpp_modelcore Feelpp::feelpp)
target_include_directories(feelpp_modelcore PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
set_property(TARGET feelpp_modelcore PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
  add_precompiled_header( feelpp_modelcore )
endif()

install(TARGETS feelpp_modelcore EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
  )
FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelcore COMPONENT Devel)
