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
target_link_libraries(feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES})
set_property(TARGET feelpp_modelcore PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
    add_precompiled_header( feelpp_modelcore )
endif()
INSTALL(TARGETS feelpp_modelcore DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)
FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelcore COMPONENT Devel)
