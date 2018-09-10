#######################################################
# mesh lib
#######################################################

set(FEELPP_MODELS_MESH_SRC
  createmesh.cpp
  fsimesh.cpp
)
add_library(
  feelpp_modelmesh
  SHARED
  ${FEELPP_MODELS_MESH_SRC}
)
target_link_libraries(feelpp_modelmesh PUBLIC feelpp_modelcore )
target_include_directories(feelpp_modelmesh PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
set_property(TARGET feelpp_modelmesh PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
    add_precompiled_header( feelpp_modelmesh )
endif()
install(TARGETS feelpp_modelmesh EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
  )
#######################################################
# meshale lib
#######################################################
if( FEELPP_MODELS_HAS_MESHALE )
  # harmonic extension src
  if( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
    set(HARMONICEXTENSION_SRC
      harmonicextension.cpp
      )
  endif()
  # winslow src
  if( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
    set(WINSLOW_SRC
      winslow_inst2d.cpp
      winslow_inst3d.cpp
      )
  endif()

  # lib src
  set(FEELPP_MODELS_MESHALE_SRC
    ${HARMONICEXTENSION_SRC}
    ${WINSLOW_SRC}
    ale2d.cpp 
    ale3d.cpp
    )
  add_library( feelpp_modelmeshale
    SHARED
    ${FEELPP_MODELS_MESHALE_SRC}
    )
  set( FEELPP_MODELS_MESHALE_LIB_DEPENDS feelpp_modelalg feelpp_modelcore )
  target_link_libraries( feelpp_modelmeshale
    ${FEELPP_MODELS_MESHALE_LIB_DEPENDS}
    )
  set_property(TARGET feelpp_modelmeshale PROPERTY MACOSX_RPATH ON)
  INSTALL(TARGETS feelpp_modelmeshale DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)
  
endif( FEELPP_MODELS_HAS_MESHALE )

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelmesh COMPONENT Devel)
