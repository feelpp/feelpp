set( FEELPP_TOOLBOXES_MATERIALS_SRC materialsproperties.cpp )

add_library(feelpp_modelmaterials
  SHARED
  ${FEELPP_TOOLBOXES_MATERIALS_SRC}
  )
target_link_libraries(feelpp_modelmaterials Feelpp::feelpp)
add_library(Feelpp::feelpp_modelmaterials ALIAS feelpp_modelmaterials )  # to match exported target
set_target_properties(feelpp_modelmaterials PROPERTIES OUTPUT_NAME "feelpp_modelmaterials")
target_include_directories(feelpp_modelmaterials PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
set_property(TARGET feelpp_modelmaterials PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
  add_precompiled_header( feelpp_modelmaterials )
endif()
include(GNUInstallDirs)
install(TARGETS feelpp_modelmaterials EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
set_target_properties(feelpp_modelmaterials PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelmaterials COMPONENT Devel)
