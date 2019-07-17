
add_library( feelpp_modelalg SHARED
  modelalgebraicfactory.cpp )
add_library(Feelpp::feelpp_modelalg ALIAS feelpp_modelalg )  # to match exported target
set_target_properties(feelpp_modelalg PROPERTIES OUTPUT_NAME "feelpp_modelalg")
target_link_libraries (feelpp_modelalg PUBLIC feelpp_modelcore )
target_include_directories(feelpp_modelalg PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
set_property(TARGET feelpp_modelalg PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
  add_precompiled_header( feelpp_modelalg )
endif()
include(GNUInstallDirs)
install(TARGETS feelpp_modelalg EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelalg COMPONENT Devel)
