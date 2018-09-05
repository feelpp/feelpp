
add_library( feelpp_modelalg SHARED
  modelalgebraicfactory.cpp )
target_link_libraries (feelpp_modelalg PUBLIC feelpp_modelcore )
target_include_directories(feelpp_modelalg PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
set_property(TARGET feelpp_modelalg PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
  add_precompiled_header( feelpp_modelalg )
endif()
install(TARGETS feelpp_modelalg EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
  )
FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelalg COMPONENT Devel)
