add_library( feelpp_modelvf INTERFACE )
add_library( Feelpp::feelpp_modelvf ALIAS feelpp_modelvf )  # to match exported target

target_include_directories(feelpp_modelvf INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
install(TARGETS feelpp_modelvf EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
  )

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelvf COMPONENT Devel)
