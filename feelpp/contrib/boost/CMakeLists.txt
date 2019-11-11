add_subdirectory( hana )
target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_BOOST_HANA )
#target_link_libraries(feelpp_contrib INTERFACE feelpp_hana )
target_include_directories(feelpp_contrib INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/boost/hana/include>
  $<INSTALL_INTERFACE:include/feelpp/boost/hana> )
set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Boost/Hana")

install(DIRECTORY fusion DESTINATION include/feelpp/contrib/boost )
