add_subdirectory( hana )
target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_BOOST_HANA )
#target_link_libraries(feelpp_contrib INTERFACE feelpp_hana )
target_include_directories(feelpp_contrib INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/boost/hana/include>
  $<INSTALL_INTERFACE:include/feelpp/boost/hana> )
set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Boost/Hana")

install(DIRECTORY fusion DESTINATION include/feelpp/contrib/boost )

foreach( _patch_boostversion "1.69.0" "1.70.0")
  if ( Boost_VERSION_STRING VERSION_EQUAL _patch_boostversion )
    target_include_directories(feelpp_contrib INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/boost/v${_patch_boostversion}>
      $<INSTALL_INTERFACE:include/feelpp/boost/v${_patch_boostversion}> )
    install(DIRECTORY v${_patch_boostversion} DESTINATION include/feelpp/boost COMPONENT Devel)
  endif()
endforeach()
