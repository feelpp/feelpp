include(FetchContent)
set(FETCHCONTENT_SOURCE_DIR_HANA ${CMAKE_CURRENT_SOURCE_DIR}/hana)
FetchContent_Declare(hana GIT_REPOSITORY https://github.com/feelpp/hana.git GIT_TAG 47b0d1f088d9d3a23dfafbc9cf228519f65c67fd) # the HEAD of feelpp branch
FetchContent_MakeAvailable(hana)

target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_BOOST_HANA )
#target_link_libraries(feelpp_contrib INTERFACE feelpp_hana )
target_include_directories(feelpp_contrib BEFORE INTERFACE
  $<BUILD_INTERFACE:${hana_SOURCE_DIR}/include>
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
