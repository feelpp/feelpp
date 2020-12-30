#add_subdirectory( hana )
include(FetchContent)
FetchContent_Declare(hana GIT_REPOSITORY https://github.com/boostorg/hana.git GIT_TAG 77b0f5ad6dad24717599e8323705149fc0f01c24) # the commit hash for develop branch Dec 28 2020
# hana requires workaround for cmake
FetchContent_GetProperties(hana)
if(NOT hana_POPULATED)
  FetchContent_Populate(hana)

  # Copy an additional/replacement file into the populated source
  file(COPY hana/CMakeLists.txt DESTINATION ${hana_SOURCE_DIR}/)

  add_subdirectory(${hana_SOURCE_DIR} ${hana_BINARY_DIR})
endif()

target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_BOOST_HANA )
#target_link_libraries(feelpp_contrib INTERFACE feelpp_hana )
target_include_directories(feelpp_contrib INTERFACE
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
