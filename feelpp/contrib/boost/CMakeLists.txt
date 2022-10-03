include(FetchContent)
set(FETCHCONTENT_SOURCE_DIR_HANA ${CMAKE_CURRENT_SOURCE_DIR}/hana)
FetchContent_Declare(hana GIT_REPOSITORY https://github.com/feelpp/hana.git GIT_TAG 4e402573f45bd78c5645359e07ad539b3c2b07a1) # the HEAD of feelpp branch
FetchContent_MakeAvailable(hana)

target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_BOOST_HANA )
#target_link_libraries(feelpp_contrib INTERFACE boost_hana )
message(STATUS "-- Boost Hana: ${hana_SOURCE_DIR}")
target_link_libraries(feelpp_contrib INTERFACE boost_hana )
target_include_directories(feelpp_contrib BEFORE INTERFACE
  $<BUILD_INTERFACE:${hana_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/feelpp/> )
install(TARGETS boost_hana
        EXPORT feelpp-contrib-export-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY hana/include/boost
        DESTINATION include/feelpp
        FILES_MATCHING PATTERN "*.hpp")
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
