#############################################################################
# create generic simple lib
#############################################################################
macro(genLibBase)
  PARSE_ARGUMENTS(FEELMOR_GENLIB_BASE
    "LIB_NAME;LIB_DIR;LIB_DEPENDS;FILES_TO_COPY;FILES_SOURCES;CONFIG_PATH;TARGET_COPY_FILES"
    ""
    ${ARGN}
    )

  set( FEELMOR_GENLIB_APPLICATION_DIR ${FEELMOR_GENLIB_BASE_LIB_DIR} )
  #CAR(APPLICATION_NAME      ${FEELMOR_GENLIB_BASE_DEFAULT_ARGS})
  #CDR(FEELMOR_GENLIB_APPLICATION_DIR       ${FEELMOR_GENLIB_BASE_DEFAULT_ARGS})

  set(LIB_DEPENDS           ${FEELMOR_GENLIB_BASE_LIB_DEPENDS})
  set(LIB_APPLICATION_NAME ${FEELMOR_GENLIB_BASE_LIB_NAME})

  set(CODEGEN_FILES_TO_COPY ${FEELMOR_GENLIB_BASE_FILES_TO_COPY})
  set(CODEGEN_SOURCES       ${FEELMOR_GENLIB_BASE_FILES_SOURCES})

  if ( FEELMOR_GENLIB_BASE_CONFIG_PATH )
    set(FEELMOR_GENLIB_CONFIG_PATH       ${FEELMOR_GENLIB_BASE_CONFIG_PATH})
    get_filename_component(FEELMOR_GENLIB_CONFIG_FILENAME_WE ${FEELMOR_GENLIB_CONFIG_PATH} NAME_WE)
    CONFIGURE_FILE( ${FEELMOR_GENLIB_CONFIG_PATH} ${FEELMOR_GENLIB_APPLICATION_DIR}/${FEELMOR_GENLIB_CONFIG_FILENAME_WE}.h  )
  endif()

  if ( FEELMOR_GENLIB_BASE_TARGET_COPY_FILES )
    set( TARGET_COPY_FILES ${FEELMOR_GENLIB_BASE_TARGET_COPY_FILES} )
  else()
    set( TARGET_COPY_FILES ${LIB_APPLICATION_NAME}_copyfiles)
    add_custom_target(${TARGET_COPY_FILES}  ALL COMMENT "Copying modified files"  )
  endif()

  # lib files
  foreach(filepath ${CODEGEN_FILES_TO_COPY})
    get_filename_component(filename ${filepath} NAME)
    if ( NOT EXISTS ${FEELMOR_GENLIB_APPLICATION_DIR}/${filename} )
      #configure_file( ${filepath} ${FEELMOR_GENLIB_APPLICATION_DIR}/${filename} COPYONLY)
      file(WRITE ${FEELMOR_GENLIB_APPLICATION_DIR}/${filename} "") #write empty file
    endif()
    add_custom_command(TARGET ${TARGET_COPY_FILES} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${filepath} ${FEELMOR_GENLIB_APPLICATION_DIR}/${filename} )
  endforeach()

   # message (CODEGEN_SOURCES   :    ${CODEGEN_SOURCES})
   # message( CODEGEN_FILES_TO_COPY  ${CODEGEN_FILES_TO_COPY})
  # generate library
  add_library(
    ${LIB_APPLICATION_NAME}
    SHARED
    ${CODEGEN_SOURCES}
    )
  add_dependencies(${LIB_APPLICATION_NAME} ${TARGET_COPY_FILES})
  target_include_directories( ${LIB_APPLICATION_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/feelpp/mor  )
  target_link_libraries(${LIB_APPLICATION_NAME} ${LIB_DEPENDS} )
  set_target_properties(${LIB_APPLICATION_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${FEELMOR_GENLIB_APPLICATION_DIR}")
  set_property(TARGET ${LIB_APPLICATION_NAME} PROPERTY MACOSX_RPATH ON)

  # install process
  set_target_properties(${LIB_APPLICATION_NAME} PROPERTIES VERSION ${FEELPP_MOR_SHARED_VERSION} SOVERSION ${FEELPP_MOR_SHARED_SOVERSION})
  INSTALL(TARGETS ${LIB_APPLICATION_NAME} EXPORT feelpp-mor-targets
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/feelpp/mor/ )


endmacro(genLibBase)


#############################################################################
#############################################################################
#############################################################################
macro( genLibThermoElectricNL )
  PARSE_ARGUMENTS(FEELMOR_APP
    "DIM;P_ORDER;GEO_ORDER;"
    ""
    ${ARGN}
    )

  if( NOT ( FEELMOR_APP_DIM OR FEELMOR_APP_P_ORDER OR FEELMOR_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_P_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(THERMOELECTRICNL_DIM ${FEELMOR_APP_DIM})
  set(THERMOELECTRICNL_ORDERPOLY ${FEELMOR_APP_P_ORDER})
  set(THERMOELECTRICNL_ORDERGEO ${FEELMOR_APP_GEO_ORDER})

  set(THERMOELECTRICNL_LIB_VARIANTS ${THERMOELECTRICNL_DIM}dP${THERMOELECTRICNL_ORDERPOLY}G${THERMOELECTRICNL_ORDERGEO} )
  set(THERMOELECTRICNL_LIB_NAME feelpp_mor_thermoelectric_nl_${THERMOELECTRICNL_LIB_VARIANTS})

  if ( NOT TARGET ${THERMOELECTRICNL_LIB_NAME} )
    # configure the lib
    set(THERMOELECTRICNL_LIB_DIR ${FEELPP_MOR_BINARY_DIR}/thermoelectric/${THERMOELECTRICNL_LIB_VARIANTS})
    set(THERMOELECTRICNL_CODEGEN_FILES_TO_COPY
      ${FEELPP_MOR_SOURCE_DIR}/thermoelectric/thermoelectric-nl_create_inst.cpp
      ${FEELPP_MOR_SOURCE_DIR}/thermoelectric/thermoelectric-nl_assemble_inst.cpp
      ${FEELPP_MOR_SOURCE_DIR}/thermoelectric/thermoelectric-nl_others_inst.cpp )
    set(THERMOELECTRICNL_CODEGEN_SOURCES
      ${THERMOELECTRICNL_LIB_DIR}/thermoelectric-nl_create_inst.cpp
      ${THERMOELECTRICNL_LIB_DIR}/thermoelectric-nl_assemble_inst.cpp
      ${THERMOELECTRICNL_LIB_DIR}/thermoelectric-nl_others_inst.cpp )
    set(THERMOELECTRICNL_LIB_DEPENDS Feelpp::feelpp )
    # generate the lib target
    genLibBase(
      LIB_NAME ${THERMOELECTRICNL_LIB_NAME}
      LIB_DIR ${THERMOELECTRICNL_LIB_DIR}
      LIB_DEPENDS ${THERMOELECTRICNL_LIB_DEPENDS}
      FILES_TO_COPY ${THERMOELECTRICNL_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${THERMOELECTRICNL_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MOR_SOURCE_DIR}/thermoelectric/thermoelectric-nl_config.h.in
      )
  endif()
endmacro(genLibThermoElectricNL)
