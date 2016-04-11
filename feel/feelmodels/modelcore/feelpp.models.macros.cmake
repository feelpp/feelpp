#############################################################################
#############################################################################
macro(resetToZeroAllPhysicalVariables)
  SET(SOLIDMECHANICS 0 )
  SET(SOLIDMECHANICS0 0 )
  SET(SOLIDMECHANICS1 0 )
  SET(SOLIDMECHANICS2 0 )
  SET(FLUIDMECHANICS 0 )
  SET(FLUIDMECHANICS0 0 )
  SET(FLUIDMECHANICS1 0 )
  SET(FLUIDMECHANICS2 0 )
  SET(THERMODYNAMICS 0 )
  SET(THERMODYNAMICS0 0 )
  SET(THERMODYNAMICS1 0 )
  SET(THERMODYNAMICS2 0 )
endmacro(resetToZeroAllPhysicalVariables)
#############################################################################
# create generic simple lib : ${PREFIX_LIB}${APPLICATION_NAME}
#############################################################################
macro(genLibBase)
  PARSE_ARGUMENTS(FEELMODELS_GENLIB_BASE
    "LIB_NAME;LIB_DIR;MARKERS;DESC;GEO;LIB_DEPENDS;PREFIX_INCLUDE_USERCONFIG;FILES_TO_COPY;FILES_SOURCES;CONFIG_PATH;ADD_CMAKE_INSTALL"
    ""
    ${ARGN}
    )
  set( FEELMODELS_GENLIB_APPLICATION_DIR ${FEELMODELS_GENLIB_BASE_LIB_DIR} )
  #CAR(APPLICATION_NAME      ${FEELMODELS_GENLIB_BASE_DEFAULT_ARGS})
  #CDR(FEELMODELS_GENLIB_APPLICATION_DIR       ${FEELMODELS_GENLIB_BASE_DEFAULT_ARGS})

  set(LIB_DEPENDS           ${FEELMODELS_GENLIB_BASE_LIB_DEPENDS})
  set(PREFIX_FILES_TO_COPY  ${FEELMODELS_GENLIB_BASE_PREFIX_INCLUDE_USERCONFIG})
  set(LIB_APPLICATION_NAME ${FEELMODELS_GENLIB_BASE_LIB_NAME})

  set(CODEGEN_FILES_TO_COPY ${FEELMODELS_GENLIB_BASE_FILES_TO_COPY})
  set(CODEGEN_SOURCES       ${FEELMODELS_GENLIB_BASE_FILES_SOURCES})

  if ( FEELMODELS_GENLIB_BASE_CONFIG_PATH )
    set(FEELMODELS_GENLIB_CONFIG_PATH       ${FEELMODELS_GENLIB_BASE_CONFIG_PATH})
    get_filename_component(FEELMODELS_GENLIB_CONFIG_FILENAME_WE ${FEELMODELS_GENLIB_CONFIG_PATH} NAME_WE)
    CONFIGURE_FILE( ${FEELMODELS_GENLIB_CONFIG_PATH} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${FEELMODELS_GENLIB_CONFIG_FILENAME_WE}.h  )
  endif()

  add_custom_target(codegen_${LIB_APPLICATION_NAME}  ALL COMMENT "Copying modified files"  )

  if ( FEELMODELS_GENLIB_BASE_MARKERS )
    # bcmarker
    set(BCMARKER_FILE ${FEELMODELS_GENLIB_BASE_MARKERS})
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/bcmarker.cpp )
      #configure_file( ${BCMARKER_FILE} ${FEELMODELS_GENLIB_APPLICATION_DIR}/bcmarker.cpp COPYONLY)
      file(WRITE ${FEELMODELS_GENLIB_APPLICATION_DIR}/bcmarker.cpp "") #write empty file
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${BCMARKER_FILE}  ${FEELMODELS_GENLIB_APPLICATION_DIR}/bcmarker.cpp )
    # bctool
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/bctool.hpp )
      #configure_file( ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/bctool.hpp ${FEELMODELS_GENLIB_APPLICATION_DIR}/bctool.hpp COPYONLY)
      file(WRITE ${FEELMODELS_GENLIB_APPLICATION_DIR}/bctool.hpp "") #write empty file
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/bctool.hpp ${FEELMODELS_GENLIB_APPLICATION_DIR}/bctool.hpp )
  endif()
  # .bc
  if ( FEELMODELS_GENLIB_BASE_DESC )
    set(BCDESC_FILE ${FEELMODELS_GENLIB_BASE_DESC})
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.bc )
      #configure_file( ${BCDESC_FILE} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.bc COPYONLY)
      file(WRITE ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.bc "") #write empty file
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${BCDESC_FILE}  ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.bc )
  endif()
  # .mesh
  if ( FEELMODELS_GENLIB_BASE_GEO )
    set(MESH_FILE ${FEELMODELS_GENLIB_BASE_GEO})
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.mesh )
      #configure_file( ${MESH_FILE} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.mesh COPYONLY)
      file(WRITE ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.mesh "") #write empty file
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${MESH_FILE}  ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.mesh )
  endif()
  # lib files
  foreach(filepath ${CODEGEN_FILES_TO_COPY})
    get_filename_component(filename ${filepath} NAME)
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} )
      #configure_file( ${filepath} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} COPYONLY)
      file(WRITE ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} "") #write empty file
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${filepath} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} )
  endforeach()

   # message (CODEGEN_SOURCES   :    ${CODEGEN_SOURCES})
   # message( CODEGEN_FILES_TO_COPY  ${CODEGEN_FILES_TO_COPY})
  # generate library
  add_library(
    ${LIB_APPLICATION_NAME}
    SHARED
    ${CODEGEN_SOURCES}
    )
  add_dependencies(${LIB_APPLICATION_NAME} codegen_${LIB_APPLICATION_NAME})
  target_link_libraries(${LIB_APPLICATION_NAME} ${LIB_DEPENDS} )
  set_target_properties(${LIB_APPLICATION_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${FEELMODELS_GENLIB_APPLICATION_DIR}")

  # install process
  if ( FEELMODELS_GENLIB_BASE_ADD_CMAKE_INSTALL )
    INSTALL(TARGETS ${LIB_APPLICATION_NAME} DESTINATION lib/ COMPONENT LibsFeelppModels-${LIB_APPLICATION_NAME})
    add_custom_target(install-${LIBBASE_NAME}
      DEPENDS ${LIBBASE_NAME}
      COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=LibsFeelppModels-${LIB_APPLICATION_NAME}
      -P "${CMAKE_BINARY_DIR}/cmake_install.cmake"
      )
  endif()

endmacro(genLibBase)
#############################################################################

#############################################################################
macro(genExecutableBase)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "SRC;LIB_DEPENDS;CONFIG_PATH;APPLICATION_DIR"
    ""
    ${ARGN}
    )
  CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})
  #CDR(APPLICATION_DIR ${FEELMODELS_APP_DEFAULT_ARGS})
  set(APPLICATION_DIR ${FEELMODELS_APP_APPLICATION_DIR})
  set(MAIN_FILE ${FEELMODELS_APP_SRC})
  set(LIB_DEPENDS ${FEELMODELS_APP_LIB_DEPENDS})

  if ( FEELMODELS_APP_CONFIG_PATH )
    set(GENEXECBASE_CONFIG_PATH       ${FEELMODELS_APP_CONFIG_PATH})
    get_filename_component(GENEXECBASE_CONFIG_FILENAME_WE ${GENEXECBASE_CONFIG_PATH} NAME_WE)
    CONFIGURE_FILE( ${GENEXECBASE_CONFIG_PATH} ${APPLICATION_DIR}/${GENEXECBASE_CONFIG_FILENAME_WE}.h  )
  endif()

  if ( NOT EXISTS ${APPLICATION_DIR}/applimanagement.hpp )
    foreach(filename applimanagement.hpp)
      #configure_file( ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/${filename} ${APPLICATION_DIR}/${filename} COPYONLY)
      file(WRITE ${APPLICATION_DIR}/${filename} "") #write empty file
    endforeach()
  endif()
  add_custom_target(codegen_env_${APPLICATION_NAME} ALL COMMENT "Copying modified files"  )
  foreach(filename applimanagement.hpp)
    add_custom_command(TARGET codegen_env_${APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/${filename} ${APPLICATION_DIR}/${filename} )
  endforeach()

  if ( NOT EXISTS ${APPLICATION_DIR}/${MAIN_FILE} )
    foreach(filename ${MAIN_FILE})
      #configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${MAIN_FILE} ${APPLICATION_DIR}/${MAIN_FILE} COPYONLY)
      file(WRITE ${APPLICATION_DIR}/${filename} "") #write empty file
    endforeach()
  endif()
  add_custom_target(codegen_src_${APPLICATION_NAME} ALL COMMENT "Copying modified files"  )
  foreach(filename ${MAIN_FILE})
    add_custom_command(TARGET codegen_src_${APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${CMAKE_CURRENT_SOURCE_DIR}/${MAIN_FILE} ${APPLICATION_DIR}/${MAIN_FILE} )
  endforeach()


  add_executable( ${APPLICATION_NAME} ${APPLICATION_DIR}/${MAIN_FILE} )
  add_dependencies(${APPLICATION_NAME} codegen_env_${APPLICATION_NAME})
  add_dependencies(${APPLICATION_NAME} codegen_src_${APPLICATION_NAME})

  target_link_libraries(${APPLICATION_NAME} ${LIB_DEPENDS} )

endmacro(genExecutableBase)
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibThermoDynamics )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;T_ORDER;GEO_ORDER;BC_MARKERS;BC_DESC;GEO_DESC;LIB_NAME;LIB_DIR;"
    "NO_UPDATE_MODEL_DEF;ADD_CMAKE_INSTALL"
    ${ARGN}
    )

  if ( NOT FEELMODELS_APP_NO_UPDATE_MODEL_DEF )
    resetToZeroAllPhysicalVariables()
    SET(THERMODYNAMICS 1 )
  endif()

  #CAR(LIB_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(THERMODYNAMICS_DIM ${FEELMODELS_APP_DIM})
  set(THERMODYNAMICS_ORDERPOLY ${FEELMODELS_APP_T_ORDER})
  set(THERMODYNAMICS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  if (0)
    MESSAGE("*** Arguments for fluid application ${LIB_NAME}")
    MESSAGE("*** DIM ${THERMODYNAMICS_DIM}")
    MESSAGE("*** ORDERPOLY ${THERMODYNAMICS_ORDERPOLY}")
    MESSAGE("*** ORDERGEO ${THERMODYNAMICS_ORDERGEO}")
  endif()

  set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${THERMODYNAMICS_DIM}dP${THERMODYNAMICS_ORDERPOLY}G${THERMODYNAMICS_ORDERGEO} )
  set(FEELMODELS_MODEL_SPECIFIC_NAME thermodyn${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})
  #set(FEELMODELS_MODEL_SPECIFIC_NAME ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})
  set(LIBBASE_DIR ${FEELPP_MODELS_BINARY_DIR}/thermodyn/${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX} )
  set(LIBBASE_CHECK_PATH ${FEELPP_MODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
  set(LIBBASE_NAME feelpp_model_${FEELMODELS_MODEL_SPECIFIC_NAME})

  if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )

    #write empty file in orter to check if this lib has already define
    file(WRITE ${LIBBASE_CHECK_PATH} "")

    # configure libmodelbase
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodynbase_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodynamics_inst.cpp )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/thermodynbase_inst.cpp
      ${LIBBASE_DIR}/thermodynamics_inst.cpp )
    set(LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 

    if ( FEELMODELS_APP_ADD_CMAKE_INSTALL )
      set( LIBBASE_ADD_CMAKE_INSTALL 1 )
    else()
      set( LIBBASE_ADD_CMAKE_INSTALL 0 )
    endif()

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodynconfig.h.in
      ADD_CMAKE_INSTALL ${LIBBASE_ADD_CMAKE_INSTALL}
      )

  endif()

  ############################
  # build lib with bctool
  ############################
  if ( FEELMODELS_APP_BC_DESC )

    get_filename_component( BCDESC_FILE ${FEELMODELS_APP_BC_DESC} ABSOLUTE )

    if ( FEELMODELS_APP_BC_MARKERS )
      get_filename_component( BCMARKER_FILE ${FEELMODELS_APP_BC_MARKERS} ABSOLUTE )
    else()
      set(BCMARKER_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.bcmarkers )
    endif()

    if ( FEELMODELS_APP_GEO_DESC )
      get_filename_component( GEO_FILE ${FEELMODELS_APP_GEO_DESC} ABSOLUTE )
    else()
      set(GEO_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.mesh )
    endif()

    if ( NOT ( FEELMODELS_APP_LIB_NAME  ) )
      message(FATAL_ERROR "miss argument! FEELMODELS_APP_LIB_NAME")
    endif()
    set(THE_LIB_NAME ${FEELMODELS_APP_LIB_NAME})

    if ( FEELMODELS_APP_LIB_DIR )
      get_filename_component( THE_LIB_DIR ${FEELMODELS_APP_LIB_DIR} ABSOLUTE )
    else()
      set(THE_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR})
    endif()


    set(THERMODYNAMICS_INCLUDE_LIBBASE_HEADER 1)

    set(PREFIX_FILES_TO_COPY thermodyn)
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/codegen_thermodyn.hpp
      ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/codegen_thermodyn.cpp )
    set(CODEGEN_SOURCES
      ${THE_LIB_DIR}/codegen_thermodyn.cpp )
    set(LIB_DEPENDS ${LIBBASE_NAME} feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 


    genLibBase(
      LIB_NAME ${THE_LIB_NAME}
      LIB_DIR ${THE_LIB_DIR}
      MARKERS ${BCMARKER_FILE}
      DESC ${BCDESC_FILE}
      GEO ${GEO_FILE}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/feelmodelscoreconfig.h.in
      )
    unset(THERMODYNAMICS_INCLUDE_LIBBASE_HEADER)

  endif( FEELMODELS_APP_BC_DESC )

endmacro(genLibThermoDynamics)
#############################################################################
#############################################################################
macro(genExecutableThermoDynamics APPLICATION_NAME FEELMODELS_DIM FEELMODELS_ORDERPOLY FEELMODELS_ORDERGEO BCMARKER_FILE BCDESC_FILE GEO_FILE MAIN_FILE)

  resetToZeroAllPhysicalVariables()
  SET(THERMODYNAMICS 1 )

  set(THE_APPLICATION_DIR ${CMAKE_CURRENT_BINARY_DIR}/codeGen_${APPLICATION_NAME})
  #set(APPLICATION_DIR codeGen_${APPLICATION_NAME})
  #genLibThermoDynamics(${APPLICATION_NAME} ${APPLICATION_DIR}/thermodyn ${FEELMODELS_DIM} ${FEELMODELS_ORDERPOLY} ${FEELMODELS_ORDERGEO} ${BCMARKER_FILE} ${BCDESC_FILE} ${MESH_FILE} )
  genLibThermoDynamics(
    LIB_NAME   feelpp_model_thermo${APPLICATION_NAME}
    LIB_DIR    ${THE_APPLICATION_DIR}/thermodyn
    DIM        ${FEELMODELS_DIM}
    T_ORDER    ${FEELMODELS_ORDERPOLY}
    GEO_ORDER  ${FEELMODELS_ORDERGEO}
    BC_DESC    ${BCDESC_FILE} 
    BC_MARKERS ${BCMARKER_FILE}
    GEO_DESC   ${GEO_FILE}
    )

  set(LIB_DEPENDS feelpp_model_thermo${APPLICATION_NAME} feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) #fsimesh 
  genExecutableBase(
    ${APPLICATION_NAME}
    APPLICATION_DIR ${THE_APPLICATION_DIR}
    SRC ${MAIN_FILE}
    LIB_DEPENDS ${LIB_DEPENDS}
    CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/feelmodelscoreconfig.h.in
    )

endmacro(genExecutableThermoDynamics)
#############################################################################
#############################################################################
macro(feelpp_add_thermo_application)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;T_ORDER;GEO_ORDER;BC_MARKERS;BC_DESC;GEO_DESC;SRC"
    ""
    ${ARGN}
    )

  CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
  set(FEELMODELS_ORDERPOLY ${FEELMODELS_APP_T_ORDER})
  set(FEELMODELS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  if ( FEELMODELS_APP_BC_MARKERS )
    set(BCMARKER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_MARKERS})
  else()
    set(BCMARKER_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.bcmarkers )
  endif()

  if ( FEELMODELS_APP_GEO_DESC )
    set(GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GEO_DESC})
  else()
    set(GEO_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.mesh )
  endif()

  set(DESC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_DESC})
  set(MAIN_FILE ${FEELMODELS_APP_SRC})

  if (0)
    MESSAGE("*** Arguments for fluid application ${APPLICATION_NAME}")
    MESSAGE("*** DIM ${FEELMODELS_DIM}")
    MESSAGE("*** ORDERPOLY ${FEELMODELS_ORDERPOLY}")
    MESSAGE("*** ORDERGEO ${FEELMODELS_ORDERGEO}")
    MESSAGE("*** MARKERS ${BCMARKER_FILE}")
    MESSAGE("*** DESC ${DESC_FILE}")
    MESSAGE("*** GEO ${GEO_FILE}")
    MESSAGE("*** SRC ${MAIN_FILE}")
  endif()

  genExecutableThermoDynamics(${APPLICATION_NAME} ${FEELMODELS_DIM} ${FEELMODELS_ORDERPOLY} ${FEELMODELS_ORDERGEO} ${BCMARKER_FILE} ${DESC_FILE} ${GEO_FILE} ${MAIN_FILE} )

endmacro(feelpp_add_thermo_application)


#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro( genLibSolidMechanics )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;DISP_ORDER;GEO_ORDER;GEO_DESC;BC_MARKERS;BC_DESC;DENSITY_COEFFLAME_TYPE;LIB_NAME;LIB_DIR;ADD_CMAKE_INSTALL_BIS"
    "NO_UPDATE_MODEL_DEF;ADD_CMAKE_INSTALL"
    ${ARGN}
    )

  if ( NOT FEELMODELS_APP_NO_UPDATE_MODEL_DEF )
    resetToZeroAllPhysicalVariables()
    SET(SOLIDMECHANICS 1 )
  endif()

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_DISP_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_DISP_ORDER OR FEELMODELS_APP_GEO_ORDER")
  endif()

  set(SOLIDMECHANICS_DIM ${FEELMODELS_APP_DIM})
  set(SOLIDMECHANICS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
  set(SOLIDMECHANICS_ORDER_DISPLACEMENT ${FEELMODELS_APP_DISP_ORDER})

  if (FEELMODELS_APP_DENSITY_COEFFLAME_TYPE)
    if ("${FEELMODELS_APP_DENSITY_COEFFLAME_TYPE}" STREQUAL "P0c" )
      set(FEELMODELS_DENSITY_COEFFLAME_TAG DCLP0c)
      set(SOLIDMECHANICS_USE_CST_DENSITY_COEFFLAME 1)
    elseif ("${FEELMODELS_APP_DENSITY_COEFFLAME_TYPE}" STREQUAL "P0d" )
      unset(FEELMODELS_DENSITY_COEFFLAME_TAG) # default tag P0d
      set(SOLIDMECHANICS_USE_CST_DENSITY_COEFFLAME 0)
    else()
      message(FATAL_ERROR "DENSITY_COEFFLAME_TYPE : ${FEELMODELS_APP_DENSITY_COEFFLAME_TYPE} is not valid! It must be P0c or P0d")
    endif()
  else()
    # default value
    unset(FEELMODELS_DENSITY_COEFFLAME_TAG) # default tag P0d
    set(SOLIDMECHANICS_USE_CST_DENSITY_COEFFLAME 0)
  endif()

  set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${SOLIDMECHANICS_DIM}dP${SOLIDMECHANICS_ORDER_DISPLACEMENT}G${SOLIDMECHANICS_ORDERGEO}${FEELMODELS_DENSITY_COEFFLAME_TAG})
  set(FEELMODELS_MODEL_SPECIFIC_NAME solidmec${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})

  set(LIBBASE_DIR ${FEELPP_MODELS_BINARY_DIR}/solid/${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX} )
  set(LIBBASE_CHECK_PATH ${FEELPP_MODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
  set(LIBBASE_NAME feelpp_model_${FEELMODELS_MODEL_SPECIFIC_NAME})

  if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )
    #write empty file in orter to check if this lib has already define
    file(WRITE ${LIBBASE_CHECK_PATH} "")

    # configure libmodelbase
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmecbasecreate_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmecbaseothers_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmecbaseupdatelinear_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmecbaseupdatelinear1dreduced_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmecbaseupdatejacobian_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmecbaseupdateresidual_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmechanics_inst.cpp
      )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/solidmecbasecreate_inst.cpp
      ${LIBBASE_DIR}/solidmecbaseothers_inst.cpp
      ${LIBBASE_DIR}/solidmecbaseupdatelinear_inst.cpp
      ${LIBBASE_DIR}/solidmecbaseupdatelinear1dreduced_inst.cpp
      ${LIBBASE_DIR}/solidmecbaseupdatejacobian_inst.cpp
      ${LIBBASE_DIR}/solidmecbaseupdateresidual_inst.cpp
      ${LIBBASE_DIR}/solidmechanics_inst.cpp
      )
    set(LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 

    if ( FEELMODELS_APP_ADD_CMAKE_INSTALL )
      set( LIBBASE_ADD_CMAKE_INSTALL 1 )
    else()
      set( LIBBASE_ADD_CMAKE_INSTALL 0 )
    endif()
    # overwrite options
    if ( FEELMODELS_APP_ADD_CMAKE_INSTALL_BIS )
      set( LIBBASE_ADD_CMAKE_INSTALL ${FEELMODELS_APP_ADD_CMAKE_INSTALL_BIS} )
    endif()

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmecconfig.h.in
      ADD_CMAKE_INSTALL ${LIBBASE_ADD_CMAKE_INSTALL}
      )

  endif()

  ############################
  # build lib with bctool
  ############################
  if ( FEELMODELS_APP_BC_DESC )

    get_filename_component( BCDESC_FILE ${FEELMODELS_APP_BC_DESC} ABSOLUTE )

    if ( FEELMODELS_APP_BC_MARKERS )
      get_filename_component( BCMARKER_FILE ${FEELMODELS_APP_BC_MARKERS} ABSOLUTE )
    else()
      set(BCMARKER_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.bcmarkers )
    endif()

    if ( FEELMODELS_APP_GEO_DESC )
      get_filename_component( GEO_FILE ${FEELMODELS_APP_GEO_DESC} ABSOLUTE )
    else()
      set(GEO_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.mesh )
    endif()

    if ( NOT ( FEELMODELS_APP_LIB_NAME  ) )
      message(FATAL_ERROR "miss argument! FEELMODELS_APP_LIB_NAME")
    endif()
    set(THE_LIB_NAME ${FEELMODELS_APP_LIB_NAME})

    if ( FEELMODELS_APP_LIB_DIR )
      get_filename_component( THE_LIB_DIR ${FEELMODELS_APP_LIB_DIR} ABSOLUTE )
    else()
      set(THE_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR})
    endif()


    set(SOLIDMECHANICS_INCLUDE_LIBBASE_HEADER 1)

    set(PREFIX_FILES_TO_COPY solid)
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/solid/codegen_solidmec.hpp
      ${FEELPP_MODELS_SOURCE_DIR}/solid/codegen_solidmec.cpp )
    set(CODEGEN_SOURCES
      ${THE_LIB_DIR}/codegen_solidmec.cpp )
    set(LIB_DEPENDS ${LIBBASE_NAME} feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 


    genLibBase(
      LIB_NAME ${THE_LIB_NAME}
      LIB_DIR ${THE_LIB_DIR}
      MARKERS ${BCMARKER_FILE}
      DESC ${BCDESC_FILE}
      GEO ${GEO_FILE}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/feelmodelscoreconfig.h.in
      )
    unset(SOLIDMECHANICS_INCLUDE_LIBBASE_HEADER)

  endif( FEELMODELS_APP_BC_DESC )

endmacro( genLibSolidMechanics )
#############################################################################
#############################################################################
macro(genExecutableSolidMechanics APPLICATION_NAME
    FEELMODELS_DIM FEELMODELS_ORDERGEO FEELMODELS_ORDER_DISP FEELMODELS_DENSITY_COEFFLAME_TYPE
    BCMARKER_FILE BCDESC_FILE GEO_FILE MAIN_FILE ADDITIONAL_LIB_DEPEND )

  resetToZeroAllPhysicalVariables()
  SET(SOLIDMECHANICS 1 )

  set(THE_APPLICATION_DIR ${CMAKE_CURRENT_BINARY_DIR}/codeGen_${APPLICATION_NAME})
  genLibSolidMechanics(
    LIB_NAME   feelpp_model_solid_${APPLICATION_NAME}
    LIB_DIR    ${THE_APPLICATION_DIR}/solid
    DIM        ${FEELMODELS_DIM}
    GEO_ORDER  ${FEELMODELS_ORDERGEO}
    DISP_ORDER ${FEELMODELS_ORDER_DISP}
    DENSITY_COEFFLAME_TYPE ${FEELMODELS_DENSITY_COEFFLAME_TYPE}
    BC_DESC    ${BCDESC_FILE} 
    BC_MARKERS ${BCMARKER_FILE}
    GEO_DESC   ${GEO_FILE}
    )

  set(LIB_DEPENDS feelpp_model_solid_${APPLICATION_NAME} feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ${ADDITIONAL_LIB_DEPEND} ) 
  genExecutableBase(
    ${APPLICATION_NAME}
    APPLICATION_DIR ${THE_APPLICATION_DIR}
    SRC ${MAIN_FILE}
    LIB_DEPENDS ${LIB_DEPENDS}
    CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/feelmodelscoreconfig.h.in
    )

endmacro(genExecutableSolidMechanics)
#############################################################################
#############################################################################
macro(feelpp_add_solid_application)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;DISP_ORDER;GEO_ORDER;GEO_DESC;BC_MARKERS;BC_DESC;SRC;DENSITY_COEFFLAME_TYPE;ADDITIONAL_LIB_DEPENDS"
    ""
    ${ARGN}
    )

  CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_DISP_ORDER OR FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_DISP_ORDER OR FEELMODELS_APP_GEO_ORDER")
  endif()
  set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
  set(FEELMODELS_ORDER_DISP ${FEELMODELS_APP_DISP_ORDER})
  set(FEELMODELS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  if ( FEELMODELS_APP_BC_MARKERS )
    set(BCMARKER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_MARKERS})
  else()
    set(BCMARKER_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.bcmarkers )
  endif()

  if ( FEELMODELS_APP_GEO_DESC )
    set(GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GEO_DESC})
  else()
    set(GEO_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.mesh )
  endif()

  set(DESC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_DESC})
  set(MAIN_FILE ${FEELMODELS_APP_SRC})


  if (FEELMODELS_APP_DENSITY_COEFFLAME_TYPE)
    set( FEELMODELS_DENSITY_COEFFLAME_TYPE ${FEELMODELS_APP_DENSITY_COEFFLAME_TYPE})
  else()
    # default value
    set(FEELMODELS_DENSITY_COEFFLAME_TYPE "P0d")
  endif()

  if ( FEELMODELS_APP_ADDITIONAL_LIB_DEPENDS )
    set( ADDITIONAL_LIB_DEPENDS ${FEELMODELS_APP_ADDITIONAL_LIB_DEPENDS} )
  else()
    # put feelpp ( duplication but allow to give an arg )
    set( ADDITIONAL_LIB_DEPENDS ${FEELPP_LIBRARY} )
  endif()

  if (0)
    MESSAGE("*** Arguments for solid application ${APPLICATION_NAME}")
    MESSAGE("*** DIM ${FEELMODELS_DIM}")
    MESSAGE("*** DISP_ORDER ${FEELMODELS_ORDER_DISP}")
    MESSAGE("*** ORDERGEO ${FEELMODELS_ORDERGEO}")
    MESSAGE("*** GEO_DESC ${GEO_FILE}")
    MESSAGE("*** BC_MARKERS ${BCMARKER_FILE}")
    MESSAGE("*** BC_DESC ${DESC_FILE}")
    MESSAGE("*** SRC ${MAIN_FILE}")
    MESSAGE("*** DENSITY_COEFFLAME_TYPE ${FEELMODELS_DENSITY_COEFFLAME_TYPE}")
  endif()

  genExecutableSolidMechanics(${APPLICATION_NAME}
    ${FEELMODELS_DIM} ${FEELMODELS_ORDERGEO} ${FEELMODELS_ORDER_DISP} ${FEELMODELS_DENSITY_COEFFLAME_TYPE}
    ${BCMARKER_FILE} ${DESC_FILE} ${GEO_FILE} ${MAIN_FILE} ${ADDITIONAL_LIB_DEPENDS} )

endmacro(feelpp_add_solid_application)
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
macro(genLibFluidMechanics)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;U_ORDER;P_ORDER;P_CONTINUITY;GEO_ORDER;GEO_DESC;BC_MARKERS;BC_DESC;DENSITY_VISCOSITY_CONTINUITY;DENSITY_VISCOSITY_ORDER;USE_PERIODICITY_BIS;LIB_NAME;LIB_DIR;ADD_CMAKE_INSTALL_BIS"
    "USE_PERIODICITY;NO_UPDATE_MODEL_DEF;ADD_CMAKE_INSTALL"
    ${ARGN}
    )

  if ( NOT FEELMODELS_APP_NO_UPDATE_MODEL_DEF )
    resetToZeroAllPhysicalVariables()
    SET(FLUIDMECHANICS 1 )
  endif()

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_GEO_ORDER OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_P_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_GEO_ORDER OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_P_ORDER")
  endif()

  set(FLUIDMECHANICS_DIM ${FEELMODELS_APP_DIM})
  set(FLUIDMECHANICS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
  set(FLUIDMECHANICS_ORDER_VELOCITY ${FEELMODELS_APP_U_ORDER})
  set(FLUIDMECHANICS_ORDER_PRESSURE ${FEELMODELS_APP_P_ORDER})

  #######################################################
  if ( FEELMODELS_APP_P_CONTINUITY )
    if ("${FEELMODELS_APP_P_CONTINUITY}" STREQUAL "Continuous" )
      set(FLUIDMECHANICS_PRESSURE_IS_CONTINUOUS 1)
      unset(FLUIDMECHANICS_PRESSURE_CONTINUITY_TAG) # default tag continuous
    elseif("${FEELMODELS_APP_P_CONTINUITY}" STREQUAL "Discontinuous" )
      set(FLUIDMECHANICS_PRESSURE_IS_CONTINUOUS 0)
      set(FLUIDMECHANICS_PRESSURE_CONTINUITY_TAG  d)
    else()
      message(FATAL_ERROR "P_CONTINUITY ${FEELMODELS_APP_P_CONTINUITY} : is not valid! It must be Continuous or Discontinuous")
    endif()
  else()
    # default value
    set(FLUIDMECHANICS_PRESSURE_IS_CONTINUOUS 1)
    unset(FLUIDMECHANICS_PRESSURE_CONTINUITY_TAG) # default tag continuous
  endif()
  #######################################################
  if (FEELMODELS_APP_DENSITY_VISCOSITY_ORDER)
    set(FLUIDMECHANICS_ORDER_DENSITY_VISCOSITY ${FEELMODELS_APP_DENSITY_VISCOSITY_ORDER} )
  else()
    # default value
    set(FLUIDMECHANICS_ORDER_DENSITY_VISCOSITY 0)
  endif()
  #######################################################
  if (FEELMODELS_APP_DENSITY_VISCOSITY_CONTINUITY)
    if ("${FEELMODELS_APP_DENSITY_VISCOSITY_CONTINUITY}" STREQUAL "Continuous" )
      set(FLUIDMECHANICS_USE_CONTINUOUS_DENSITY_VISCOSITY 1)
      set(FLUIDMECHANICS_DENSITY_VISCOSITY_TAG DVP${FLUIDMECHANICS_ORDER_DENSITY_VISCOSITY}c)
    elseif ("${FEELMODELS_APP_DENSITY_VISCOSITY_CONTINUITY}" STREQUAL "Discontinuous" )
      set(FLUIDMECHANICS_USE_CONTINUOUS_DENSITY_VISCOSITY 0)
      if ( "${FLUIDMECHANICS_ORDER_DENSITY_VISCOSITY}" STREQUAL "0" )
        unset(FLUIDMECHANICS_DENSITY_VISCOSITY_TAG) # default value P0d
      else()
        set(FLUIDMECHANICS_DENSITY_VISCOSITY_TAG DVP${FLUIDMECHANICS_ORDER_DENSITY_VISCOSITY}d)
      endif()
    else()
      message(FATAL_ERROR "DENSITY_VISCOSITY_CONTINUITY ${FEELMODELS_APP_DENSITY_VISCOSITY}_CONTINUITY : is not valid! It must be Continuous or Discontinuous")
    endif()
  else()
    # default value
    set(FLUIDMECHANICS_USE_CONTINUOUS_DENSITY_VISCOSITY 0)
    unset(FLUIDMECHANICS_DENSITY_VISCOSITY_TAG) # default value P0d
  endif()
  #######################################################
   if ( FEELMODELS_APP_USE_PERIODICITY )
     set(FLUIDMECHANICS_USE_PERIODICITY 1)
     set(FLUIDMECHANICS_USE_PERIODICITY_TAG Periodic)
   else()
     set(FLUIDMECHANICS_USE_PERIODICITY 0)
     unset(FLUIDMECHANICS_USE_PERIODICITY_TAG)
   endif()
   # allow to owerwrite USE_PERIODICITY (used in genExecutable )
   if( FEELMODELS_APP_USE_PERIODICITY_BIS )
     set(FLUIDMECHANICS_USE_PERIODICITY ${FEELMODELS_APP_USE_PERIODICITY_BIS} )
     if ( FLUIDMECHANICS_USE_PERIODICITY )
       set(FLUIDMECHANICS_USE_PERIODICITY_TAG Periodic)
     else()
       unset(FLUIDMECHANICS_USE_PERIODICITY_TAG)
     endif()
   endif()
  #######################################################


  set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${FLUIDMECHANICS_DIM}dP${FLUIDMECHANICS_ORDER_VELOCITY}P${FLUIDMECHANICS_ORDER_PRESSURE}${FLUIDMECHANICS_PRESSURE_CONTINUITY_TAG}G${FLUIDMECHANICS_ORDERGEO}${FLUIDMECHANICS_DENSITY_VISCOSITY_TAG}${FLUIDMECHANICS_USE_PERIODICITY_TAG})
  set(FEELMODELS_MODEL_SPECIFIC_NAME fluidmec${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})

  set(LIBBASE_DIR ${FEELPP_MODELS_BINARY_DIR}/fluid/${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX} )
  set(LIBBASE_CHECK_PATH ${FEELPP_MODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
  set(LIBBASE_NAME feelpp_model_${FEELMODELS_MODEL_SPECIFIC_NAME})

  if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )
    #write empty file in orter to check if this lib has already define
    file(WRITE ${LIBBASE_CHECK_PATH} "")

    # configure libmodelbase
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbasecreate_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbaseothers_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbaseupdatelinear_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbaseupdatelinearweakbc_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbaseupdatejacobian_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbaseupdatejacobianstresstensorlaw_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbaseupdateresidual_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbaseupdateresidualstresstensorlaw_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecbaseupdatestabilisation_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanics_inst.cpp
      )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/fluidmecbasecreate_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseothers_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseupdatelinear_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseupdatelinearweakbc_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseupdatejacobian_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseupdatejacobianstresstensorlaw_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseupdateresidual_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseupdateresidualstresstensorlaw_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseupdateresidualstresstensorlaw_inst.cpp
      ${LIBBASE_DIR}/fluidmecbaseupdatestabilisation_inst.cpp
      ${LIBBASE_DIR}/fluidmechanics_inst.cpp
      )
    set(LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} )
    if (FEELPP_MODELS_ENABLE_MESHALE )
      set(LIB_DEPENDS feelpp_modelmeshale ${LIB_DEPENDS})
    endif()
    # thermodynamcis depend
    set(LIB_DEPENDS feelpp_model_thermodyn${FLUIDMECHANICS_DIM}dP${FLUIDMECHANICS_ORDERGEO}G${FLUIDMECHANICS_ORDERGEO} ${LIB_DEPENDS})


    if ( FEELMODELS_APP_ADD_CMAKE_INSTALL )
      set( LIBBASE_ADD_CMAKE_INSTALL 1 )
    else()
      set( LIBBASE_ADD_CMAKE_INSTALL 0 )
    endif()
    # overwrite options
    if ( FEELMODELS_APP_ADD_CMAKE_INSTALL_BIS )
      set( LIBBASE_ADD_CMAKE_INSTALL ${FEELMODELS_APP_ADD_CMAKE_INSTALL_BIS} )
    endif()

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmecconfig.h.in
      ADD_CMAKE_INSTALL ${LIBBASE_ADD_CMAKE_INSTALL}
      )

  endif()

  ############################
  # build lib with bctool
  ############################
  if ( FEELMODELS_APP_BC_DESC )

    get_filename_component( BCDESC_FILE ${FEELMODELS_APP_BC_DESC} ABSOLUTE )

    if ( FEELMODELS_APP_BC_MARKERS )
      get_filename_component( BCMARKER_FILE ${FEELMODELS_APP_BC_MARKERS} ABSOLUTE )
    else()
      set(BCMARKER_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.bcmarkers )
    endif()

    if ( FEELMODELS_APP_GEO_DESC )
      get_filename_component( GEO_FILE ${FEELMODELS_APP_GEO_DESC} ABSOLUTE )
    else()
      set(GEO_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.mesh )
    endif()

    if ( NOT ( FEELMODELS_APP_LIB_NAME  ) )
      message(FATAL_ERROR "miss argument! FEELMODELS_APP_LIB_NAME")
    endif()
    set(THE_LIB_NAME ${FEELMODELS_APP_LIB_NAME})

    if ( FEELMODELS_APP_LIB_DIR )
      get_filename_component( THE_LIB_DIR ${FEELMODELS_APP_LIB_DIR} ABSOLUTE )
    else()
      set(THE_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR})
    endif()


    set(FLUIDMECHANICS_INCLUDE_LIBBASE_HEADER 1)

    set(PREFIX_FILES_TO_COPY fluid)
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/codegen_fluidmec.hpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/codegen_fluidmec.cpp )
    set(CODEGEN_SOURCES
      ${THE_LIB_DIR}/codegen_fluidmec.cpp )
    set(LIB_DEPENDS ${LIBBASE_NAME} feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 

    genLibBase(
      LIB_NAME ${THE_LIB_NAME}
      LIB_DIR ${THE_LIB_DIR}
      MARKERS ${BCMARKER_FILE}
      DESC ${BCDESC_FILE}
      GEO ${GEO_FILE}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/feelmodelscoreconfig.h.in
      )
    unset(SOLIDMECHANICS_INCLUDE_LIBBASE_HEADER)

  endif( FEELMODELS_APP_BC_DESC )


 endmacro( genLibFluidMechanics )
#############################################################################
#############################################################################
macro(genExecutableFluidMechanics APPLICATION_NAME
    FEELMODELS_DIM FEELMODELS_ORDERGEO
    FEELMODELS_ORDER_U FEELMODELS_ORDER_P FEELMODELS_P_CONTINUITY
    FEELMODELS_DENSITY_VISCOSITY_ORDER FEELMODELS_DENSITY_VISCOSITY_CONTINUITY FEELMODELS_USE_PERIODICITY
    BCMARKER_FILE BCDESC_FILE GEO_FILE MAIN_FILE )

  resetToZeroAllPhysicalVariables()
  SET(FLUIDMECHANICS 1 )

  set(THE_APPLICATION_DIR ${CMAKE_CURRENT_BINARY_DIR}/codeGen_${APPLICATION_NAME})

  genLibFluidMechanics(
    LIB_NAME   feelpp_model_fluid_${APPLICATION_NAME}
    LIB_DIR    ${THE_APPLICATION_DIR}/fluid
    DIM        ${FEELMODELS_DIM}
    GEO_ORDER  ${FEELMODELS_ORDERGEO}
    U_ORDER    ${FEELMODELS_ORDER_U}
    P_ORDER    ${FEELMODELS_ORDER_P}
    P_CONTINUITY            ${FEELMODELS_P_CONTINUITY}
    DENSITY_VISCOSITY_CONTINUITY  ${FEELMODELS_DENSITY_VISCOSITY_CONTINUITY}
    DENSITY_VISCOSITY_ORDER ${FEELMODELS_DENSITY_VISCOSITY_ORDER}
    USE_PERIODICITY_BIS     ${FEELMODELS_USE_PERIODICITY}
    BC_DESC    ${BCDESC_FILE} 
    BC_MARKERS ${BCMARKER_FILE}
    GEO_DESC   ${GEO_FILE}
    )

  set(LIB_DEPENDS feelpp_model_fluid_${APPLICATION_NAME} feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 
  genExecutableBase(
    ${APPLICATION_NAME}
    APPLICATION_DIR ${THE_APPLICATION_DIR}
    SRC ${MAIN_FILE}
    LIB_DEPENDS ${LIB_DEPENDS}
    CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/feelmodelscoreconfig.h.in
    )

endmacro( genExecutableFluidMechanics )
#############################################################################
#############################################################################

macro(feelpp_add_fluid_application)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;U_ORDER;P_ORDER;P_CONTINUITY;GEO_ORDER;GEO_DESC;BC_MARKERS;BC_DESC;SRC;CFG;DENSITY_VISCOSITY_CONTINUITY;DENSITY_VISCOSITY_ORDER"
    "USE_PERIODICITY"
    ${ARGN}
    )

  CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_GEO_ORDER")
  endif()
  set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
  set(FEELMODELS_ORDER_U ${FEELMODELS_APP_U_ORDER})
  set(FEELMODELS_ORDER_P ${FEELMODELS_APP_P_ORDER})
  set(FEELMODELS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  if ( FEELMODELS_APP_BC_MARKERS )
    set(BCMARKER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_MARKERS})
  else()
    set(BCMARKER_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.bcmarkers )
  endif()

  if ( FEELMODELS_APP_GEO_DESC )
    set(GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GEO_DESC})
  else()
    set(GEO_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/codegen/default.mesh )
  endif()

  if ( NOT FEELMODELS_APP_BC_DESC )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_BC_DESC")
  endif()
  set(DESC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_DESC})
  set(MAIN_FILE ${FEELMODELS_APP_SRC})


  if ( FEELMODELS_APP_P_CONTINUITY )
    set( FEELMODELS_P_CONTINUITY ${FEELMODELS_APP_P_CONTINUITY})
  else()
    set( FEELMODELS_P_CONTINUITY "Continuous")
  endif()

  if (FEELMODELS_APP_DENSITY_VISCOSITY_CONTINUITY)
    set(FEELMODELS_DENSITY_VISCOSITY_CONTINUITY ${FEELMODELS_APP_DENSITY_VISCOSITY_CONTINUITY} )
  else()
    set(FEELMODELS_DENSITY_VISCOSITY_CONTINUITY "Discontinuous" )
  endif()

  if (FEELMODELS_APP_DENSITY_VISCOSITY_ORDER)
    set(FEELMODELS_DENSITY_VISCOSITY_ORDER ${FEELMODELS_APP_DENSITY_VISCOSITY_ORDER} )
  else()
    set(FEELMODELS_DENSITY_VISCOSITY_ORDER 0)
  endif()

  if ( FEELMODELS_APP_USE_PERIODICITY )
    set( FEELMODELS_USE_PERIODICITY 1 )
  else()
    set( FEELMODELS_USE_PERIODICITY 0 )
  endif()

  if ( FEELMODELS_APP_CFG )
    foreach(  cfg ${FEELMODELS_APP_CFG} )
      #      if ( EXISTS ${cfg} )
      # extract cfg filename  to be copied in binary dir
      get_filename_component( CFG_NAME ${cfg} NAME )
      configure_file( ${cfg} ${CFG_NAME} )
      #INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
      #      else()
      #        message(WARNING "Executable ${FEELMODELS_APP_NAME}: configuration file ${cfg} does not exist")
      #      endif()
    endforeach()
  endif(FEELMODELS_APP_CFG)


  if (0)
    MESSAGE("*** Arguments for fluid application ${APPLICATION_NAME}")
    MESSAGE("*** DIM ${FEELMODELS_DIM}")
    MESSAGE("*** ORDER_U ${FEELMODELS_ORDER_U}")
    MESSAGE("*** ORDER_P ${FEELMODELS_ORDER_P}")
    MESSAGE("*** ORDERGEO ${FEELMODELS_ORDERGEO}")
    MESSAGE("*** P_CONTINUITY ${FEELMODELS_P_CONTINUITY}")
    MESSAGE("*** DENSITY_VISCOSITY_CONTINUITY ${FEELMODELS_DENSITY_VISCOSITY_CONTINUITY}")
    MESSAGE("*** MARKERS ${BCMARKER_FILE}")
    MESSAGE("*** DESC ${DESC_FILE}")
    MESSAGE("*** GEO ${GEO_FILE}")
    MESSAGE("*** SRC ${MAIN_FILE}")
  endif()

  genExecutableFluidMechanics(${APPLICATION_NAME}
    ${FEELMODELS_DIM} ${FEELMODELS_ORDERGEO}
    ${FEELMODELS_ORDER_U} ${FEELMODELS_ORDER_P} ${FEELMODELS_P_CONTINUITY}
    ${FEELMODELS_DENSITY_VISCOSITY_ORDER} ${FEELMODELS_DENSITY_VISCOSITY_CONTINUITY}  ${FEELMODELS_USE_PERIODICITY}
    ${BCMARKER_FILE} ${DESC_FILE} ${GEO_FILE} ${MAIN_FILE}
    )

endmacro(feelpp_add_fluid_application)
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
macro(genLibFSI)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;BC_MARKERS;FLUID_U_ORDER;FLUID_P_ORDER;FLUID_P_CONTINUITY;FLUID_GEO_ORDER;FLUID_GEO_DESC;FLUID_BC_DESC;FLUID_DENSITY_VISCOSITY_CONTINUITY;FLUID_DENSITY_VISCOSITY_ORDER;SOLID_DISP_ORDER;SOLID_GEO_ORDER;SOLID_BC_DESC;SOLID_GEO_DESC;SOLID_DENSITY_COEFFLAME_TYPE"
    "FLUID_USE_PERIODICITY;ADD_CMAKE_INSTALL"
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_FLUID_GEO_ORDER OR FEELMODELS_APP_FLUID_U_ORDER OR FEELMODELS_APP_FLUID_P_ORDER OR
        FEELMODELS_APP_SOLID_DISP_ORDER OR FEELMODELS_APP_SOLID_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument!")
  endif()

  ###############################################################
  if ( FEELMODELS_APP_ADD_CMAKE_INSTALL )
    set( LIBBASE_ADD_CMAKE_INSTALL 1 )
  else()
    set( LIBBASE_ADD_CMAKE_INSTALL 0 )
  endif()
  ###############################################################
  # fluid lib
  if ( FEELMODELS_APP_FLUID_USE_PERIODICITY )
    set( FEELMODELS_FLUID_USE_PERIODICITY 1 )
  else()
    set( FEELMODELS_FLUID_USE_PERIODICITY 0 )
  endif()
  genLibFluidMechanics(
    DIM          ${FEELMODELS_APP_DIM}
    GEO_ORDER    ${FEELMODELS_APP_FLUID_GEO_ORDER}
    U_ORDER      ${FEELMODELS_APP_FLUID_U_ORDER}
    P_ORDER      ${FEELMODELS_APP_FLUID_P_ORDER}
    P_CONTINUITY ${FEELMODELS_APP_FLUID_P_CONTINUITY}
    DENSITY_VISCOSITY_CONTINUITY ${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_CONTINUITY}
    DENSITY_VISCOSITY_ORDER      ${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_ORDER}
    USE_PERIODICITY_BIS ${FEELMODELS_FLUID_USE_PERIODICITY}
    ADD_CMAKE_INSTALL_BIS ${LIBBASE_ADD_CMAKE_INSTALL}
    )
  set(FLUID_LIB_NAME ${LIBBASE_NAME})
  set(FLUID_MODEL_SPECIFIC_NAME_SUFFIX ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX}  )
  ###############################################################
  # solid lib
  genLibSolidMechanics(
    DIM ${FEELMODELS_APP_DIM}
    DISP_ORDER ${FEELMODELS_APP_SOLID_DISP_ORDER}
    GEO_ORDER ${FEELMODELS_APP_SOLID_GEO_ORDER}
    DENSITY_COEFFLAME_TYPE ${FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE}
    ADD_CMAKE_INSTALL_BIS ${LIBBASE_ADD_CMAKE_INSTALL}
    )
  set(SOLID_LIB_NAME ${LIBBASE_NAME})
  set(SOLID_MODEL_SPECIFIC_NAME_SUFFIX ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX}  )
  ###############################################################
  # fsi lib base
  set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${FLUID_MODEL_SPECIFIC_NAME_SUFFIX}_${SOLID_MODEL_SPECIFIC_NAME_SUFFIX}  )
  set(FEELMODELS_MODEL_SPECIFIC_NAME fsi_${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})

  set(LIBBASE_DIR ${FEELPP_MODELS_BINARY_DIR}/fsi/${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX} )
  set(LIBBASE_CHECK_PATH ${FEELPP_MODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
  set(LIBBASE_NAME feelpp_model_${FEELMODELS_MODEL_SPECIFIC_NAME})

  set(LIB_DEPENDS ${FLUID_LIB_NAME} ${SOLID_LIB_NAME})

  if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )
    #write empty file in orter to check if this lib has already define
    file(WRITE ${LIBBASE_CHECK_PATH} "")

    # configure libmodelbase
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/fsi/fsi_inst.cpp
      )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/fsi_inst.cpp
      )

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/fsi/fsiconfig.h.in
      ADD_CMAKE_INSTALL ${LIBBASE_ADD_CMAKE_INSTALL}
      )

    # fluid and solid dependencies in install process
    if ( ${LIBBASE_ADD_CMAKE_INSTALL} )
      add_dependencies(install-${LIBBASE_NAME} install-${FLUID_LIB_NAME} install-${SOLID_LIB_NAME} )
    endif()

  endif()

endmacro( genLibFSI )

# macro(genExecutableFluidStructureInteraction APPLICATION_NAME BCMARKER_FILE FEELMODELS_DIM
#     FLUID_ORDER_U FLUID_ORDER_P FLUID_ORDERGEO FLUID_PRESSURE_IS_CONTINUOUS FLUID_USE_CST_DENSITY_VISCOSITY BCDESC_FLUID_FILE MESH_FLUID_FILE
#     SOLID_ORDER_U SOLID_ORDERGEO SOLID_USE_CST_DENSITY_COEFFLAME BCDESC_SOLID_FILE MESH_SOLID_FILE
#     MAIN_FILE)

#   resetToZeroAllPhysicalVariables()
#   SET(FLUIDMECHANICS 1 )
#   SET(SOLIDMECHANICS 1 )

#   set(APPLICATION_DIR codeGen_${APPLICATION_NAME})
#   set(FEELMODELS_DENSITY_VISCOSITY_ORDER 0)
#   set(FEELMODELS_USE_PERIODICITY 0)
#   genLibFluidMechanics(${APPLICATION_NAME} ${APPLICATION_DIR}/fluid ${BCMARKER_FILE} ${FEELMODELS_DIM}
#      ${FLUID_ORDER_U} ${FLUID_ORDER_P} ${FLUID_ORDERGEO} ${FLUID_PRESSURE_IS_CONTINUOUS}
#      ${FLUID_USE_CST_DENSITY_VISCOSITY} ${FEELMODELS_DENSITY_VISCOSITY_ORDER} ${FEELMODELS_USE_PERIODICITY} ${BCDESC_FLUID_FILE} ${MESH_FLUID_FILE} )
#   genLibSolidMechanics(${APPLICATION_NAME} ${APPLICATION_DIR}/solid ${BCMARKER_FILE} ${FEELMODELS_DIM}
#     ${SOLID_ORDER_U} ${SOLID_ORDERGEO} ${SOLID_USE_CST_DENSITY_COEFFLAME} ${BCDESC_SOLID_FILE} ${MESH_SOLID_FILE})

#   set(LIB_DEPENDS fm${APPLICATION_NAME} sm${APPLICATION_NAME} fsimesh fsialg fsicore ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} )
#   genExecutableBase(
#     ${APPLICATION_NAME} ${APPLICATION_DIR}
#     SRC ${MAIN_FILE}
#     LIB_DEPENDS ${LIB_DEPENDS}
#     CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in
#     )

# endmacro( genExecutableFluidStructureInteraction )
# #############################################################################
# #############################################################################

# macro(fsi_add_fsi_application)
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "BC_MARKERS;DIM;FLUID_U_ORDER;FLUID_P_ORDER;FLUID_GEO_ORDER;FLUID_BC_DESC;FLUID_GEO_DESC;FLUID_P_CONTINUITY;FLUID_DENSITY_VISCOSITY_TYPE;SOLID_U_ORDER;SOLID_GEO_ORDER;SOLID_BC_DESC;SOLID_GEO_DESC;SOLID_DENSITY_COEFFLAME_TYPE;SRC"
#     ""
#     ${ARGN}
#     )

#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

#   if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_FLUID_U_ORDER OR FEELMODELS_APP_FLUID_P_ORDER OR FEELMODELS_APP_FLUID_GEO_ORDER ) )
#     message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_FLUID_U_ORDER OR FEELMODELS_APP_FLUID_P_ORDER OR FEELMODELS_APP_FLUID_GEO_ORDER")
#   endif()

#   set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
#   set(FEELMODELS_FLUID_ORDER_U ${FEELMODELS_APP_FLUID_U_ORDER})
#   set(FEELMODELS_FLUID_ORDER_P ${FEELMODELS_APP_FLUID_P_ORDER})
#   set(FEELMODELS_FLUID_ORDERGEO ${FEELMODELS_APP_FLUID_GEO_ORDER})
#   set(FEELMODELS_SOLID_ORDER_U ${FEELMODELS_APP_SOLID_U_ORDER})
#   set(FEELMODELS_SOLID_ORDERGEO ${FEELMODELS_APP_SOLID_GEO_ORDER})

#   if ( FEELMODELS_APP_BC_MARKERS )
#     set(BCMARKER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_MARKERS})
#   else()
#     set(BCMARKER_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.bcmarkers )
#   endif()

#   if ( FEELMODELS_APP_FLUID_GEO_DESC )
#     set(FLUID_GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_FLUID_GEO_DESC})
#   else()
#     set(FLUID_GEO_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#   endif()
#   if ( FEELMODELS_APP_SOLID_GEO_DESC )
#     set(SOLID_GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_SOLID_GEO_DESC})
#   else()
#     set(SOLID_GEO_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#   endif()

#   set(FLUID_DESC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_FLUID_BC_DESC})
#   set(SOLID_DESC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_SOLID_BC_DESC})
#   set(MAIN_FILE ${FEELMODELS_APP_SRC})

#   if ( FEELMODELS_APP_FLUID_P_CONTINUITY )
#     if ("${FEELMODELS_APP_FLUID_P_CONTINUITY}" STREQUAL "Continuous" )
#       set(FLUID_PRESSURE_IS_CONTINUOUS 1)
#     elseif("${FEELMODELS_APP_FLUID_P_CONTINUITY}" STREQUAL "Discontinuous" )
#       set(FLUID_PRESSURE_IS_CONTINUOUS 0)
#     else()
#       message(FATAL_ERROR "FLUID_P_CONTINUITY ${FEELMODELS_APP_FLUID_P_CONTINUITY} : is not valid! It must be Continuous or Discontinuous")
#     endif()
#   else()
#     # default value
#     set(FLUID_PRESSURE_IS_CONTINUOUS 1)
#   endif()

#   if (FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_TYPE)
#     if ("${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_TYPE}" STREQUAL "P0c" )
#       set(FLUID_USE_CST_DENSITY_VISCOSITY 1)
#     elseif ("${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_TYPE}" STREQUAL "P0d" )
#       set(FLUID_USE_CST_DENSITY_VISCOSITY 0)
#     else()
#       message(FATAL_ERROR "FLUID_DENSITY_VISCOSITY_TYPE ${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_TYPE} : is not valid! It must be P0c or P0d")
#     endif()
#   else()
#     # default value
#     set(FLUID_USE_CST_DENSITY_VISCOSITY 1)
#   endif()

#   if (FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE)
#     if ("${FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE}" STREQUAL "P0c" )
#       set(SOLID_USE_CST_DENSITY_COEFFLAME 1)
#     elseif ("${FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE}" STREQUAL "P0d" )
#       set(SOLID_USE_CST_DENSITY_COEFFLAME 0)
#     else()
#       message(FATAL_ERROR "SOLID_DENSITY_COEFFLAME_TYPE ${FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE} : is not valid! It must be P0c or P0d")
#     endif()
#   else()
#     # default value
#     set(SOLID_USE_CST_DENSITY_COEFFLAME 1)
#   endif()

#   if (0)
#     MESSAGE("*** Arguments for fsi application ${APPLICATION_NAME}")
#     MESSAGE("*** MARKERS ${BCMARKER_FILE}")
#     MESSAGE("*** DIM ${FEELMODELS_DIM}")
#     MESSAGE("*** FLUID_ORDER_U ${FEELMODELS_FLUID_ORDER_U}")
#     MESSAGE("*** FLUID_ORDER_P ${FEELMODELS_FLUID_ORDER_P}")
#     MESSAGE("*** FLUID_ORDERGEO ${FEELMODELS_FLUID_ORDERGEO}")
#     MESSAGE("*** FLUID_PRESSURE_IS_CONTINUOUS ${FLUID_PRESSURE_IS_CONTINUOUS}")
#     MESSAGE("*** FLUID_USE_CST_DENSITY_VISCOSITY ${FLUID_USE_CST_DENSITY_VISCOSITY}")
#     MESSAGE("*** FLUID_DESC ${FLUID_DESC_FILE}")
#     MESSAGE("*** FLUID_GEO ${FLUID_GEO_FILE}")
#     MESSAGE("*** SOLID_DESC ${SOLID_DESC_FILE}")
#     MESSAGE("*** SOLID_GEO ${SOLID_GEO_FILE}")
#     MESSAGE("*** SRC ${MAIN_FILE}")
#   endif()

#   genExecutableFluidStructureInteraction(${APPLICATION_NAME} ${BCMARKER_FILE} ${FEELMODELS_DIM}
#     ${FEELMODELS_FLUID_ORDER_U} ${FEELMODELS_FLUID_ORDER_P} ${FEELMODELS_FLUID_ORDERGEO} ${FLUID_PRESSURE_IS_CONTINUOUS} ${FLUID_USE_CST_DENSITY_VISCOSITY}  ${FLUID_DESC_FILE} ${FLUID_GEO_FILE}
#     ${FEELMODELS_SOLID_ORDER_U} ${FEELMODELS_SOLID_ORDERGEO} ${SOLID_USE_CST_DENSITY_COEFFLAME} ${SOLID_DESC_FILE} ${SOLID_GEO_FILE}
#     ${MAIN_FILE} )

# endmacro(fsi_add_fsi_application)


