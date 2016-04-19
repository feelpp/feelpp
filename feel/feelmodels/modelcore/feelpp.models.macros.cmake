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
    "DIM;T_ORDER;GEO_ORDER;"
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


endmacro(genLibThermoDynamics)

#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro( genLibSolidMechanics )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;DISP_ORDER;GEO_ORDER;DENSITY_COEFFLAME_TYPE;ADD_CMAKE_INSTALL_BIS"
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

endmacro( genLibSolidMechanics )
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
macro(genLibFluidMechanics)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;U_ORDER;P_ORDER;P_CONTINUITY;GEO_ORDER;DENSITY_VISCOSITY_CONTINUITY;DENSITY_VISCOSITY_ORDER;USE_PERIODICITY_BIS;ADD_CMAKE_INSTALL_BIS"
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

 endmacro( genLibFluidMechanics )
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

