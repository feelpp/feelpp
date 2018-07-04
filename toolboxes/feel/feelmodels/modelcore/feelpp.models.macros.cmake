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
  SET(THERMOELECTRIC 0 )
endmacro(resetToZeroAllPhysicalVariables)
#############################################################################
# create generic simple lib : ${PREFIX_LIB}${APPLICATION_NAME}
#############################################################################
macro(genLibBase)
  PARSE_ARGUMENTS(FEELMODELS_GENLIB_BASE
    "LIB_NAME;LIB_DIR;LIB_DEPENDS;PREFIX_INCLUDE_USERCONFIG;FILES_TO_COPY;FILES_SOURCES;CONFIG_PATH"
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
  set_property(TARGET ${LIB_APPLICATION_NAME} PROPERTY MACOSX_RPATH ON)

  if( FEELPP_ENABLE_PCH_MODELS )
      add_precompiled_header( ${LIB_APPLICATION_NAME} )
  endif()

  # install process
  INSTALL(TARGETS ${LIB_APPLICATION_NAME} DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)

endmacro(genLibBase)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibThermoDynamics )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;T_ORDER;GEO_ORDER;"
    "NO_UPDATE_MODEL_DEF"
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

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodynconfig.h.in
      )

  endif()


endmacro(genLibThermoDynamics)

#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro( genLibSolidMechanics )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;DISP_ORDER;GEO_ORDER;DENSITY_COEFFLAME_TYPE"
    "NO_UPDATE_MODEL_DEF"
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

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/solid/solidmecconfig.h.in
      )

  endif()

endmacro( genLibSolidMechanics )

#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro(genLibFluidMechanics)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;U_ORDER;P_ORDER;P_CONTINUITY;GEO_ORDER;DENSITY_VISCOSITY_CONTINUITY;DENSITY_VISCOSITY_ORDER;USE_PERIODICITY_BIS"
    "USE_PERIODICITY;NO_UPDATE_MODEL_DEF"
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
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicscreate_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsothers_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsupdatelinear_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsupdatelinearbc_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsupdatejacobian_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsupdatejacobianbc_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsupdateresidual_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsupdateresidualbc_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsupdatestabilisation_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsupdatestabilisationgls_inst.cpp
      )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/fluidmechanicscreate_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsothers_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsupdatelinear_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsupdatelinearbc_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsupdatejacobian_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsupdatejacobianbc_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsupdateresidual_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsupdateresidualbc_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsupdatestabilisation_inst.cpp
      ${LIBBASE_DIR}/fluidmechanicsupdatestabilisationgls_inst.cpp
      )
    set(LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} )
    if (FEELPP_MODELS_ENABLE_MESHALE )
      set(LIB_DEPENDS feelpp_modelmeshale ${LIB_DEPENDS})
    endif()
    # thermodynamcis depend
    set(LIB_DEPENDS feelpp_model_thermodyn${FLUIDMECHANICS_DIM}dP${FLUIDMECHANICS_ORDERGEO}G${FLUIDMECHANICS_ORDERGEO} ${LIB_DEPENDS})

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/fluid/fluidmechanicsconfig.h.in
      )

  endif()

endmacro( genLibFluidMechanics )

#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro(genLibFSI)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;BC_MARKERS;FLUID_U_ORDER;FLUID_P_ORDER;FLUID_P_CONTINUITY;FLUID_GEO_ORDER;FLUID_GEO_DESC;FLUID_BC_DESC;FLUID_DENSITY_VISCOSITY_CONTINUITY;FLUID_DENSITY_VISCOSITY_ORDER;SOLID_DISP_ORDER;SOLID_GEO_ORDER;SOLID_BC_DESC;SOLID_GEO_DESC;SOLID_DENSITY_COEFFLAME_TYPE"
    "FLUID_USE_PERIODICITY"
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_FLUID_GEO_ORDER OR FEELMODELS_APP_FLUID_U_ORDER OR FEELMODELS_APP_FLUID_P_ORDER OR
        FEELMODELS_APP_SOLID_DISP_ORDER OR FEELMODELS_APP_SOLID_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument!")
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
      )


  endif()

endmacro( genLibFSI )

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibAdvection )
  PARSE_ARGUMENTS(FEELMODELS_APP
      "DIM;POLY_ORDER;CONTINUITY;POLY_SET;GEO_ORDER;DIFFUSION_REACTION_ORDER;DIFFUSIONCOEFF_POLY_SET;REACTIONCOEFF_POLY_SET;DIFFUSION_REACTION_CONTINUITY"
      "NO_UPDATE_MODEL_DEF"
      ${ARGN}
    )

  if ( NOT FEELMODELS_APP_NO_UPDATE_MODEL_DEF )
    resetToZeroAllPhysicalVariables()
    SET(ADVECTION 1 )
  endif()

  #CAR(LIB_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(ADVECTION_DIM ${FEELMODELS_APP_DIM})
  set(ADVECTION_ORDERPOLY ${FEELMODELS_APP_POLY_ORDER})

  # continuity
  if (FEELMODELS_APP_CONTINUITY)
      if ("${FEELMODELS_APP_CONTINUITY}" STREQUAL "Continuous" )
          set(ADVECTION_USE_CONTINUOUS 1)
          unset(ADVECTION_CONTINUITY_TAG)
      elseif ("${FEELMODELS_APP_CONTINUITY}" STREQUAL "Discontinuous" )
          set(ADVECTION_USE_CONTINUOUS 0)
          set(ADVECTION_CONTINUITY_TAG d)
      else()
          message(FATAL_ERROR "CONTINUITY ${FEELMODELS_APP_CONTINUITY} : is not valid! It must be Continuous or Discontinuous")
      endif()
  else()
      # default value (Continuous)
    set(ADVECTION_USE_CONTINUOUS 1)
    unset(ADVECTION_CONTINUITY_TAG) 
  endif()

  # geometric order
  set(ADVECTION_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  # polynomial set
  if(FEELMODELS_APP_POLY_SET)
      set(ADVECTION_POLYSET ${FEELMODELS_APP_POLY_SET})
      if("${FEELMODELS_APP_POLY_SET}" STREQUAL "Scalar" )
          unset(ADVECTION_POLY_SET_TAG)
      elseif("${FEELMODELS_APP_POLY_SET}" STREQUAL "Vectorial" )
          set(ADVECTION_POLY_SET_TAG Vec)
      else()
          message(FATAL_ERROR "POLY_SET ${FEELMODELS_APP_POLY_SET} is not valid! It must be Scalar or Vectorial")
      endif()
  else()
      #default value
      set(ADVECTION_POLYSET Scalar)
      unset(ADVECTION_POLY_SET_TAG)
  endif()
  #######################################################
  if (FEELMODELS_APP_DIFFUSION_REACTION_ORDER)
      set(ADVECTION_ORDER_DIFFUSION_REACTION ${FEELMODELS_APP_DIFFUSION_REACTION_ORDER} )
  else()
    # default value
    set(ADVECTION_ORDER_DIFFUSION_REACTION ${FEELMODELS_APP_POLY_ORDER} )
  endif()
  if(FEELMODELS_APP_DIFFUSIONCOEFF_POLY_SET)
      set(ADVECTION_DIFFUSIONCOEFF_POLYSET ${FEELMODELS_APP_DIFFUSIONCOEFF_POLY_SET})
      if("${FEELMODELS_APP_DIFFUSIONCOEFF_POLY_SET}" STREQUAL "Scalar" )
          unset(ADVECTION_DIFFUSIONCOEFF_POLY_SET_TAG)
      elseif("${FEELMODELS_APP_DIFFUSIONCOEFF_POLY_SET}" STREQUAL "Tensor2" )
          set(ADVECTION_DIFFUSIONCOEFF_POLY_SET_TAG t)
      elseif("${FEELMODELS_APP_DIFFUSIONCOEFF_POLY_SET}" STREQUAL "Tensor2Symm" )
          set(ADVECTION_DIFFUSIONCOEFF_POLY_SET_TAG ts)
      else()
          message(FATAL_ERROR "DIFFUSIONCOEFF_POLY_SET ${FEELMODELS_APP_DIFFUSIONCOEFF_POLY_SET} is not valid! It must be Scalar or Tensor2 or Tensor2Symm")
      endif()
  else()
      #default value
      set(ADVECTION_DIFFUSIONCOEFF_POLYSET Scalar)
      unset(ADVECTION_DIFFUSIONCOEFF_POLY_SET_TAG)
  endif()
  if(FEELMODELS_APP_REACTIONCOEFF_POLY_SET)
      set(ADVECTION_REACTIONCOEFF_POLYSET ${FEELMODELS_APP_REACTIONCOEFF_POLY_SET})
      if("${FEELMODELS_APP_REACTIONCOEFF_POLY_SET}" STREQUAL "Scalar" )
          unset(ADVECTION_REACTIONCOEFF_POLY_SET_TAG)
      elseif("${FEELMODELS_APP_REACTIONCOEFF_POLY_SET}" STREQUAL "Tensor2" )
          set(ADVECTION_REACTIONCOEFF_POLY_SET_TAG t)
      elseif("${FEELMODELS_APP_REACTIONCOEFF_POLY_SET}" STREQUAL "Tensor2Symm" )
          set(ADVECTION_REACTIONCOEFF_POLY_SET_TAG ts)
      else()
          message(FATAL_ERROR "REACTIONCOEFF_POLY_SET ${FEELMODELS_APP_REACTIONCOEFF_POLY_SET} is not valid! It must be Scalar or Tensor2 or Tensor2Symm")
      endif()
  else()
      #default value
      set(ADVECTION_REACTIONCOEFF_POLYSET Scalar)
      unset(ADVECTION_REACTIONCOEFF_POLY_SET_TAG)
  endif()
  #######################################################
  if (FEELMODELS_APP_DIFFUSION_REACTION_CONTINUITY)
      if ("${FEELMODELS_APP_DIFFUSION_REACTION_CONTINUITY}" STREQUAL "Continuous" )
          set(ADVECTION_USE_CONTINUOUS_DIFFUSION_REACTION 1)
          unset(ADVECTION_DIFFUSION_REACTION_CONTINUITY_TAG)
      elseif ("${FEELMODELS_APP_DIFFUSION_REACTION_CONTINUITY}" STREQUAL "Discontinuous" )
          set(ADVECTION_USE_CONTINUOUS_DIFFUSION_REACTION 0)
          set(ADVECTION_DIFFUSION_REACTION_CONTINUITY_TAG d)
      else()
          message(FATAL_ERROR "DIFFUSION_REACTION_CONTINUITY ${FEELMODELS_APP_DIFFUSION_REACTION_CONTINUITY} : is not valid! It must be Continuous or Discontinuous")
      endif()
  else()
    # default value
    set(ADVECTION_USE_CONTINUOUS_DIFFUSION_REACTION 1)
    unset(ADVECTION_DIFFUSION_REACTION_CONTINUITY_TAG)
  endif()
  set(ADVECTION_DIFFUSION_REACTION_TAG DRP${ADVECTION_ORDER_DIFFUSION_REACTION}${ADVECTION_DIFFUSIONCOEFF_POLY_SET_TAG}${ADVECTION_DIFFUSION_REACTION_CONTINUITY_TAG}) # default value P${ORDER_POLY}${POLY_SET}
  #######################################################

  if (0)
    MESSAGE("*** Arguments for advection application ${LIB_NAME}")
    MESSAGE("*** DIM ${ADVECTION_DIM}")
    MESSAGE("*** ORDERPOLY ${ADVECTION_ORDERPOLY}")
    MESSAGE("*** ORDERGEO ${ADVECTION_ORDERGEO}")
  endif()

  set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${ADVECTION_DIM}dP${ADVECTION_ORDERPOLY}${ADVECTION_CONTINUITY_TAG}${ADVECTION_POLY_SET_TAG}G${ADVECTION_ORDERGEO}${ADVECTION_DIFFUSION_REACTION_TAG} )
  set(FEELMODELS_MODEL_SPECIFIC_NAME advection${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})
  #set(FEELMODELS_MODEL_SPECIFIC_NAME ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})
  set(LIBBASE_DIR ${FEELPP_MODELS_BINARY_DIR}/advection/${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX} )
  set(LIBBASE_CHECK_PATH ${FEELPP_MODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
  set(LIBBASE_NAME feelpp_model_${FEELMODELS_MODEL_SPECIFIC_NAME})

  if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )

    #write empty file in orter to check if this lib has already define
    file(WRITE ${LIBBASE_CHECK_PATH} "")

    # configure libmodelbase
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/advection/advectionbase_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/advection/advection_inst.cpp )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/advectionbase_inst.cpp
      ${LIBBASE_DIR}/advection_inst.cpp )
    set(LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/advection/advectionconfig.h.in
      )

  endif()


endmacro(genLibAdvection)
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibLevelset )
    PARSE_ARGUMENTS(FEELMODELS_APP
        "DIM;LEVELSET_ORDER;LEVELSET_PN_ORDER;GEO_ORDER;ADVECTION_DIFFUSION_REACTION_ORDER;ADVECTION_DIFFUSION_REACTION_CONTINUITY"
        "NO_UPDATE_MODEL_DEF"
        ${ARGN}
        )

    if ( NOT FEELMODELS_APP_NO_UPDATE_MODEL_DEF )
        resetToZeroAllPhysicalVariables()
        SET(LEVELSET 1 )
    endif()

    #CAR(LIB_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

    if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_LEVELSET_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
        message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_LEVELSET_ORDER OR  FEELMODELS_APP_GEO_ORDER")
    endif()

    set(LEVELSET_DIM ${FEELMODELS_APP_DIM})
    set(LEVELSET_ORDERPOLY ${FEELMODELS_APP_LEVELSET_ORDER})
    set(LEVELSET_PN_ORDERPOLY ${FEELMODELS_APP_LEVELSET_PN_ORDER})
    set(LEVELSET_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
    #######################################################
    if(FEELMODELS_APP_ADVECTION_DIFFUSION_REACTION_ORDER)
        set(LEVELSET_ADVECTION_DIFFUSION_REACTION_ORDER ${FEELMODELS_APP_ADVECTION_DIFFUSION_REACTION_ORDER})
    else()
        set(LEVELSET_ADVECTION_DIFFUSION_REACTION_ORDER ${LEVELSET_ORDERPOLY})
    endif()
    if(FEELMODELS_APP_ADVECTION_DIFFUSION_REACTION_CONTINUITY)
        set(LEVELSET_ADVECTION_DIFFUSION_REACTION_CONTINUITY ${FEELMODELS_APP_ADVECTION_DIFFUSION_REACTION_CONTINUITY})
    else()
        set(LEVELSET_ADVECTION_DIFFUSION_REACTION_CONTINUITY Continuous)
    endif()
    #######################################################
    # advection
    genLibAdvection(
        DIM          ${FEELMODELS_APP_DIM}
        GEO_ORDER    ${FEELMODELS_APP_GEO_ORDER}
        POLY_ORDER      ${FEELMODELS_APP_LEVELSET_ORDER}
        DIFFUSION_REACTION_ORDER      ${LEVELSET_ADVECTION_DIFFUSION_REACTION_ORDER}
        DIFFUSION_REACTION_CONTINUITY ${LEVELSET_ADVECTION_DIFFUSION_REACTION_CONTINUITY}
        )
    set(ADVECTION_LIB_NAME ${LIBBASE_NAME})
    set(ADVECTION_MODEL_SPECIFIC_NAME_SUFFIX ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX}  )
    # vectorial advection (backward characteristics)
    genLibAdvection(
        DIM          ${FEELMODELS_APP_DIM}
        GEO_ORDER    ${FEELMODELS_APP_GEO_ORDER}
        POLY_ORDER      ${FEELMODELS_APP_LEVELSET_ORDER}
        POLY_SET        Vectorial
        DIFFUSION_REACTION_ORDER      ${LEVELSET_ADVECTION_DIFFUSION_REACTION_ORDER}
        DIFFUSION_REACTION_CONTINUITY ${LEVELSET_ADVECTION_DIFFUSION_REACTION_CONTINUITY}
        )
    set(ADVECTIONVEC_LIB_NAME ${LIBBASE_NAME})
    set(ADVECTIONVEC_MODEL_SPECIFIC_NAME_SUFFIX ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX}  )
    #######################################################

    set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${LEVELSET_DIM}dP${LEVELSET_ORDERPOLY}G${LEVELSET_ORDERGEO} )
    set(FEELMODELS_MODEL_SPECIFIC_NAME levelset${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})

    set(LIBBASE_DIR ${FEELPP_MODELS_BINARY_DIR}/levelset/${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX} )
    set(LIBBASE_CHECK_PATH ${FEELPP_MODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
    set(LIBBASE_NAME feelpp_model_${FEELMODELS_MODEL_SPECIFIC_NAME})

    set(LIB_DEPENDS ${ADVECTION_LIB_NAME} ${ADVECTIONVEC_LIB_NAME} ${FEELPP_LIBRARIES} )

    if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )

        #write empty file in orter to check if this lib has already define
        file(WRITE ${LIBBASE_CHECK_PATH} "")

        # configure libmodelbase
        set(CODEGEN_FILES_TO_COPY
            ${FEELPP_MODELS_SOURCE_DIR}/levelset/levelset_inst.cpp
            ${FEELPP_MODELS_SOURCE_DIR}/levelset/levelsetadvection_inst.cpp
            ${FEELPP_MODELS_SOURCE_DIR}/levelset/parameter_map.cpp )
        set(CODEGEN_SOURCES
            ${LIBBASE_DIR}/levelset_inst.cpp
            ${LIBBASE_DIR}/levelsetadvection_inst.cpp
            ${LIBBASE_DIR}/parameter_map.cpp )
        #set(LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 

        # generate libmodelbase
        genLibBase(
            LIB_NAME ${LIBBASE_NAME}
            LIB_DIR ${LIBBASE_DIR}
            LIB_DEPENDS ${LIB_DEPENDS}
            PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
            FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
            FILES_SOURCES ${CODEGEN_SOURCES}
            CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/levelset/levelsetconfig.h.in
            )

    endif()
endmacro(genLibLevelset)
#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro(genLibMultiFluid)
  PARSE_ARGUMENTS(FEELMODELS_APP
      "DIM;FLUID_U_ORDER;FLUID_P_ORDER;FLUID_P_CONTINUITY;GEO_ORDER;FLUID_DENSITY_VISCOSITY_CONTINUITY;FLUID_DENSITY_VISCOSITY_ORDER;LEVELSET_ORDER;LEVELSET_PN_ORDER;LEVELSET_DIFFUSION_REACTION_ORDER;LEVELSET_DIFFUSION_REACTION_CONTINUITY"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_FLUID_GEO_ORDER OR FEELMODELS_APP_FLUID_U_ORDER OR FEELMODELS_APP_FLUID_P_ORDER OR
      LEVELSET_ORDER ) )
    message(FATAL_ERROR "miss argument!")
  endif()

  ###############################################################
  # fluid lib
  genLibFluidMechanics(
    DIM          ${FEELMODELS_APP_DIM}
    GEO_ORDER    ${FEELMODELS_APP_GEO_ORDER}
    U_ORDER      ${FEELMODELS_APP_FLUID_U_ORDER}
    P_ORDER      ${FEELMODELS_APP_FLUID_P_ORDER}
    P_CONTINUITY ${FEELMODELS_APP_FLUID_P_CONTINUITY}
    DENSITY_VISCOSITY_CONTINUITY ${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_CONTINUITY}
    DENSITY_VISCOSITY_ORDER      ${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_ORDER}
    )
  set(FLUID_LIB_NAME ${LIBBASE_NAME})
  set(FLUID_MODEL_SPECIFIC_NAME_SUFFIX ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX}  )
  ###############################################################
  # levelset lib
  genLibLevelset(
    DIM ${FEELMODELS_APP_DIM}
    LEVELSET_ORDER ${FEELMODELS_APP_LEVELSET_ORDER}
    LEVELSET_PN_ORDER ${FEELMODELS_APP_LEVELSET_PN_ORDER}
    GEO_ORDER ${FEELMODELS_APP_GEO_ORDER}
    ADVECTION_DIFFUSION_REACTION_ORDER ${FEELMODELS_APP_LEVELSET_DIFFUSION_REACTION_ORDER}
    ADVECTION_DIFFUSION_REACTION_CONTINUITY ${FEELMODELS_APP_LEVELSET_DIFFUSION_REACTION_CONTINUITY}
    )
  set(LEVELSET_LIB_NAME ${LIBBASE_NAME})
  set(LEVELSET_MODEL_SPECIFIC_NAME_SUFFIX ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX}  )
  ###############################################################
  # multifluid lib base
  set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${FLUID_MODEL_SPECIFIC_NAME_SUFFIX}_${LEVELSET_MODEL_SPECIFIC_NAME_SUFFIX}  )
  set(FEELMODELS_MODEL_SPECIFIC_NAME multifluid_${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})

  set(LIBBASE_DIR ${FEELPP_MODELS_BINARY_DIR}/multifluid/${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX} )
  set(LIBBASE_CHECK_PATH ${FEELPP_MODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
  set(LIBBASE_NAME feelpp_model_${FEELMODELS_MODEL_SPECIFIC_NAME})

  set(LIB_DEPENDS ${FLUID_LIB_NAME} ${LEVELSET_LIB_NAME})

  if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )
    #write empty file in orter to check if this lib has already define
    file(WRITE ${LIBBASE_CHECK_PATH} "")

    # configure libmodelbase
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/multifluid/multifluid_inst.cpp
      )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/multifluid_inst.cpp
      )

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/multifluid/multifluidconfig.h.in
      )


  endif()

endmacro( genLibMultiFluid )

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibThermoElectric )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;T_ORDER;GEO_ORDER;"
    "NO_UPDATE_MODEL_DEF"
    ${ARGN}
    )

  if ( NOT FEELMODELS_APP_NO_UPDATE_MODEL_DEF )
    resetToZeroAllPhysicalVariables()
    SET(THERMOELECTRIC 1 )
  endif()

  #CAR(LIB_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(THERMOELECTRIC_DIM ${FEELMODELS_APP_DIM})
  set(THERMOELECTRIC_ORDERPOLY ${FEELMODELS_APP_T_ORDER})
  set(THERMOELECTRIC_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  if (0)
    MESSAGE("*** Arguments for fluid application ${LIB_NAME}")
    MESSAGE("*** DIM ${THERMOELECTRIC_DIM}")
    MESSAGE("*** ORDERPOLY ${THERMOELECTRIC_ORDERPOLY}")
    MESSAGE("*** ORDERGEO ${THERMOELECTRIC_ORDERGEO}")
  endif()

  set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${THERMOELECTRIC_DIM}dP${THERMOELECTRIC_ORDERPOLY}G${THERMOELECTRIC_ORDERGEO} )
  set(FEELMODELS_MODEL_SPECIFIC_NAME thermoelectric${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})
  #set(FEELMODELS_MODEL_SPECIFIC_NAME ${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})
  set(LIBBASE_DIR ${FEELPP_MODELS_BINARY_DIR}/thermoelectric/${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX} )
  set(LIBBASE_CHECK_PATH ${FEELPP_MODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
  set(LIBBASE_NAME feelpp_model_${FEELMODELS_MODEL_SPECIFIC_NAME})

  if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )

    #write empty file in orter to check if this lib has already define
    file(WRITE ${LIBBASE_CHECK_PATH} "")

    # configure libmodelbase
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/thermoelectric/thermoelectric_inst.cpp )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/thermoelectric_inst.cpp )
    set(LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 
    # thermodynamcis depend
    set(LIB_DEPENDS feelpp_model_thermodyn${THERMOELECTRIC_DIM}dP${THERMOELECTRIC_ORDERPOLY}G${THERMOELECTRIC_ORDERGEO} ${LIB_DEPENDS})

    # generate libmodelbase
    genLibBase(
      LIB_NAME ${LIBBASE_NAME}
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/thermoelectric/thermoelectricconfig.h.in
      )

  endif()


endmacro(genLibThermoElectric)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
