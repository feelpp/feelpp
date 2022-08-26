macro(feelpp_toolboxes_add_library)
  PARSE_ARGUMENTS(FEELPP_TOOLBOXES_ADD_LIB
    "SRCS;LINK_LIBRARIES;DEPS"
    ""
    ${ARGN}
    )

  CAR(FEELPP_TOOLBOXES_LIB_BASE_NAME ${FEELPP_TOOLBOXES_ADD_LIB_DEFAULT_ARGS})
  set(FEELPP_TOOLBOXES_TARGET_NAME feelpp_toolbox_${FEELPP_TOOLBOXES_LIB_BASE_NAME})
  add_library(
    ${FEELPP_TOOLBOXES_TARGET_NAME} SHARED ${FEELPP_TOOLBOXES_ADD_LIB_SRCS}
   )
  target_link_libraries(${FEELPP_TOOLBOXES_TARGET_NAME} PUBLIC ${FEELPP_TOOLBOXES_ADD_LIB_LINK_LIBRARIES} )
  add_library(Feelpp::${FEELPP_TOOLBOXES_TARGET_NAME} ALIAS ${FEELPP_TOOLBOXES_TARGET_NAME} )  # to match exported target
  set_target_properties(${FEELPP_TOOLBOXES_TARGET_NAME} PROPERTIES OUTPUT_NAME "${FEELPP_TOOLBOXES_TARGET_NAME}")
  target_include_directories(${FEELPP_TOOLBOXES_TARGET_NAME} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/feelpp/toolboxes>  )
  set_property(TARGET ${FEELPP_TOOLBOXES_TARGET_NAME} PROPERTY MACOSX_RPATH ON)
  if ( FEELPP_TOOLBOXES_ADD_LIB_DEPS )
    add_dependencies(${FEELPP_TOOLBOXES_TARGET_NAME} ${FEELPP_TOOLBOXES_ADD_LIB_DEPS})
  endif()
  if( FEELPP_ENABLE_PCH_MODELS )
    add_precompiled_header( ${FEELPP_TOOLBOXES_TARGET_NAME} )
  endif()
  #include(GNUInstallDirs)
  set_target_properties(${FEELPP_TOOLBOXES_TARGET_NAME} PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
  install(TARGETS ${FEELPP_TOOLBOXES_TARGET_NAME} EXPORT feelpp-toolboxes-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

endmacro(feelpp_toolboxes_add_library)

#############################################################################
# create generic simple lib
#############################################################################
macro(genLibBase)
  PARSE_ARGUMENTS(FEELMODELS_GENLIB_BASE
    "LIB_NAME;LIB_DIR;LIB_DEPENDS;FILES_TO_COPY;FILES_SOURCES;CONFIG_PATH;TARGET_COPY_FILES"
    ""
    ${ARGN}
    )

  set( FEELMODELS_GENLIB_APPLICATION_DIR ${FEELMODELS_GENLIB_BASE_LIB_DIR} )
  #CAR(APPLICATION_NAME      ${FEELMODELS_GENLIB_BASE_DEFAULT_ARGS})
  #CDR(FEELMODELS_GENLIB_APPLICATION_DIR       ${FEELMODELS_GENLIB_BASE_DEFAULT_ARGS})

  set(LIB_DEPENDS           ${FEELMODELS_GENLIB_BASE_LIB_DEPENDS})
  set(LIB_APPLICATION_NAME ${FEELMODELS_GENLIB_BASE_LIB_NAME})

  set(CODEGEN_FILES_TO_COPY ${FEELMODELS_GENLIB_BASE_FILES_TO_COPY})
  set(CODEGEN_SOURCES       ${FEELMODELS_GENLIB_BASE_FILES_SOURCES})

  if ( FEELMODELS_GENLIB_BASE_CONFIG_PATH )
    set(FEELMODELS_GENLIB_CONFIG_PATH       ${FEELMODELS_GENLIB_BASE_CONFIG_PATH})
    get_filename_component(FEELMODELS_GENLIB_CONFIG_FILENAME_WE ${FEELMODELS_GENLIB_CONFIG_PATH} NAME_WE)
    CONFIGURE_FILE( ${FEELMODELS_GENLIB_CONFIG_PATH} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${FEELMODELS_GENLIB_CONFIG_FILENAME_WE}.h  )
  endif()

  if ( FEELMODELS_GENLIB_BASE_TARGET_COPY_FILES )
    set( TARGET_COPY_FILES ${FEELMODELS_GENLIB_BASE_TARGET_COPY_FILES} )
  else()
    set( TARGET_COPY_FILES ${LIB_APPLICATION_NAME}_copyfiles)
    add_custom_target(${TARGET_COPY_FILES}  ALL COMMENT "Copying modified files"  )
  endif()

  # lib files
  foreach(filepath ${CODEGEN_FILES_TO_COPY})
    get_filename_component(filename ${filepath} NAME)
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} )
      #configure_file( ${filepath} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} COPYONLY)
      file(WRITE ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} "") #write empty file
    endif()
    add_custom_command(TARGET ${TARGET_COPY_FILES} COMMAND ${CMAKE_COMMAND} -E copy_if_different
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
  add_dependencies(${LIB_APPLICATION_NAME} ${TARGET_COPY_FILES})
  target_link_libraries(${LIB_APPLICATION_NAME} ${LIB_DEPENDS} )
  set_target_properties(${LIB_APPLICATION_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${FEELMODELS_GENLIB_APPLICATION_DIR}")
  set_property(TARGET ${LIB_APPLICATION_NAME} PROPERTY MACOSX_RPATH ON)

  if( FEELPP_ENABLE_PCH_MODELS )
      add_precompiled_header( ${LIB_APPLICATION_NAME} )
  endif()

  # install process
  set_target_properties(${LIB_APPLICATION_NAME} PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
  INSTALL(TARGETS ${LIB_APPLICATION_NAME} EXPORT feelpp-toolboxes-targets
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )


endmacro(genLibBase)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibHeat )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;T_ORDER;GEO_ORDER"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(HEAT_DIM ${FEELMODELS_APP_DIM})
  set(HEAT_ORDERPOLY ${FEELMODELS_APP_T_ORDER})
  set(HEAT_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  set(HEAT_LIB_VARIANTS ${HEAT_DIM}dP${HEAT_ORDERPOLY}G${HEAT_ORDERGEO} )
  set(HEAT_LIB_NAME feelpp_toolbox_heat_lib_${HEAT_LIB_VARIANTS})

  if ( NOT TARGET ${HEAT_LIB_NAME} )
    # configure the lib
    set(HEAT_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/heat/${HEAT_LIB_VARIANTS})
    set(HEAT_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heat/heat_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heat/heatassemblylinear_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heat/heatassemblyjacobian_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heat/heatassemblyresidual_inst.cpp
      )
    set(HEAT_CODEGEN_SOURCES
      ${HEAT_LIB_DIR}/heat_inst.cpp
      ${HEAT_LIB_DIR}/heatassemblylinear_inst.cpp
      ${HEAT_LIB_DIR}/heatassemblyjacobian_inst.cpp
      ${HEAT_LIB_DIR}/heatassemblyresidual_inst.cpp
      )
    set(HEAT_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials feelpp_toolbox_heatbase ) 
    # generate the lib target
    genLibBase(
      LIB_NAME ${HEAT_LIB_NAME}
      LIB_DIR ${HEAT_LIB_DIR}
      LIB_DEPENDS ${HEAT_LIB_DEPENDS}
      FILES_TO_COPY ${HEAT_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${HEAT_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heat/heatconfig.h.in
      )
  endif()
endmacro(genLibHeat)


#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibElectric )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;P_ORDER;GEO_ORDER;"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_P_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_P_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(ELECTRIC_DIM ${FEELMODELS_APP_DIM})
  set(ELECTRIC_ORDERPOLY ${FEELMODELS_APP_P_ORDER})
  set(ELECTRIC_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  set(ELECTRIC_LIB_VARIANTS ${ELECTRIC_DIM}dP${ELECTRIC_ORDERPOLY}G${ELECTRIC_ORDERGEO} )
  set(ELECTRIC_LIB_NAME feelpp_toolbox_electric_lib_${ELECTRIC_LIB_VARIANTS})

  if ( NOT TARGET ${ELECTRIC_LIB_NAME} )
    # configure the lib
    set(ELECTRIC_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/electric/${ELECTRIC_LIB_VARIANTS})
    set(ELECTRIC_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/electric/electric_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/electric/electricassemblylinear_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/electric/electricassemblyjacobian_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/electric/electricassemblyresidual_inst.cpp
      )
    set(ELECTRIC_CODEGEN_SOURCES
      ${ELECTRIC_LIB_DIR}/electric_inst.cpp
      ${ELECTRIC_LIB_DIR}/electricassemblylinear_inst.cpp
      ${ELECTRIC_LIB_DIR}/electricassemblyjacobian_inst.cpp
      ${ELECTRIC_LIB_DIR}/electricassemblyresidual_inst.cpp
      )
    set(ELECTRIC_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials feelpp_toolbox_electricbase )
    # generate the lib target
    genLibBase(
      LIB_NAME ${ELECTRIC_LIB_NAME}
      LIB_DIR ${ELECTRIC_LIB_DIR}
      LIB_DEPENDS ${ELECTRIC_LIB_DEPENDS}
      FILES_TO_COPY ${ELECTRIC_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${ELECTRIC_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/electric/electricconfig.h.in
      )
  endif()
endmacro(genLibElectric)

#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro( genLibSolidMechanics )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;DISP_ORDER;GEO_ORDER"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_DISP_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_DISP_ORDER OR FEELMODELS_APP_GEO_ORDER")
  endif()

  set(SOLIDMECHANICS_DIM ${FEELMODELS_APP_DIM})
  set(SOLIDMECHANICS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
  set(SOLIDMECHANICS_ORDER_DISPLACEMENT ${FEELMODELS_APP_DISP_ORDER})

  set(SOLIDMECHANICS_LIB_VARIANTS ${SOLIDMECHANICS_DIM}dP${SOLIDMECHANICS_ORDER_DISPLACEMENT}G${SOLIDMECHANICS_ORDERGEO})
  set(SOLIDMECHANICS_LIB_NAME feelpp_toolbox_solid_lib_${SOLIDMECHANICS_LIB_VARIANTS})

  if ( NOT TARGET ${SOLIDMECHANICS_LIB_NAME} )
    # configure the lib
    set(SOLIDMECHANICS_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/solid/${SOLIDMECHANICS_LIB_VARIANTS})
    set(SOLIDMECHANICS_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/solid/solidmechanicscreate_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/solid/solidmechanicsothers_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/solid/solidmechanicsupdatelinear_inst.cpp
      #${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/solid/solidmechanicsupdatelinear1dreduced_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/solid/solidmechanics1dreduced_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/solid/solidmechanicsupdatejacobian_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/solid/solidmechanicsupdateresidual_inst.cpp
      )
    set(SOLIDMECHANICS_CODEGEN_SOURCES
      ${SOLIDMECHANICS_LIB_DIR}/solidmechanicscreate_inst.cpp
      ${SOLIDMECHANICS_LIB_DIR}/solidmechanicsothers_inst.cpp
      ${SOLIDMECHANICS_LIB_DIR}/solidmechanicsupdatelinear_inst.cpp
      #${SOLIDMECHANICS_LIB_DIR}/solidmechanicsupdatelinear1dreduced_inst.cpp
      ${SOLIDMECHANICS_LIB_DIR}/solidmechanics1dreduced_inst.cpp
      ${SOLIDMECHANICS_LIB_DIR}/solidmechanicsupdatejacobian_inst.cpp
      ${SOLIDMECHANICS_LIB_DIR}/solidmechanicsupdateresidual_inst.cpp
      )
    set(SOLIDMECHANICS_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials feelpp_toolbox_solidbase ) 
    # generate the lib target
    genLibBase(
      LIB_NAME ${SOLIDMECHANICS_LIB_NAME}
      LIB_DIR ${SOLIDMECHANICS_LIB_DIR}
      LIB_DEPENDS ${SOLIDMECHANICS_LIB_DEPENDS}
      FILES_TO_COPY ${SOLIDMECHANICS_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${SOLIDMECHANICS_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/solid/solidmechanicsconfig.h.in
      )
  endif()

endmacro( genLibSolidMechanics )

#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro(genLibFluidMechanics)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;U_ORDER;P_ORDER;P_CONTINUITY;GEO_ORDER"
    ""
    ${ARGN}
    )

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


  set(FLUIDMECHANICS_LIB_VARIANTS ${FLUIDMECHANICS_DIM}dP${FLUIDMECHANICS_ORDER_VELOCITY}P${FLUIDMECHANICS_ORDER_PRESSURE}${FLUIDMECHANICS_PRESSURE_CONTINUITY_TAG}G${FLUIDMECHANICS_ORDERGEO})
  set(FLUIDMECHANICS_LIB_NAME feelpp_toolbox_fluid_lib_${FLUIDMECHANICS_LIB_VARIANTS})

  if ( NOT TARGET ${FLUIDMECHANICS_LIB_NAME} )
    # configure the lib
    set(FLUIDMECHANICS_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/fluid/${FLUIDMECHANICS_LIB_VARIANTS})
    set(FLUIDMECHANICS_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicscreate_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicsothers_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicsassemblylinear_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicsassemblyjacobian_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicsassemblyresidual_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicsupdatestabilisation_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicsassemblyturbulence_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicsmaterialproperties.cpp
      )
    set(FLUIDMECHANICS_CODEGEN_SOURCES
      ${FLUIDMECHANICS_LIB_DIR}/fluidmechanicscreate_inst.cpp
      ${FLUIDMECHANICS_LIB_DIR}/fluidmechanicsothers_inst.cpp
      ${FLUIDMECHANICS_LIB_DIR}/fluidmechanicsassemblylinear_inst.cpp
      ${FLUIDMECHANICS_LIB_DIR}/fluidmechanicsassemblyjacobian_inst.cpp
      ${FLUIDMECHANICS_LIB_DIR}/fluidmechanicsassemblyresidual_inst.cpp
      ${FLUIDMECHANICS_LIB_DIR}/fluidmechanicsupdatestabilisation_inst.cpp
      ${FLUIDMECHANICS_LIB_DIR}/fluidmechanicsassemblyturbulence_inst.cpp
      ${FLUIDMECHANICS_LIB_DIR}/fluidmechanicsmaterialproperties.cpp
      )
    set(FLUIDMECHANICS_LIB_DEPENDS feelpp_toolbox_fluidbase feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials feelpp_toolbox_coefficientformpdes_${FLUIDMECHANICS_DIM}dG${FLUIDMECHANICS_ORDERGEO} )
    # if ( FEELPP_TOOLBOXES_ENABLE_MESHALE )
    #   set(FLUIDMECHANICS_LIB_DEPENDS feelpp_modelmeshale ${FLUIDMECHANICS_LIB_DEPENDS})
    # endif()

    # generate the lib target
    genLibBase(
      LIB_NAME ${FLUIDMECHANICS_LIB_NAME}
      LIB_DIR ${FLUIDMECHANICS_LIB_DIR}
      LIB_DEPENDS ${FLUIDMECHANICS_LIB_DEPENDS}
      FILES_TO_COPY ${FLUIDMECHANICS_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${FLUIDMECHANICS_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fluid/fluidmechanicsconfig.h.in
      )

  endif()

endmacro( genLibFluidMechanics )

#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro(genLibFSI)
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;BC_MARKERS;FLUID_U_ORDER;FLUID_P_ORDER;FLUID_P_CONTINUITY;FLUID_GEO_ORDER;FLUID_GEO_DESC;FLUID_BC_DESC;SOLID_DISP_ORDER;SOLID_GEO_ORDER;SOLID_BC_DESC;SOLID_GEO_DESC"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_FLUID_GEO_ORDER OR FEELMODELS_APP_FLUID_U_ORDER OR FEELMODELS_APP_FLUID_P_ORDER OR
        FEELMODELS_APP_SOLID_DISP_ORDER OR FEELMODELS_APP_SOLID_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument!")
  endif()

  # fluid lib
  genLibFluidMechanics(
    DIM          ${FEELMODELS_APP_DIM}
    GEO_ORDER    ${FEELMODELS_APP_FLUID_GEO_ORDER}
    U_ORDER      ${FEELMODELS_APP_FLUID_U_ORDER}
    P_ORDER      ${FEELMODELS_APP_FLUID_P_ORDER}
    P_CONTINUITY ${FEELMODELS_APP_FLUID_P_CONTINUITY}
    )
  # solid lib
  genLibSolidMechanics(
    DIM ${FEELMODELS_APP_DIM}
    DISP_ORDER ${FEELMODELS_APP_SOLID_DISP_ORDER}
    GEO_ORDER ${FEELMODELS_APP_SOLID_GEO_ORDER}
    )

  set(FSI_LIB_VARIANTS ${FLUIDMECHANICS_LIB_VARIANTS}_${SOLIDMECHANICS_LIB_VARIANTS})
  set(FSI_LIB_NAME feelpp_toolbox_fsi_lib_${FSI_LIB_VARIANTS})

  if ( NOT TARGET ${FSI_LIB_NAME} )
    # configure the lib
    set(FSI_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/fsi/${FSI_LIB_VARIANTS})
    set(FSI_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fsi/fsi_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fsi/fsibc_inst.cpp
      )
    set(FSI_CODEGEN_SOURCES
      ${FSI_LIB_DIR}/fsi_inst.cpp
      ${FSI_LIB_DIR}/fsibc_inst.cpp
      )
    set(FSI_LIB_DEPENDS ${FLUIDMECHANICS_LIB_NAME} ${SOLIDMECHANICS_LIB_NAME})
    # generate the lib target
    genLibBase(
      LIB_NAME ${FSI_LIB_NAME}
      LIB_DIR ${FSI_LIB_DIR}
      LIB_DEPENDS ${FSI_LIB_DEPENDS}
      FILES_TO_COPY ${FSI_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${FSI_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/fsi/fsiconfig.h.in
      )
  endif()

endmacro( genLibFSI )

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro( genLibCoefficientFormPDE )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;UNKNOWN_BASIS_TYPE;UNKNOWN_BASIS_TAG;GEO_ORDER"
    ""
    ${ARGN}
    )

    if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_UNKNOWN_BASIS_TYPE OR FEELMODELS_APP_UNKNOWN_BASIS_TAG OR  FEELMODELS_APP_GEO_ORDER ) )
        message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR  FEELMODELS_APP_UNKNOWN_BASIS_TYPE OR FEELMODELS_APP_UNKNOWN_BASIS_TAG OR  FEELMODELS_APP_GEO_ORDER")
    endif()

  set(COEFFICIENTFORMPDE_DIM ${FEELMODELS_APP_DIM})
  set(COEFFICIENTFORMPDE_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
  set(COEFFICIENTFORMPDE_GEOSHAPE Simplex<${COEFFICIENTFORMPDE_DIM},${COEFFICIENTFORMPDE_ORDERGEO},${COEFFICIENTFORMPDE_DIM}>)
  set(COEFFICIENTFORMPDE_UNKNOWN_BASIS_TYPE ${FEELMODELS_APP_UNKNOWN_BASIS_TYPE})
  set(COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG ${FEELMODELS_APP_UNKNOWN_BASIS_TAG})

  set(COEFFICIENTFORMPDE_LIB_VARIANTS ${COEFFICIENTFORMPDE_DIM}d${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG}G${COEFFICIENTFORMPDE_ORDERGEO})
  set(COEFFICIENTFORMPDE_LIB_NAME feelpp_toolbox_coefficientformpde_${COEFFICIENTFORMPDE_LIB_VARIANTS})

  if ( NOT TARGET ${COEFFICIENTFORMPDE_LIB_NAME} )
      # configure the lib
      set(COEFFICIENTFORMPDE_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/coefficientformpdes/pde_${COEFFICIENTFORMPDE_LIB_VARIANTS})
      set(COEFFICIENTFORMPDE_CODEGEN_FILES_TO_COPY
          ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpde_inst.cpp
          )
      set(COEFFICIENTFORMPDE_CODEGEN_SOURCES
          ${COEFFICIENTFORMPDE_LIB_DIR}/coefficientformpde_inst.cpp
          )
      set(COEFFICIENTFORMPDE_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials feelpp_toolbox_coefficientformpdebase  ) 
      # generate the lib target
      genLibBase(
          LIB_NAME ${COEFFICIENTFORMPDE_LIB_NAME}
          LIB_DIR ${COEFFICIENTFORMPDE_LIB_DIR}
          LIB_DEPENDS ${COEFFICIENTFORMPDE_LIB_DEPENDS}
          FILES_TO_COPY ${COEFFICIENTFORMPDE_CODEGEN_FILES_TO_COPY}
          FILES_SOURCES ${COEFFICIENTFORMPDE_CODEGEN_SOURCES}
          CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpdeconfig.h.in
          )
  endif()

endmacro(genLibCoefficientFormPDE)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro( genLibCoefficientFormPDEs )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;UNKNOWN_BASIS_TYPE;UNKNOWN_BASIS_TAG;GEO_ORDER"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_UNKNOWN_BASIS_TYPE OR FEELMODELS_APP_UNKNOWN_BASIS_TAG OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR  FEELMODELS_APP_UNKNOWN_BASIS_TYPE OR FEELMODELS_APP_UNKNOWN_BASIS_TAG OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(COEFFICIENTFORMPDES_DIM ${FEELMODELS_APP_DIM})
  set(COEFFICIENTFORMPDES_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
  set(COEFFICIENTFORMPDES_GEOSHAPE Simplex<${COEFFICIENTFORMPDES_DIM},${COEFFICIENTFORMPDES_ORDERGEO},${COEFFICIENTFORMPDES_DIM}>)
  set(COEFFICIENTFORMPDES_UNKNOWN_BASIS_TYPE ${FEELMODELS_APP_UNKNOWN_BASIS_TYPE})
  set(COEFFICIENTFORMPDES_UNKNOWN_BASIS_TAG ${FEELMODELS_APP_UNKNOWN_BASIS_TAG})

  list(LENGTH COEFFICIENTFORMPDES_UNKNOWN_BASIS_TYPE count)
  list(LENGTH COEFFICIENTFORMPDES_UNKNOWN_BASIS_TAG count2)
  if ( NOT (count EQUAL count2) )
    message( FATAL_ERROR "UNKNOWN_BASIS_TYPE and UNKNOWN_BASIS_TAG should be same size" )
  endif()

  unset( COEFFICIENTFORMPDES_LIB_DEPENDS )
  math(EXPR count "${count}-1")
  foreach(i RANGE ${count})
    list(GET COEFFICIENTFORMPDES_UNKNOWN_BASIS_TYPE ${i} COEFFICIENTFORMPDE_UNKNOWN_BASIS_TYPE)
    list(GET COEFFICIENTFORMPDES_UNKNOWN_BASIS_TAG ${i} COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG)

    if ( i EQUAL 0 )
      set( COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TYPE "${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TYPE}")
      set( COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TAG "\"${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG}\"")
    else()
      set( COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TYPE "${COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TYPE} , ${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TYPE}")
      set( COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TAG "${COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TAG} , \"${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG}\"")
    endif()
    #set(COEFFICIENTFORMPDE_LIB_VARIANTS ${COEFFICIENTFORMPDES_DIM}d${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG}G${COEFFICIENTFORMPDES_ORDERGEO} )
    #set(COEFFICIENTFORMPDE_LIB_NAME feelpp_toolbox_coefficientformpde_${COEFFICIENTFORMPDE_LIB_VARIANTS})

    #if ( NOT TARGET ${COEFFICIENTFORMPDE_LIB_NAME} )
      ## configure the lib
      #set(COEFFICIENTFORMPDE_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/coefficientformpdes/pde_${COEFFICIENTFORMPDE_LIB_VARIANTS})
      #set(COEFFICIENTFORMPDE_CODEGEN_FILES_TO_COPY
        #${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpde_inst.cpp
        #)
      #set(COEFFICIENTFORMPDE_CODEGEN_SOURCES
        #${COEFFICIENTFORMPDE_LIB_DIR}/coefficientformpde_inst.cpp
        #)
      #set(COEFFICIENTFORMPDE_LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore feelpp_toolbox_coefficientformpdebase  ) 
      ## generate the lib target
      #genLibBase(
        #LIB_NAME ${COEFFICIENTFORMPDE_LIB_NAME}
        #LIB_DIR ${COEFFICIENTFORMPDE_LIB_DIR}
        #LIB_DEPENDS ${COEFFICIENTFORMPDE_LIB_DEPENDS}
        #FILES_TO_COPY ${COEFFICIENTFORMPDE_CODEGEN_FILES_TO_COPY}
        #FILES_SOURCES ${COEFFICIENTFORMPDE_CODEGEN_SOURCES}
        #CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpdeconfig.h.in
        #)
    #endif()
    genLibCoefficientFormPDE(
        DIM                     ${COEFFICIENTFORMPDES_DIM}
        UNKNOWN_BASIS_TYPE      ${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TYPE}
        UNKNOWN_BASIS_TAG       ${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG}
        GEO_ORDER               ${COEFFICIENTFORMPDES_ORDERGEO}
        )
    set(COEFFICIENTFORMPDES_LIB_DEPENDS ${COEFFICIENTFORMPDES_LIB_DEPENDS} ${COEFFICIENTFORMPDE_LIB_NAME})
  endforeach()

  set(COEFFICIENTFORMPDES_LIB_VARIANTS ${COEFFICIENTFORMPDES_DIM}dG${COEFFICIENTFORMPDES_ORDERGEO} )
  set(COEFFICIENTFORMPDES_LIB_NAME feelpp_toolbox_coefficientformpdes_${COEFFICIENTFORMPDES_LIB_VARIANTS})

  if ( NOT TARGET ${COEFFICIENTFORMPDES_LIB_NAME} )
    # configure the lib
    set(COEFFICIENTFORMPDES_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/coefficientformpdes/pdes_${COEFFICIENTFORMPDES_LIB_VARIANTS})
    set(COEFFICIENTFORMPDES_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpdes_inst.cpp
      )
    set(COEFFICIENTFORMPDES_CODEGEN_SOURCES
      ${COEFFICIENTFORMPDES_LIB_DIR}/coefficientformpdes_inst.cpp
      )

    set(COEFFICIENTFORMPDES_TARGET_COPY_FILES ${COEFFICIENTFORMPDES_LIB_NAME}_copyfiles)
    add_custom_target(${COEFFICIENTFORMPDES_TARGET_COPY_FILES}  ALL COMMENT "Copying modified files"  )

    # specialisation
    foreach(i RANGE ${count})
      list(GET COEFFICIENTFORMPDES_UNKNOWN_BASIS_TYPE ${i} COEFFICIENTFORMPDE_UNKNOWN_BASIS_TYPE)
      list(GET COEFFICIENTFORMPDES_UNKNOWN_BASIS_TAG ${i} COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG)

      set( COEFFICIENTFORMPDES_UNKNOWN_BASIS_SPECIALISATION ${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TYPE} )
      set( COEFFICIENTFORMPDES_LIB_SPECIALISATION_DIR  ${COEFFICIENTFORMPDES_LIB_DIR}/${COEFFICIENTFORMPDE_UNKNOWN_BASIS_TAG} )

      set(FEELMODELS_GENLIB_CONFIG_BASISSPEC_PATH  ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpdesbasisspecialisation.h.in)
      get_filename_component(FEELMODELS_GENLIB_CONFIG_BASISSPEC_FILENAME_WE ${FEELMODELS_GENLIB_CONFIG_BASISSPEC_PATH} NAME_WE)
      CONFIGURE_FILE( ${FEELMODELS_GENLIB_CONFIG_BASISSPEC_PATH} ${COEFFICIENTFORMPDES_LIB_SPECIALISATION_DIR}/${FEELMODELS_GENLIB_CONFIG_BASISSPEC_FILENAME_WE}.h )


      set( COEFFICIENTFORMPDES_SPECIALISATION_CODEGEN_FILES_TO_COPY
        ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpdesassemblylinear_spec.cpp
        ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpdesassemblyjacobian_spec.cpp
        ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpdesassemblyresidual_spec.cpp
        )
      # lib files
      foreach(filepath ${COEFFICIENTFORMPDES_SPECIALISATION_CODEGEN_FILES_TO_COPY})
        get_filename_component(filename ${filepath} NAME)
        if ( NOT EXISTS ${COEFFICIENTFORMPDES_LIB_SPECIALISATION_DIR}/${filename} )
          file(WRITE ${COEFFICIENTFORMPDES_LIB_SPECIALISATION_DIR}/${filename} "") #write empty file
        endif()
        add_custom_command(TARGET ${COEFFICIENTFORMPDES_TARGET_COPY_FILES} COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${filepath} ${COEFFICIENTFORMPDES_LIB_SPECIALISATION_DIR}/${filename} )

        set(COEFFICIENTFORMPDES_CODEGEN_SOURCES ${COEFFICIENTFORMPDES_CODEGEN_SOURCES} ${COEFFICIENTFORMPDES_LIB_SPECIALISATION_DIR}/${filename} )
      endforeach()

    endforeach()


    # generate the lib target
    genLibBase(
      LIB_NAME ${COEFFICIENTFORMPDES_LIB_NAME}
      LIB_DIR ${COEFFICIENTFORMPDES_LIB_DIR}
      LIB_DEPENDS ${COEFFICIENTFORMPDES_LIB_DEPENDS}
      FILES_TO_COPY ${COEFFICIENTFORMPDES_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${COEFFICIENTFORMPDES_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/coefficientformpdes/coefficientformpdesconfig.h.in
      TARGET_COPY_FILES ${COEFFICIENTFORMPDES_TARGET_COPY_FILES}
      )

    set(FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ENTRY_CLASS_TYPE "boost::mpl::pair< ${COEFFICIENTFORMPDES_GEOSHAPE}, Feel::FeelModels::CoefficientFormPDEs< ${COEFFICIENTFORMPDES_GEOSHAPE}, ${COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TYPE} > >" )
    
  endif()



endmacro(genLibCoefficientFormPDEs)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro( genLibAdvection )
  PARSE_ARGUMENTS(FEELMODELS_APP
      "DIM;POLY_ORDER;CONTINUITY;POLY_SET;GEO_ORDER;DIFFUSION_REACTION_ORDER;DIFFUSIONCOEFF_POLY_SET;REACTIONCOEFF_POLY_SET;DIFFUSION_REACTION_CONTINUITY"
      ""
      ${ARGN}
    )

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

  set(ADVECTION_LIB_VARIANTS ${ADVECTION_DIM}dP${ADVECTION_ORDERPOLY}${ADVECTION_CONTINUITY_TAG}${ADVECTION_POLY_SET_TAG}G${ADVECTION_ORDERGEO}${ADVECTION_DIFFUSION_REACTION_TAG} )
  set(ADVECTION_LIB_NAME feelpp_toolbox_advection_lib_${ADVECTION_LIB_VARIANTS})

  if ( NOT TARGET ${ADVECTION_LIB_NAME} )
    # configure the lib
    set(ADVECTION_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/advection/${ADVECTION_LIB_VARIANTS})    
    set(ADVECTION_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/advection/advection_inst.cpp )
    set(ADVECTION_CODEGEN_SOURCES
      ${ADVECTION_LIB_DIR}/advection_inst.cpp )
    set(ADVECTION_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials ) 
    # generate the lib target
    genLibBase(
      LIB_NAME ${ADVECTION_LIB_NAME}
      LIB_DIR ${ADVECTION_LIB_DIR}
      LIB_DEPENDS ${ADVECTION_LIB_DEPENDS}
      FILES_TO_COPY ${ADVECTION_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${ADVECTION_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/advection/advectionconfig.h.in
      )
  endif()


endmacro(genLibAdvection)
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibLevelsetBase )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;LEVELSET_ORDER;LEVELSET_PN_ORDER;GEO_ORDER"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_LEVELSET_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_LEVELSET_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(LEVELSETBASE_DIM ${FEELMODELS_APP_DIM})
  set(LEVELSETBASE_ORDERPOLY ${FEELMODELS_APP_LEVELSET_ORDER})
  set(LEVELSETBASE_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
  #######################################################
  if(FEELMODELS_APP_LEVELSET_PN_ORDER)
      set(LEVELSETBASE_PN_ORDERPOLY ${FEELMODELS_APP_LEVELSET_PN_ORDER})
  else()
      set(LEVELSETBASE_PN_ORDERPOLY ${FEELMODELS_APP_LEVELSET_ORDER})
  endif()
  #######################################################

  set(LEVELSETBASE_LIB_VARIANTS ${LEVELSETBASE_DIM}dP${LEVELSETBASE_ORDERPOLY}G${LEVELSETBASE_ORDERGEO} )
  set(LEVELSETBASE_LIB_NAME feelpp_toolbox_levelsetbase_lib_${LEVELSETBASE_LIB_VARIANTS})

  if ( NOT TARGET ${LEVELSETBASE_LIB_NAME} )
      # configure the lib
      set(LEVELSETBASE_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/levelset/${LEVELSETBASE_LIB_VARIANTS})
      set(LEVELSETBASE_CODEGEN_FILES_TO_COPY
          ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/levelsetbase_inst.cpp
          ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/levelsetspacemanager_inst.cpp
          ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/levelsetredistanciation_hj_inst.cpp
          ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/levelsetredistanciation_fm_inst.cpp
          ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/parameter_map.cpp )
      set(LEVELSETBASE_CODEGEN_SOURCES
          ${LEVELSETBASE_LIB_DIR}/levelsetbase_inst.cpp
          ${LEVELSETBASE_LIB_DIR}/levelsetspacemanager_inst.cpp
          ${LEVELSETBASE_LIB_DIR}/levelsetredistanciation_hj_inst.cpp
          ${LEVELSETBASE_LIB_DIR}/levelsetredistanciation_fm_inst.cpp
          ${LEVELSETBASE_LIB_DIR}/parameter_map.cpp )
      set(LEVELSETBASE_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore )
      # generate the lib target
      genLibBase(
          LIB_NAME ${LEVELSETBASE_LIB_NAME}
          LIB_DIR ${LEVELSETBASE_LIB_DIR}
          LIB_DEPENDS ${LEVELSETBASE_LIB_DEPENDS}
          FILES_TO_COPY ${LEVELSETBASE_CODEGEN_FILES_TO_COPY}
          FILES_SOURCES ${LEVELSETBASE_CODEGEN_SOURCES}
          CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/levelsetbaseconfig.h.in
          )
  endif()
endmacro(genLibLevelsetBase)

#############################################################################
#############################################################################
#############################################################################
macro( genLibLevelset )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;LEVELSET_ORDER;LEVELSET_PN_ORDER;GEO_ORDER"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_LEVELSET_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_LEVELSET_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(LEVELSET_DIM ${FEELMODELS_APP_DIM})
  set(LEVELSET_ORDERPOLY ${FEELMODELS_APP_LEVELSET_ORDER})
  set(LEVELSET_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
  #######################################################
  if(FEELMODELS_APP_LEVELSET_PN_ORDER)
      set(LEVELSET_PN_ORDERPOLY ${FEELMODELS_APP_LEVELSET_PN_ORDER})
  else()
      set(LEVELSET_PN_ORDERPOLY ${FEELMODELS_APP_LEVELSET_ORDER})
  endif()
  #######################################################
  # set levelset lib name
  set(LEVELSET_LIB_VARIANTS ${LEVELSET_DIM}dP${LEVELSET_ORDERPOLY}G${LEVELSET_ORDERGEO} )
  set(LEVELSET_LIB_NAME feelpp_toolbox_levelset_lib_${LEVELSET_LIB_VARIANTS})

  unset(LEVELSET_LIB_DEPENDS)
  #######################################################
  # generate levelsetbase target
  genLibLevelsetBase(
      DIM               ${FEELMODELS_APP_DIM}
      LEVELSET_ORDER    ${FEELMODELS_APP_LEVELSET_ORDER}
      LEVELSET_PN_ORDER ${FEELMODELS_APP_LEVELSET_PN_ORDER}
      GEO_ORDER         ${FEELMODELS_APP_GEO_ORDER}
      )
  set(LEVELSET_LIB_DEPENDS ${LEVELSET_LIB_DEPENDS} ${LEVELSETBASE_LIB_NAME})
  #######################################################
  # generate coefficientformpde targets (scalar+vectorial) (dependency of levelset lib)
  set( SCALAR_UNKNOWN_BASIS_TYPE Lagrange<${LEVELSET_ORDERPOLY},Scalar,Continuous,PointSetFekete> )
  set( SCALAR_UNKNOWN_BASIS_TAG Pch${LEVELSET_ORDERPOLY} )
  genLibCoefficientFormPDE(
      DIM                       ${FEELMODELS_APP_DIM}
      UNKNOWN_BASIS_TYPE        ${SCALAR_UNKNOWN_BASIS_TYPE}
      UNKNOWN_BASIS_TAG         ${SCALAR_UNKNOWN_BASIS_TAG}
      GEO_ORDER                 ${LEVELSET_ORDERGEO}
      )
  set(LEVELSET_LIB_DEPENDS ${LEVELSET_LIB_DEPENDS} ${COEFFICIENTFORMPDE_LIB_NAME})

  set( VECTORIAL_UNKNOWN_BASIS_TYPE Lagrange<${LEVELSET_ORDERPOLY},Vectorial,Continuous,PointSetFekete> )
  set( VECTORIAL_UNKNOWN_BASIS_TAG Pchv${LEVELSET_ORDERPOLY} )
  genLibCoefficientFormPDE(
      DIM                       ${FEELMODELS_APP_DIM}
      UNKNOWN_BASIS_TYPE        ${VECTORIAL_UNKNOWN_BASIS_TYPE}
      UNKNOWN_BASIS_TAG         ${VECTORIAL_UNKNOWN_BASIS_TAG}
      GEO_ORDER                 ${LEVELSET_ORDERGEO}
      )
  set(LEVELSET_LIB_DEPENDS ${LEVELSET_LIB_DEPENDS} ${COEFFICIENTFORMPDE_LIB_NAME})

  #######################################################
  # generate levelset target
  if ( NOT TARGET ${LEVELSET_LIB_NAME} )
      # configure the lib
      set(LEVELSET_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/levelset/${LEVELSET_LIB_VARIANTS})
      set(LEVELSET_CODEGEN_FILES_TO_COPY
          ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/levelset_inst.cpp
          #${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/levelsetadvection_inst.cpp
          ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/parameter_map.cpp )
      set(LEVELSET_CODEGEN_SOURCES
          ${LEVELSET_LIB_DIR}/levelset_inst.cpp
          #${LEVELSET_LIB_DIR}/levelsetadvection_inst.cpp
          ${LEVELSET_LIB_DIR}/parameter_map.cpp )
      #set(LEVELSET_LIB_DEPENDS feelpp_modelalg feelpp_modelmesh feelpp_modelcore ${LEVELSETBASE_LIB_NAME} feelpp_toolbox_coefficientformpde_2dPch1G1 )
      # generate the lib target
      genLibBase(
          LIB_NAME ${LEVELSET_LIB_NAME}
          LIB_DIR ${LEVELSET_LIB_DIR}
          LIB_DEPENDS ${LEVELSET_LIB_DEPENDS}
          FILES_TO_COPY ${LEVELSET_CODEGEN_FILES_TO_COPY}
          FILES_SOURCES ${LEVELSET_CODEGEN_SOURCES}
          CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/levelset/levelsetconfig.h.in
          )
  endif()

endmacro(genLibLevelset)
#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro(genLibMultiFluid)
    PARSE_ARGUMENTS(FEELMODELS_APP
        "DIM;FLUID_U_ORDER;FLUID_P_ORDER;FLUID_P_CONTINUITY;GEO_ORDER;LEVELSET_ORDER;LEVELSET_PN_ORDER"
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
    )
  ###############################################################
  # levelset lib
  genLibLevelset(
    DIM ${FEELMODELS_APP_DIM}
    LEVELSET_ORDER ${FEELMODELS_APP_LEVELSET_ORDER}
    LEVELSET_PN_ORDER ${FEELMODELS_APP_LEVELSET_PN_ORDER}
    GEO_ORDER ${FEELMODELS_APP_GEO_ORDER}
    )
  ###############################################################
  # multifluid lib
  set(MULTIFLUID_LIB_VARIANTS ${FLUIDMECHANICS_LIB_VARIANTS}_${LEVELSET_LIB_VARIANTS}  )
  set(MULTIFLUID_LIB_NAME feelpp_toolbox_multifluid_lib_${MULTIFLUID_LIB_VARIANTS})

  if ( NOT TARGET ${MULTIFLUID_LIB_NAME} )
    # configure the lib
    set(MULTIFLUID_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/multifluid/${MULTIFLUID_LIB_VARIANTS})
    set(MULTIFLUID_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/multifluid/multifluid_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/multifluid/multifluidassemblylinear_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/multifluid/multifluidassemblyjacobian_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/multifluid/multifluidassemblyresidual_inst.cpp
      )
    set(MULTIFLUID_CODEGEN_SOURCES
      ${MULTIFLUID_LIB_DIR}/multifluid_inst.cpp
      ${MULTIFLUID_LIB_DIR}/multifluidassemblylinear_inst.cpp
      ${MULTIFLUID_LIB_DIR}/multifluidassemblyjacobian_inst.cpp
      ${MULTIFLUID_LIB_DIR}/multifluidassemblyresidual_inst.cpp
      )
    set(MULTIFLUID_LIB_DEPENDS ${FLUIDMECHANICS_LIB_NAME} ${LEVELSET_LIB_NAME})
    # generate the lib target
    genLibBase(
      LIB_NAME ${MULTIFLUID_LIB_NAME}
      LIB_DIR ${MULTIFLUID_LIB_DIR}
      LIB_DEPENDS ${MULTIFLUID_LIB_DEPENDS}
      FILES_TO_COPY ${MULTIFLUID_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${MULTIFLUID_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/multifluid/multifluidconfig.h.in
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
    "DIM;T_ORDER;P_ORDER;GEO_ORDER;"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_T_ORDER OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_GEO_ORDER")
  endif()

  set(THERMOELECTRIC_DIM ${FEELMODELS_APP_DIM})
  set(THERMOELECTRIC_T_ORDER ${FEELMODELS_APP_T_ORDER})
  set(THERMOELECTRIC_P_ORDER ${FEELMODELS_APP_P_ORDER})
  set(THERMOELECTRIC_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  genLibHeat(
    DIM     ${THERMOELECTRIC_DIM}
    T_ORDER ${THERMOELECTRIC_T_ORDER}
    GEO_ORDER ${THERMOELECTRIC_ORDERGEO}
    )
  genLibElectric(
    DIM     ${THERMOELECTRIC_DIM}
    P_ORDER ${THERMOELECTRIC_P_ORDER}
    GEO_ORDER ${THERMOELECTRIC_ORDERGEO}
    )

  set(THERMOELECTRIC_LIB_VARIANTS ${HEAT_LIB_VARIANTS}_${ELECTRIC_LIB_VARIANTS})
  set(THERMOELECTRIC_LIB_NAME feelpp_toolbox_thermoelectric_lib_${THERMOELECTRIC_LIB_VARIANTS})

  if ( NOT TARGET ${THERMOELECTRIC_LIB_NAME} )
    # configure the lib
    set(THERMOELECTRIC_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/thermoelectric/${THERMOELECTRIC_LIB_VARIANTS})
    set(THERMOELECTRIC_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/thermoelectric/thermoelectric_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/thermoelectric/thermoelectricassemblylinear_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/thermoelectric/thermoelectricassemblyjacobian_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/thermoelectric/thermoelectricassemblyresidual_inst.cpp
      )
    set(THERMOELECTRIC_CODEGEN_SOURCES
      ${THERMOELECTRIC_LIB_DIR}/thermoelectric_inst.cpp
      ${THERMOELECTRIC_LIB_DIR}/thermoelectricassemblylinear_inst.cpp
      ${THERMOELECTRIC_LIB_DIR}/thermoelectricassemblyjacobian_inst.cpp
      ${THERMOELECTRIC_LIB_DIR}/thermoelectricassemblyresidual_inst.cpp
      )
    set(THERMOELECTRIC_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials)
    set(THERMOELECTRIC_LIB_DEPENDS ${HEAT_LIB_NAME} ${ELECTRIC_LIB_NAME} ${THERMOELECTRIC_LIB_DEPENDS} )
    # generate the lib target
    genLibBase(
      LIB_NAME ${THERMOELECTRIC_LIB_NAME}
      LIB_DIR ${THERMOELECTRIC_LIB_DIR}
      LIB_DEPENDS ${THERMOELECTRIC_LIB_DEPENDS}
      FILES_TO_COPY ${THERMOELECTRIC_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${THERMOELECTRIC_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/thermoelectric/thermoelectricconfig.h.in
      )
  endif()

endmacro(genLibThermoElectric)

#############################################################################
#############################################################################
#############################################################################
#############################################################################

macro( genLibHeatFluid )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;U_ORDER;P_ORDER;T_ORDER;GEO_ORDER;"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_T_ORDER OR FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_T_ORDER OR FEELMODELS_APP_GEO_ORDER")
  endif()

  set(HEATFLUID_DIM ${FEELMODELS_APP_DIM})
  set(HEATFLUID_U_ORDER ${FEELMODELS_APP_U_ORDER})
  set(HEATFLUID_P_ORDER ${FEELMODELS_APP_P_ORDER})
  set(HEATFLUID_T_ORDER ${FEELMODELS_APP_T_ORDER})
  set(HEATFLUID_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  genLibHeat(
    DIM     ${HEATFLUID_DIM}
    T_ORDER ${HEATFLUID_T_ORDER}
    GEO_ORDER ${HEATFLUID_ORDERGEO}
    )
  genLibFluidMechanics(
    DIM     ${HEATFLUID_DIM}
    U_ORDER ${HEATFLUID_U_ORDER}
    P_ORDER ${HEATFLUID_P_ORDER}
    GEO_ORDER ${HEATFLUID_ORDERGEO}
    )

  set(HEATFLUID_LIB_VARIANTS ${HEAT_LIB_VARIANTS}_${FLUIDMECHANICS_LIB_VARIANTS})
  set(HEATFLUID_LIB_NAME feelpp_toolbox_heatfluid_lib_${HEATFLUID_LIB_VARIANTS})

  if ( NOT TARGET ${HEATFLUID_LIB_NAME} )
    # configure the lib
    set(HEATFLUID_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/heatfluid/${HEATFLUID_LIB_VARIANTS})
    set(HEATFLUID_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heatfluid/heatfluid_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heatfluid/heatfluidassemblylinear_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heatfluid/heatfluidassemblyjacobian_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heatfluid/heatfluidassemblyresidual_inst.cpp
      )
    set(HEATFLUID_CODEGEN_SOURCES
      ${HEATFLUID_LIB_DIR}/heatfluid_inst.cpp
      ${HEATFLUID_LIB_DIR}/heatfluidassemblylinear_inst.cpp
      ${HEATFLUID_LIB_DIR}/heatfluidassemblyjacobian_inst.cpp
      ${HEATFLUID_LIB_DIR}/heatfluidassemblyresidual_inst.cpp
      )
    set(HEATFLUID_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials)
    set(HEATFLUID_LIB_DEPENDS ${HEAT_LIB_NAME} ${FLUIDMECHANICS_LIB_NAME} ${HEATFLUID_LIB_DEPENDS} )
    # generate the lib target
    genLibBase(
      LIB_NAME ${HEATFLUID_LIB_NAME}
      LIB_DIR ${HEATFLUID_LIB_DIR}
      LIB_DEPENDS ${HEATFLUID_LIB_DEPENDS}
      FILES_TO_COPY ${HEATFLUID_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${HEATFLUID_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/heatfluid/heatfluidconfig.h.in
      )
  endif()

endmacro(genLibHeatFluid)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibMaxwell )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;P_ORDER;GEO_ORDER;"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_P_ORDER OR  FEELMODELS_APP_GEO_ORDER ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_P_ORDER OR  FEELMODELS_APP_GEO_ORDER")
  endif()

  set(MAXWELL_DIM ${FEELMODELS_APP_DIM})
  set(MAXWELL_ORDERPOLY ${FEELMODELS_APP_P_ORDER})
  set(MAXWELL_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

  set(MAXWELL_LIB_VARIANTS ${MAXWELL_DIM}dP${MAXWELL_ORDERPOLY}G${MAXWELL_ORDERGEO} )
  set(MAXWELL_LIB_NAME feelpp_toolbox_maxwell_lib_${MAXWELL_LIB_VARIANTS})

  if ( NOT TARGET ${MAXWELL_LIB_NAME} )
    # configure the lib
    set(MAXWELL_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/maxwell/${MAXWELL_LIB_VARIANTS})
    set(MAXWELL_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/maxwell/maxwell_inst.cpp )
    set(MAXWELL_CODEGEN_SOURCES
      ${MAXWELL_LIB_DIR}/maxwell_inst.cpp )
    set(MAXWELL_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials)
    # generate the lib target
    genLibBase(
      LIB_NAME ${MAXWELL_LIB_NAME}
      LIB_DIR ${MAXWELL_LIB_DIR}
      LIB_DEPENDS ${MAXWELL_LIB_DEPENDS}
      FILES_TO_COPY ${MAXWELL_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${MAXWELL_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/maxwell/maxwellconfig.h.in
      )
  endif()
endmacro(genLibMaxwell)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
macro( genLibHdg )
  PARSE_ARGUMENTS(FEELMODELS_APP
    "DIM;P_ORDER;GEO_ORDER;POLYSET;"
    ""
    ${ARGN}
    )

  if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_GEO_ORDER OR FEELMODELS_APP_POLYSET ) )
    message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_GEO_ORDER OR FEELMODELS_APP_POLYSET")
  endif()

  set(HDG_DIM ${FEELMODELS_APP_DIM})
  set(HDG_ORDERPOLY ${FEELMODELS_APP_P_ORDER})
  set(HDG_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})
  set(HDG_POLYSET ${FEELMODELS_APP_POLYSET})

  set(HDG_LIB_VARIANTS ${HDG_DIM}dP${HDG_ORDERPOLY}G${HDG_ORDERGEO}_${HDG_POLYSET} )
  set(HDG_LIB_NAME feelpp_toolbox_hdg_lib_${HDG_LIB_VARIANTS})

  if ( NOT TARGET ${HDG_LIB_NAME} )
    # configure the lib
    set(HDG_LIB_DIR ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/hdg/${HDG_LIB_VARIANTS})
    set(HDG_CODEGEN_FILES_TO_COPY
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/hdg/mixedpoisson_inst.cpp
      ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/hdg/mixedpoissonassemblylinear_inst.cpp )
    set(HDG_CODEGEN_SOURCES
      ${HDG_LIB_DIR}/mixedpoisson_inst.cpp
      ${HDG_LIB_DIR}/mixedpoissonassemblylinear_inst.cpp )
    set(HDG_LIB_DEPENDS feelpp_modelmesh feelpp_modelcore feelpp_modelmaterials feelpp_toolbox_hdgbase )
    # generate the lib target
    genLibBase(
      LIB_NAME ${HDG_LIB_NAME}
      LIB_DIR ${HDG_LIB_DIR}
      LIB_DEPENDS ${HDG_LIB_DEPENDS}
      FILES_TO_COPY ${HDG_CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${HDG_CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/hdg/mixedpoissonconfig.h.in
      )
  endif()
endmacro(genLibHdg)

