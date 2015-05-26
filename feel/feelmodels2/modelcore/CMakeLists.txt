#############################################################################
option(FEELMODELS_ENABLE_FLUID_MECHANICS "Enable/Disable FluidMechanics" ON)
option(FEELMODELS_ENABLE_STOKES_VARIANTS "Enable/Disable StokesVariants" OFF)

option(FEELMODELS_FLUID_BUILD_LINEAR_CODE "Enable/Disable" ON)
option(FEELMODELS_FLUID_BUILD_JACOBIAN_CODE "Enable/Disable" ON)
option(FEELMODELS_FLUID_BUILD_RESIDUAL_CODE "Enable/Disable" ON)
option(FEELMODELS_FLUID_BUILD_PTFIXED_CODE "Enable/Disable" OFF)
option(FEELMODELS_FLUID_BUILD_STABILISATION_CODE "Enable/Disable" ON)

option(FEELMODELS_SOLID_BUILD_LINEAR_CODE "Enable/Disable" ON)
option(FEELMODELS_SOLID_BUILD_JACOBIAN_CODE "Enable/Disable" ON)
option(FEELMODELS_SOLID_BUILD_RESIDUAL_CODE "Enable/Disable" ON)

SET(FEELMODELS_FSIBASE_NUMBER_OF_PARAMETERS "10" CACHE STRING "number of parameters" )
SET(FEELMODELS_FSIBASE_NUMBER_OF_GEOPARAMETERS "10" CACHE STRING "number of geometric parameters" )
SET(FEELMODELS_FSIBASE_NUMBER_OF_GINACEXPR "10" CACHE STRING "number of ginac expr" )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/feelmodelscoreconstconfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/feelmodelscoreconstconfig.h )

set(FEELMODELS_LIBBASE_CHECK_DIR "${CMAKE_CURRENT_BINARY_DIR}/checkmodelsbuilt" CACHE STRING "FEELMODELS_LIBBASE_CHECK_DIR")
file(REMOVE_RECURSE ${FEELMODELS_LIBBASE_CHECK_DIR})


#############################################################################
#############################################################################

add_library(feelpp_modelcore
  SHARED
  options.cpp log.cpp timertool.cpp
  applibase.cpp applibasenumericalsimulationtransitory.cpp
  markermanagement.cpp
  )
target_link_libraries(feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES})


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
    "LIB_DIR;MARKERS;DESC;GEO;LIB_DEPENDS;PREFIX_INCLUDE_USERCONFIG;PREFIX_LIBNAME;FILES_TO_COPY;FILES_SOURCES;CONFIG_PATH"
    ""
    ${ARGN}
    )
  set( FEELMODELS_GENLIB_APPLICATION_DIR ${FEELMODELS_GENLIB_BASE_LIB_DIR} )
  #CAR(APPLICATION_NAME      ${FEELMODELS_GENLIB_BASE_DEFAULT_ARGS})
  #CDR(FEELMODELS_GENLIB_APPLICATION_DIR       ${FEELMODELS_GENLIB_BASE_DEFAULT_ARGS})

  set(LIB_DEPENDS           ${FEELMODELS_GENLIB_BASE_LIB_DEPENDS})
  set(PREFIX_FILES_TO_COPY  ${FEELMODELS_GENLIB_BASE_PREFIX_INCLUDE_USERCONFIG})
  set(LIB_APPLICATION_NAME ${FEELMODELS_GENLIB_BASE_PREFIX_LIBNAME})

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
      configure_file( ${BCMARKER_FILE} ${FEELMODELS_GENLIB_APPLICATION_DIR}/bcmarker.cpp COPYONLY)
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${BCMARKER_FILE}  ${FEELMODELS_GENLIB_APPLICATION_DIR}/bcmarker.cpp )
    # bctool
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/bctool.hpp )
      configure_file( ${FEELPP_MODELS_SOURCE_DIR}/modelcore/bctool.hpp ${FEELMODELS_GENLIB_APPLICATION_DIR}/bctool.hpp COPYONLY)
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${FEELPP_MODELS_SOURCE_DIR}/modelcore/bctool.hpp ${FEELMODELS_GENLIB_APPLICATION_DIR}/bctool.hpp )
  endif()
  # .bc
  if ( FEELMODELS_GENLIB_BASE_DESC )
    set(BCDESC_FILE ${FEELMODELS_GENLIB_BASE_DESC})
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.bc )
      configure_file( ${BCDESC_FILE} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.bc COPYONLY)
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${BCDESC_FILE}  ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.bc )
  endif()
  # .mesh
  if ( FEELMODELS_GENLIB_BASE_GEO )
    set(MESH_FILE ${FEELMODELS_GENLIB_BASE_GEO})
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.mesh )
      configure_file( ${MESH_FILE} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.mesh COPYONLY)
    endif()
    add_custom_command(TARGET codegen_${LIB_APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${MESH_FILE}  ${FEELMODELS_GENLIB_APPLICATION_DIR}/${PREFIX_FILES_TO_COPY}.mesh )
  endif()
  # lib files
  foreach(filepath ${CODEGEN_FILES_TO_COPY})
    get_filename_component(filename ${filepath} NAME)
    if ( NOT EXISTS ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} )
      configure_file( ${filepath} ${FEELMODELS_GENLIB_APPLICATION_DIR}/${filename} COPYONLY)
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
      configure_file( ${FEELPP_MODELS_SOURCE_DIR}/modelcore/${filename} ${APPLICATION_DIR}/${filename} COPYONLY)
    endforeach()
  endif()
  add_custom_target(codegen_env_${APPLICATION_NAME} ALL COMMENT "Copying modified files"  )
  foreach(filename applimanagement.hpp)
    add_custom_command(TARGET codegen_env_${APPLICATION_NAME} COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${FEELPP_MODELS_SOURCE_DIR}/modelcore/${filename} ${APPLICATION_DIR}/${filename} )
  endforeach()

  if ( NOT EXISTS ${APPLICATION_DIR}/${MAIN_FILE} )
    foreach(filename ${MAIN_FILE})
      configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${MAIN_FILE} ${APPLICATION_DIR}/${MAIN_FILE} COPYONLY)
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
  set(THERMODYNAMICS_ORDERGEO ${FEELMODELS_APP_T_ORDER})
  set(THERMODYNAMICS_ORDERPOLY ${FEELMODELS_APP_GEO_ORDER})

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
  set(LIBBASE_CHECK_PATH ${FEELMODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
  set(LIBBASE_NAME feelmodels_${FEELMODELS_MODEL_SPECIFIC_NAME})

  if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )

    #write empty file in orter to check if this lib has already define
    file(WRITE ${LIBBASE_CHECK_PATH} "")

    # configure libmodelbase
    # set(CODEGEN_FILES_TO_COPY
    #   ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodynbase.hpp
    #   ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodynbase.cpp
    #   ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodyn.hpp
    #   ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodyn.cpp )
    # set(CODEGEN_SOURCES
    #   ${LIBBASE_DIR}/thermodynbase.cpp
    #   ${LIBBASE_DIR}/thermodyn.cpp )
    set(CODEGEN_FILES_TO_COPY
      ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodynbase_inst.cpp
      ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodyn_inst.cpp )
    set(CODEGEN_SOURCES
      ${LIBBASE_DIR}/thermodynbase_inst.cpp
      ${LIBBASE_DIR}/thermodyn_inst.cpp )
    set(LIB_DEPENDS feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) 

    # generate libmodelbase
    genLibBase(
      LIB_DIR ${LIBBASE_DIR}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      PREFIX_LIBNAME ${LIBBASE_NAME}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/thermodyn/thermodynconfig.h.in
      )
  endif()

  ############################
  # build lib with bctool
  ############################
  if ( FEELMODELS_APP_BC_DESC )

    #set(BCDESC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_DESC})
    get_filename_component( BCDESC_FILE ${FEELMODELS_APP_BC_DESC} ABSOLUTE )

    if ( FEELMODELS_APP_BC_MARKERS )
      get_filename_component( BCMARKER_FILE ${FEELMODELS_APP_BC_MARKERS} ABSOLUTE )
      #set(BCMARKER_FILE ${FEELMODELS_APP_BC_MARKERS})
    else()
      set(BCMARKER_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/default.bcmarkers )
    endif()

    if ( FEELMODELS_APP_GEO_DESC )
      get_filename_component( GEO_FILE ${FEELMODELS_APP_GEO_DESC} ABSOLUTE )
      #set(GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GEO_DESC})
    else()
      set(GEO_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/default.mesh )
    endif()

    #set(PREFIX_LIB thermo)

    if ( NOT ( FEELMODELS_APP_LIB_NAME  ) )
      message(FATAL_ERROR "miss argument! FEELMODELS_APP_LIB_NAME")
    endif()
    set(LIB_NAME ${FEELMODELS_APP_LIB_NAME})

    if ( FEELMODELS_APP_LIB_DIR )
      get_filename_component( THE_LIB_DIR ${FEELMODELS_APP_LIB_DIR} ABSOLUTE )
      #set(THE_LIB_DIR ${FEELMODELS_APP_LIB_DIR})
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
      LIB_DIR ${THE_LIB_DIR}
      MARKERS ${BCMARKER_FILE}
      DESC ${BCDESC_FILE}
      GEO ${GEO_FILE}
      LIB_DEPENDS ${LIB_DEPENDS}
      PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
      PREFIX_LIBNAME ${LIB_NAME} #${PREFIX_LIB}${APPLICATION_NAME}
      FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
      FILES_SOURCES ${CODEGEN_SOURCES}
      CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/modelcore/feelmodelscoreconfig.h.in
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
    LIB_NAME   feelmodels_thermo${APPLICATION_NAME}
    LIB_DIR    ${THE_APPLICATION_DIR}/thermodyn
    DIM        ${FEELMODELS_DIM}
    T_ORDER    ${FEELMODELS_ORDERPOLY}
    GEO_ORDER  ${FEELMODELS_ORDERGEO}
    BC_DESC    ${BCDESC_FILE} 
    BC_MARKERS ${BCMARKER_FILE}
    GEO_DESC   ${GEO_FILE}
    )

  set(LIB_DEPENDS feelmodels_thermo${APPLICATION_NAME} feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ) #fsimesh 
  genExecutableBase(
    ${APPLICATION_NAME}
    APPLICATION_DIR ${THE_APPLICATION_DIR}
    SRC ${MAIN_FILE}
    LIB_DEPENDS ${LIB_DEPENDS}
    CONFIG_PATH ${FEELPP_MODELS_SOURCE_DIR}/modelcore/feelmodelscoreconfig.h.in
    )

endmacro(genExecutableThermoDynamics)
#############################################################################
#############################################################################
macro(fsi_add_thermo_application)
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
    set(BCMARKER_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/default.bcmarkers )
  endif()

  if ( FEELMODELS_APP_GEO_DESC )
    set(GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GEO_DESC})
  else()
    set(GEO_FILE ${FEELPP_MODELS_SOURCE_DIR}/modelcore/default.mesh )
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

endmacro(fsi_add_thermo_application)
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
# macro( genLibSolidMechanics APPLICATION_NAME APPLICATION_DIR BCMARKER_FILE FEELMODELS_DIM
#     FEELMODELS_ORDER_U FEELMODELS_ORDERGEO FEELMODELS_USE_CST_DENSITY_COEFFLAME BCDESC_FILE MESH_FILE )

#   set(SOLIDMECHANICS_DIM ${FEELMODELS_DIM})
#   set(SOLIDMECHANICS_ORDERGEO ${FEELMODELS_ORDERGEO})
#   set(SOLIDMECHANICS_ORDER_DISPLACEMENT ${FEELMODELS_ORDER_U})

#   if ("${FEELMODELS_USE_CST_DENSITY_COEFFLAME}" STREQUAL "1" ) 
#     set(SOLIDMECHANICS_USE_CST_DENSITY_COEFFLAME 1)
#     set(FEELMODELS_USE_CST_DENSITY_COEFFLAME_TAG DCLCst)
#   else()
#     set(SOLIDMECHANICS_USE_CST_DENSITY_COEFFLAME 0)
#     set(FEELMODELS_USE_CST_DENSITY_COEFFLAME_TAG DCLP0d)
#    endif()

#   set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${SOLIDMECHANICS_DIM}dP${SOLIDMECHANICS_ORDER_DISPLACEMENT}G${SOLIDMECHANICS_ORDERGEO}${FEELMODELS_USE_CST_DENSITY_COEFFLAME_TAG})
#   set(FEELMODELS_MODEL_SPECIFIC_NAME solidmec${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})
#   set(LIBBASE_DIR ${PROJECT_BINARY_DIR}/fsi/fsimodel/${FEELMODELS_MODEL_SPECIFIC_NAME} )

#   set(LIBBASE_CHECK_PATH ${FEELMODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
#   set(LIBBASE_NAME fsi_modelbase_${FEELMODELS_MODEL_SPECIFIC_NAME})

#   if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )
#     #write empty file in orter to check if this lib has already define
#     file(WRITE ${LIBBASE_CHECK_PATH} "")

#     set(CODEGEN_FILES_TO_COPY
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmecbase.hpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmecbasecreate.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmecbaseothers.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmecbaseupdateresidual.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmecbaseupdatejacobian.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmecbaseupdatelinear.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmecbaseupdatelinearaxi.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmecbaseupdatelinear1dreduced.cpp )
#     set(CODEGEN_SOURCES
#       ${LIBBASE_DIR}/solidmecbasecreate.cpp
#       ${LIBBASE_DIR}/solidmecbaseothers.cpp
#       ${LIBBASE_DIR}/solidmecbaseupdateresidual.cpp
#       ${LIBBASE_DIR}/solidmecbaseupdatejacobian.cpp
#       ${LIBBASE_DIR}/solidmecbaseupdatelinear.cpp
#       ${LIBBASE_DIR}/solidmecbaseupdatelinearaxi.cpp
#       ${LIBBASE_DIR}/solidmecbaseupdatelinear1dreduced.cpp
#       )
#     set(LIB_DEPENDS fsialg fsicore fsifilters ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} )

#     # generate libmodelbase
#     genLibBase( ${APPLICATION_NAME} ${LIBBASE_DIR}
#       LIB_DEPENDS ${LIB_DEPENDS}
#       PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
#       PREFIX_LIBNAME ${LIBBASE_NAME}
#       FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
#       FILES_SOURCES ${CODEGEN_SOURCES}
#       CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fsisolidmecconfig.h.in
#       )


#   endif()

#   set(PREFIX_FILES_TO_COPY solid)
#   set(PREFIX_LIB sm)
#   set(CODEGEN_FILES_TO_COPY
#     ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmec.hpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsimodel/solidmec.cpp
#     )
#   set(CODEGEN_SOURCES
#     ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}/solidmec.cpp
#     )
#   set(LIB_DEPENDS ${LIBBASE_NAME} fsialg fsicore ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} )

#   set(SOLIDMECHANICS_INCLUDE_LIBBASE_HEADER 1)

#   genLibBase( ${APPLICATION_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}
#     MARKERS ${BCMARKER_FILE}
#     DESC ${BCDESC_FILE}
#     GEO ${MESH_FILE}
#     LIB_DEPENDS ${LIB_DEPENDS}
#     PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
#     PREFIX_LIBNAME ${PREFIX_LIB}${APPLICATION_NAME}
#     FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
#     FILES_SOURCES ${CODEGEN_SOURCES}
#     CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in
#     )

#   unset(SOLIDMECHANICS_INCLUDE_LIBBASE_HEADER)

# endmacro( genLibSolidMechanics )
# #############################################################################
# #############################################################################
# macro(genExecutableSolidMechanics APPLICATION_NAME BCMARKER_FILE FEELMODELS_DIM
#     FEELMODELS_ORDER_U FEELMODELS_ORDERGEO FEELMODELS_USE_CST_DENSITY_COEFFLAME BCDESC_FILE MESH_FILE
#     MAIN_FILE ADDITIONAL_LIB_DEPEND )
#   resetToZeroAllPhysicalVariables()
#   SET(SOLIDMECHANICS 1 )

#   set(APPLICATION_DIR codeGen_${APPLICATION_NAME})
#   genLibSolidMechanics(${APPLICATION_NAME} ${APPLICATION_DIR}/solid ${BCMARKER_FILE} ${FEELMODELS_DIM}
#     ${FEELMODELS_ORDER_U} ${FEELMODELS_ORDERGEO} ${FEELMODELS_USE_CST_DENSITY_COEFFLAME} ${BCDESC_FILE} ${MESH_FILE} )

#   set(LIB_DEPENDS sm${APPLICATION_NAME} fsialg fsicore ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} ${ADDITIONAL_LIB_DEPEND} )
#   genExecutableBase(
#     ${APPLICATION_NAME} ${APPLICATION_DIR}
#     SRC ${MAIN_FILE}
#     LIB_DEPENDS ${LIB_DEPENDS}
#     CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in
#     )

# endmacro( genExecutableSolidMechanics )
# #############################################################################
# #############################################################################
# macro(fsi_add_solid_application)
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "DIM;U_ORDER;GEO_ORDER;GEO_DESC;BC_MARKERS;BC_DESC;SRC;DENSITY_COEFFLAME_TYPE;ADDITIONAL_LIB_DEPENDS"
#     ""
#     ${ARGN}
#     )

#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

#   if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_GEO_ORDER ) )
#     message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_GEO_ORDER")
#   endif()
#   set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
#   set(FEELMODELS_ORDER_U ${FEELMODELS_APP_U_ORDER})
#   set(FEELMODELS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

#   if ( FEELMODELS_APP_BC_MARKERS )
#     set(BCMARKER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_MARKERS})
#   else()
#     set(BCMARKER_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.bcmarkers )
#   endif()

#   if ( FEELMODELS_APP_GEO_DESC )
#     set(GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GEO_DESC})
#   else()
#     set(GEO_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#   endif()

#   set(DESC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_DESC})
#   set(MAIN_FILE ${FEELMODELS_APP_SRC})

#   if (FEELMODELS_APP_DENSITY_COEFFLAME_TYPE)
#     if ("${FEELMODELS_APP_DENSITY_COEFFLAME_TYPE}" STREQUAL "Cst" )
#       set(FEELMODELS_USE_CST_DENSITY_COEFFLAME 1)
#     elseif ("${FEELMODELS_APP_DENSITY_COEFFLAME_TYPE}" STREQUAL "P0d" )
#       set(FEELMODELS_USE_CST_DENSITY_COEFFLAME 0)
#     else()
#       message(FATAL_ERROR "DENSITY_COEFFLAME_TYPE ${FEELMODELS_APP_DENSITY_COEFFLAME_TYPE} : is not valid! It must be Cst or P0d")
#     endif()
#   else()
#     # default value
#     set(FEELMODELS_USE_CST_DENSITY_COEFFLAME 1)
#   endif()

#   if ( FEELMODELS_APP_ADDITIONAL_LIB_DEPENDS )
#     set( ADDITIONAL_LIB_DEPENDS ${FEELMODELS_APP_ADDITIONAL_LIB_DEPENDS} )
#   else()
#     # put feelpp ( duplication but allow to give an arg )
#     set( ADDITIONAL_LIB_DEPENDS ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} )
#   endif()

#   if (0)
#     MESSAGE("*** Arguments for solid application ${APPLICATION_NAME}")
#     MESSAGE("*** DIM ${FEELMODELS_DIM}")
#     MESSAGE("*** U_ORDER ${FEELMODELS_ORDER_U}")
#     MESSAGE("*** ORDERGEO ${FEELMODELS_ORDERGEO}")
#     MESSAGE("*** GEO_DESC ${GEO_FILE}")
#     MESSAGE("*** BC_MARKERS ${BCMARKER_FILE}")
#     MESSAGE("*** BC_DESC ${DESC_FILE}")
#     MESSAGE("*** SRC ${MAIN_FILE}")
#     MESSAGE("*** USE_CST_DENSITY_COEFFLAME ${FEELMODELS_USE_CST_DENSITY_COEFFLAME}")
#   endif()

#   genExecutableSolidMechanics(${APPLICATION_NAME} ${BCMARKER_FILE} ${FEELMODELS_DIM}
#     ${FEELMODELS_ORDER_U} ${FEELMODELS_ORDERGEO} ${FEELMODELS_USE_CST_DENSITY_COEFFLAME} ${DESC_FILE} ${GEO_FILE}
#     ${MAIN_FILE} ${ADDITIONAL_LIB_DEPENDS} )

# endmacro(fsi_add_solid_application)
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
# macro( genLibFluidMechanics APPLICATION_NAME APPLICATION_DIR BCMARKER_FILE FEELMODELS_DIM
#     FEELMODELS_ORDER_U FEELMODELS_ORDER_P FEELMODELS_ORDERGEO FEELMODELS_PRESSURE_IS_CONTINUOUS FEELMODELS_USE_CST_DENSITY_VISCOSITY FEELMODELS_DENSITY_VISCOSITY_ORDER FEELMODELS_USE_PERIODICITY BCDESC_FILE MESH_FILE )

#   set(FLUIDMECHANICS_DIM ${FEELMODELS_DIM})
#   set(FLUIDMECHANICS_ORDERGEO ${FEELMODELS_ORDERGEO})
#   set(FLUIDMECHANICS_ORDER_VELOCITY ${FEELMODELS_ORDER_U})
#   set(FLUIDMECHANICS_ORDER_PRESSURE ${FEELMODELS_ORDER_P})
#   if ("${FEELMODELS_PRESSURE_IS_CONTINUOUS}" STREQUAL "1" ) 
#     set(FLUIDMECHANICS_PRESSURE_IS_CONTINUOUS 1)
#     set(FEELMODELS_PRESSURE_CONTINUITY_TAG  c)
#   else()
#     set(FLUIDMECHANICS_PRESSURE_IS_CONTINUOUS 0)
#     set(FEELMODELS_PRESSURE_CONTINUITY_TAG  d)
#    endif()

#    set(FLUIDMECHANICS_ORDER_DENSITY_VISCOSITY ${FEELMODELS_DENSITY_VISCOSITY_ORDER})

#   if ("${FEELMODELS_USE_CST_DENSITY_VISCOSITY}" STREQUAL "1" ) 
#     set(FLUIDMECHANICS_USE_CST_DENSITY_VISCOSITY 1)
#     if ("${FEELMODELS_DENSITY_VISCOSITY_ORDER}" STREQUAL "0" ) 
#       set(FEELMODELS_USE_CST_DENSITY_VISCOSITY_TAG DVCst)
#     else()
#       set(FEELMODELS_USE_CST_DENSITY_VISCOSITY_TAG DVP${FEELMODELS_DENSITY_VISCOSITY_ORDER}c)
#     endif()
#   else()
#     set(FLUIDMECHANICS_USE_CST_DENSITY_VISCOSITY 0)
#     set(FEELMODELS_USE_CST_DENSITY_VISCOSITY_TAG DVP${FEELMODELS_DENSITY_VISCOSITY_ORDER}d)
#    endif()

#    if ( FEELMODELS_USE_PERIODICITY )
#      set(FLUIDMECHANICS_USE_PERIODICITY 1)
#      set(FLUIDMECHANICS_USE_PERIODICITY_TAG Periodic)
#    else()
#      set(FLUIDMECHANICS_USE_PERIODICITY 0)
#      unset(FLUIDMECHANICS_USE_PERIODICITY_TAG)
#    endif()

#   set(FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX ${FLUIDMECHANICS_DIM}dP${FLUIDMECHANICS_ORDER_VELOCITY}P${FLUIDMECHANICS_ORDER_PRESSURE}${FEELMODELS_PRESSURE_CONTINUITY_TAG}G${FLUIDMECHANICS_ORDERGEO}${FEELMODELS_USE_CST_DENSITY_VISCOSITY_TAG}${FLUIDMECHANICS_USE_PERIODICITY_TAG})
#   set(FEELMODELS_MODEL_SPECIFIC_NAME fluidmec${FEELMODELS_MODEL_SPECIFIC_NAME_SUFFIX})
#   set(LIBBASE_DIR ${PROJECT_BINARY_DIR}/fsi/fsimodel/${FEELMODELS_MODEL_SPECIFIC_NAME} )

#   set(LIBBASE_CHECK_PATH ${FEELMODELS_LIBBASE_CHECK_DIR}/${FEELMODELS_MODEL_SPECIFIC_NAME}.txt )
#   set(LIBBASE_NAME fsi_modelbase_${FEELMODELS_MODEL_SPECIFIC_NAME})

#   if ( NOT EXISTS ${LIBBASE_CHECK_PATH} )
#     #write empty file in orter to check if this lib has already define
#     file(WRITE ${LIBBASE_CHECK_PATH} "")

#     set(CODEGEN_FILES_TO_COPY
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbase.hpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbasecreate.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbaseothers.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbaseupdatelinear.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbaseupdatelinearweakbc.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbaseupdateresidual.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbaseupdateresidualstresstensorlaw.cpp
#       # ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecupdateptfixed.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbaseupdatejacobian.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbaseupdatejacobianstresstensorlaw.cpp
#       ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmecbaseupdatestabilisation.cpp )
#     set(CODEGEN_SOURCES
#       ${LIBBASE_DIR}/fluidmecbasecreate.cpp
#       ${LIBBASE_DIR}/fluidmecbaseothers.cpp
#       ${LIBBASE_DIR}/fluidmecbaseupdatelinear.cpp
#       ${LIBBASE_DIR}/fluidmecbaseupdatelinearweakbc.cpp
#       ${LIBBASE_DIR}/fluidmecbaseupdateresidual.cpp
#       ${LIBBASE_DIR}/fluidmecbaseupdateresidualstresstensorlaw.cpp
#       ${LIBBASE_DIR}/fluidmecbaseupdatejacobian.cpp
#       ${LIBBASE_DIR}/fluidmecbaseupdatejacobianstresstensorlaw.cpp
#       ${LIBBASE_DIR}/fluidmecbaseupdatestabilisation.cpp )
#     set(LIB_DEPENDS fsimesh fsialg fsicore fsifilters ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} )

#     # generate libmodelbase
#     genLibBase( ${APPLICATION_NAME} ${LIBBASE_DIR}
#       LIB_DEPENDS ${LIB_DEPENDS}
#       PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
#       PREFIX_LIBNAME ${LIBBASE_NAME}
#       FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
#       FILES_SOURCES ${CODEGEN_SOURCES}
#       CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fsifluidmecconfig.h.in
#       )

#   endif()


#   set(PREFIX_FILES_TO_COPY fluid)
#   set(PREFIX_LIB fm)
#   set(CODEGEN_FILES_TO_COPY
#     ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmec.hpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsimodel/fluidmec.cpp
#     )
#   set(CODEGEN_SOURCES
#     ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}/fluidmec.cpp
#     )
#   set(LIB_DEPENDS ${LIBBASE_NAME} fsimesh fsialg fsicore ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} )

#   set(FLUIDMECHANICS_INCLUDE_LIBBASE_HEADER 1)

#   genLibBase( ${APPLICATION_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}
#     MARKERS ${BCMARKER_FILE}
#     DESC ${BCDESC_FILE}
#     GEO ${MESH_FILE}
#     LIB_DEPENDS ${LIB_DEPENDS}
#     PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
#     PREFIX_LIBNAME ${PREFIX_LIB}${APPLICATION_NAME}
#     FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
#     FILES_SOURCES ${CODEGEN_SOURCES}
#     CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in
#     )

#   unset(FLUIDMECHANICS_INCLUDE_LIBBASE_HEADER)

# endmacro( genLibFluidMechanics )
# #############################################################################
# #############################################################################
# macro(genExecutableFluidMechanics APPLICATION_NAME BCMARKER_FILE FEELMODELS_DIM
#     FEELMODELS_ORDER_U FEELMODELS_ORDER_P FEELMODELS_ORDERGEO FEELMODELS_PRESSURE_IS_CONTINUOUS FEELMODELS_USE_CST_DENSITY_VISCOSITY FEELMODELS_DENSITY_VISCOSITY_ORDER FEELMODELS_USE_PERIODICITY BCDESC_FILE MESH_FILE MAIN_FILE )

#   resetToZeroAllPhysicalVariables()
#   SET(FLUIDMECHANICS 1 )

#   set(APPLICATION_DIR codeGen_${APPLICATION_NAME})
#   genLibFluidMechanics(${APPLICATION_NAME} ${APPLICATION_DIR}/fluid ${BCMARKER_FILE} ${FEELMODELS_DIM}
#     ${FEELMODELS_ORDER_U} ${FEELMODELS_ORDER_P} ${FEELMODELS_ORDERGEO} ${FEELMODELS_PRESSURE_IS_CONTINUOUS} ${FEELMODELS_USE_CST_DENSITY_VISCOSITY} ${FEELMODELS_DENSITY_VISCOSITY_ORDER} ${FEELMODELS_USE_PERIODICITY} ${BCDESC_FILE} ${MESH_FILE} )

#   set(LIB_DEPENDS fm${APPLICATION_NAME} fsimesh fsialg fsicore ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} )
#   genExecutableBase(
#     ${APPLICATION_NAME} ${APPLICATION_DIR}
#     SRC ${MAIN_FILE}
#     LIB_DEPENDS ${LIB_DEPENDS}
#     CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in
#     )

# endmacro( genExecutableFluidMechanics )
# #############################################################################
# #############################################################################
# macro(fsi_add_fluid_application)
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "DIM;U_ORDER;P_ORDER;P_CONTINUITY;GEO_ORDER;GEO_DESC;BC_MARKERS;BC_DESC;SRC;CFG;DENSITY_VISCOSITY_TYPE;DENSITY_VISCOSITY_ORDER"
#     "USE_PERIODICITY"
#     ${ARGN}
#     )

#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

#   if ( NOT ( FEELMODELS_APP_DIM OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_GEO_ORDER ) )
#     message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_U_ORDER OR FEELMODELS_APP_P_ORDER OR FEELMODELS_APP_GEO_ORDER")
#   endif()
#   set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
#   set(FEELMODELS_ORDER_U ${FEELMODELS_APP_U_ORDER})
#   set(FEELMODELS_ORDER_P ${FEELMODELS_APP_P_ORDER})
#   set(FEELMODELS_ORDERGEO ${FEELMODELS_APP_GEO_ORDER})

#   if ( FEELMODELS_APP_BC_MARKERS )
#     set(BCMARKER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_MARKERS})
#   else()
#     set(BCMARKER_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.bcmarkers )
#   endif()

#   if ( FEELMODELS_APP_GEO_DESC )
#     set(GEO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GEO_DESC})
#   else()
#     set(GEO_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#   endif()

#   if ( NOT FEELMODELS_APP_BC_DESC )
#     message(FATAL_ERROR "miss argument! FEELMODELS_APP_BC_DESC")
#   endif()
#   set(DESC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_DESC})
#   set(MAIN_FILE ${FEELMODELS_APP_SRC})

#   if ( FEELMODELS_APP_P_CONTINUITY )
#     if ("${FEELMODELS_APP_P_CONTINUITY}" STREQUAL "Continuous" )
#       set(FEELMODELS_PRESSURE_IS_CONTINUOUS 1)
#     elseif("${FEELMODELS_APP_P_CONTINUITY}" STREQUAL "Discontinuous" )
#       set(FEELMODELS_PRESSURE_IS_CONTINUOUS 0)
#     else()
#       message(FATAL_ERROR "P_CONTINUITY ${FEELMODELS_APP_P_CONTINUITY} : is not valid! It must be Continuous or Discontinuous")
#     endif()
#   else()
#     # default value
#     set(FEELMODELS_PRESSURE_IS_CONTINUOUS 1)
#   endif()

#   if (FEELMODELS_APP_DENSITY_VISCOSITY_TYPE)
#     if ("${FEELMODELS_APP_DENSITY_VISCOSITY_TYPE}" STREQUAL "Cst" )
#       set(FEELMODELS_USE_CST_DENSITY_VISCOSITY 1)
#     elseif ("${FEELMODELS_APP_DENSITY_VISCOSITY_TYPE}" STREQUAL "P0d" )
#       set(FEELMODELS_USE_CST_DENSITY_VISCOSITY 0)
#     else()
#       message(FATAL_ERROR "DENSITY_VISCOSITY_TYPE ${FEELMODELS_APP_DENSITY_VISCOSITY_TYPE} : is not valid! It must be Cst or P0d")
#     endif()
#   else()
#     # default value
#     set(FEELMODELS_USE_CST_DENSITY_VISCOSITY 1)
#   endif()

#   if (FEELMODELS_APP_DENSITY_VISCOSITY_ORDER)
#     set(FEELMODELS_DENSITY_VISCOSITY_ORDER ${FEELMODELS_APP_DENSITY_VISCOSITY_ORDER} )
#   else()
#     # default value
#     set(FEELMODELS_DENSITY_VISCOSITY_ORDER 0)
#   endif()

#   if ( FEELMODELS_APP_CFG )
#     foreach(  cfg ${FEELMODELS_APP_CFG} )
#       #      if ( EXISTS ${cfg} )
#       # extract cfg filename  to be copied in binary dir
#       get_filename_component( CFG_NAME ${cfg} NAME )
#       configure_file( ${cfg} ${CFG_NAME} )
#       #INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
#       #      else()
#       #        message(WARNING "Executable ${FEELMODELS_APP_NAME}: configuration file ${cfg} does not exist")
#       #      endif()
#     endforeach()
#   endif(FEELMODELS_APP_CFG)

#   if ( FEELMODELS_APP_USE_PERIODICITY )
#     set( FEELMODELS_USE_PERIODICITY 1 )
#   else()
#     set( FEELMODELS_USE_PERIODICITY 0 )
#   endif()

#   if (0)
#     MESSAGE("*** Arguments for fluid application ${APPLICATION_NAME}")
#     MESSAGE("*** DIM ${FEELMODELS_DIM}")
#     MESSAGE("*** ORDER_U ${FEELMODELS_ORDER_U}")
#     MESSAGE("*** ORDER_P ${FEELMODELS_ORDER_P}")
#     MESSAGE("*** ORDERGEO ${FEELMODELS_ORDERGEO}")
#     MESSAGE("*** PRESSURE_IS_CONTINUOUS ${FEELMODELS_PRESSURE_IS_CONTINUOUS}")
#     MESSAGE("*** USE_CST_DENSITY_VISCOSITY ${FEELMODELS_USE_CST_DENSITY_VISCOSITY}")
#     MESSAGE("*** MARKERS ${BCMARKER_FILE}")
#     MESSAGE("*** DESC ${DESC_FILE}")
#     MESSAGE("*** GEO ${GEO_FILE}")
#     MESSAGE("*** SRC ${MAIN_FILE}")
#   endif()

#   genExecutableFluidMechanics(${APPLICATION_NAME} ${BCMARKER_FILE} ${FEELMODELS_DIM}
#     ${FEELMODELS_ORDER_U} ${FEELMODELS_ORDER_P} ${FEELMODELS_ORDERGEO} ${FEELMODELS_PRESSURE_IS_CONTINUOUS} ${FEELMODELS_USE_CST_DENSITY_VISCOSITY} ${FEELMODELS_DENSITY_VISCOSITY_ORDER} ${FEELMODELS_USE_PERIODICITY}  ${DESC_FILE} ${GEO_FILE}
#     ${MAIN_FILE}
#      )

# endmacro(fsi_add_fluid_application)
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
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
#     if ("${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_TYPE}" STREQUAL "Cst" )
#       set(FLUID_USE_CST_DENSITY_VISCOSITY 1)
#     elseif ("${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_TYPE}" STREQUAL "P0d" )
#       set(FLUID_USE_CST_DENSITY_VISCOSITY 0)
#     else()
#       message(FATAL_ERROR "FLUID_DENSITY_VISCOSITY_TYPE ${FEELMODELS_APP_FLUID_DENSITY_VISCOSITY_TYPE} : is not valid! It must be Cst or P0d")
#     endif()
#   else()
#     # default value
#     set(FLUID_USE_CST_DENSITY_VISCOSITY 1)
#   endif()

#   if (FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE)
#     if ("${FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE}" STREQUAL "Cst" )
#       set(SOLID_USE_CST_DENSITY_COEFFLAME 1)
#     elseif ("${FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE}" STREQUAL "P0d" )
#       set(SOLID_USE_CST_DENSITY_COEFFLAME 0)
#     else()
#       message(FATAL_ERROR "SOLID_DENSITY_COEFFLAME_TYPE ${FEELMODELS_APP_SOLID_DENSITY_COEFFLAME_TYPE} : is not valid! It must be Cst or P0d")
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
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
# macro( genLibFBMFluidMechanics )
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "MARKERS;LIB_DEPENDS"
#     ""
#     ${ARGN}
#     )
#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})
#   CDR(APPLICATION_DIR  ${FEELMODELS_APP_DEFAULT_ARGS})
#   set(BCMARKER_FILE    ${FEELMODELS_APP_MARKERS})
#   set(LIB_DEPENDS      ${FEELMODELS_APP_LIB_DEPENDS})

#   set(PREFIX_FILES_TO_COPY fbm)
#   set(PREFIX_LIB fbm)
#   set(CODEGEN_FILES_TO_COPY
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbm.hpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbm.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmmonolitic.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmfixpoint.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmcouplingmonolitic.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmcouplingmonoliticlinear.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmcouplingmonoliticjacobian.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmcouplingmonoliticresidual.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmcouplingfixpoint.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmcouplingfixpointlinear.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmcouplingfixpointjacobian.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmcouplingfixpointresidual.cpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbminstantiation1.cpp )
#   set(CODEGEN_SOURCES
#       ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}/fbminstantiation1.cpp
#       ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}/fbmcouplingmonolitic.cpp
#       ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}/fbmcouplingfixpoint.cpp
#     )

#   genLibBase( ${APPLICATION_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}
#     MARKERS ${BCMARKER_FILE}
#     LIB_DEPENDS ${LIB_DEPENDS}
#     PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
#     PREFIX_LIBNAME ${PREFIX_LIB}${APPLICATION_NAME}
#     FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
#     FILES_SOURCES ${CODEGEN_SOURCES}
#     CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in
#     )
# endmacro( genLibFBMFluidMechanics )
# #############################################################################
# #############################################################################
# macro(genExecutableFBMFluidMechanics)
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "BC_MARKERS;DIM;GLOB_U_ORDER;GLOB_P_ORDER;GLOB_GEO_ORDER;GLOB_BC_DESC;GLOB_GEO_DESC;LOC_U_ORDER;LOC_P_ORDER;LOC_GEO_ORDER;LOC_BC_DESC;LOC_GEO_DESC;PERFO_U_ORDER;PERFO_P_ORDER;PERFO_GEO_ORDER;PERFO_BC_DESC;PERFO_GEO_DESC;SOLID_PERFO_U_ORDER;SOLID_PERFO_GEO_ORDER;SOLID_PERFO_BC_DESC;SOLID_PERFO_GEO_DESC;SRC"
#     ""
#     ${ARGN}
#     )

#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})
#   set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
#   set(FEELMODELS_ORDER_U_GLOB ${FEELMODELS_APP_GLOB_U_ORDER})
#   set(FEELMODELS_ORDER_P_GLOB ${FEELMODELS_APP_GLOB_P_ORDER})
#   set(FEELMODELS_ORDERGEO_GLOB ${FEELMODELS_APP_GLOB_GEO_ORDER})
#   set(FEELMODELS_ORDER_U_LOC ${FEELMODELS_APP_LOC_U_ORDER})
#   set(FEELMODELS_ORDER_P_LOC ${FEELMODELS_APP_LOC_P_ORDER})
#   set(FEELMODELS_ORDERGEO_LOC ${FEELMODELS_APP_LOC_GEO_ORDER})
#   #warning defalut!!
#   set(FLUID_PRESSURE_IS_CONTINUOUS 1)
#   set(FLUID_USE_CST_DENSITY_VISCOSITY 1)
#   set(SOLID_PERFO_USE_CST_DENSITY_COEFFLAME 1)

#   set(BCMARKER_FILE    ${FEELMODELS_APP_BC_MARKERS})
#   set(BCDESC_GLOB_FILE ${FEELMODELS_APP_GLOB_BC_DESC})
#   set(MESH_GLOB_FILE   ${FEELMODELS_APP_GLOB_GEO_DESC})
#   set(BCDESC_LOC_FILE  ${FEELMODELS_APP_LOC_BC_DESC})
#   set(MESH_LOC_FILE    ${FEELMODELS_APP_LOC_GEO_DESC})
#   set(MAIN_FILE        ${FEELMODELS_APP_SRC})

#   resetToZeroAllPhysicalVariables()
#   set(FBM_APPLICATION_DIR codeGen_${APPLICATION_NAME})
#   set(FBM_APPLICATION_NAME ${APPLICATION_NAME})

#   # generate the fluid0 library
#   SET(FLUIDMECHANICS0 1 )
#   SET(FLUIDMECHANICS1 0 )
#   set(FEELMODELS_DENSITY_VISCOSITY_ORDER 0)
#   set(FEELMODELS_USE_PERIODICITY 0)
#   genLibFluidMechanics(${FBM_APPLICATION_NAME}_NUM0 ${FBM_APPLICATION_DIR}/fluidmec0 ${BCMARKER_FILE} ${FEELMODELS_DIM}
#     ${FEELMODELS_ORDER_U_GLOB} ${FEELMODELS_ORDER_P_GLOB} ${FEELMODELS_ORDERGEO_GLOB} ${FLUID_PRESSURE_IS_CONTINUOUS}
#     ${FLUID_USE_CST_DENSITY_VISCOSITY} ${FEELMODELS_DENSITY_VISCOSITY_ORDER} ${FEELMODELS_USE_PERIODICITY} ${BCDESC_GLOB_FILE} ${MESH_GLOB_FILE})
#   # generate the fluid1 library
#   SET(FLUIDMECHANICS0 0 )
#   SET(FLUIDMECHANICS1 1 )
#   genLibFluidMechanics(${FBM_APPLICATION_NAME}_NUM1 ${FBM_APPLICATION_DIR}/fluidmec1 ${BCMARKER_FILE}  ${FEELMODELS_DIM}
#     ${FEELMODELS_ORDER_U_LOC} ${FEELMODELS_ORDER_P_LOC} ${FEELMODELS_ORDERGEO_LOC} ${FLUID_PRESSURE_IS_CONTINUOUS}
#     ${FLUID_USE_CST_DENSITY_VISCOSITY} ${FEELMODELS_DENSITY_VISCOSITY_ORDER} ${FEELMODELS_USE_PERIODICITY} ${BCDESC_LOC_FILE} ${MESH_LOC_FILE})
#   # generate the fluid2 library
#   if ( FEELMODELS_APP_PERFO_BC_DESC )
#     SET(FLUIDMECHANICS0 0 )
#     SET(FLUIDMECHANICS1 0 )
#     SET(FLUIDMECHANICS2 1 )
#     set(FEELMODELS_ORDER_U_PERFO ${FEELMODELS_APP_PERFO_U_ORDER})
#     set(FEELMODELS_ORDER_P_PERFO ${FEELMODELS_APP_PERFO_P_ORDER})
#     set(FEELMODELS_ORDERGEO_PERFO ${FEELMODELS_APP_PERFO_GEO_ORDER})
#     set(BCDESC_PERFO_FILE ${FEELMODELS_APP_PERFO_BC_DESC})
#     set(MESH_PERFO_FILE ${FEELMODELS_APP_PERFO_GEO_DESC})
#     genLibFluidMechanics(${FBM_APPLICATION_NAME}_NUM2 ${FBM_APPLICATION_DIR}/fluidmec2 ${BCMARKER_FILE}  ${FEELMODELS_DIM}
#       ${FEELMODELS_ORDER_U_PERFO} ${FEELMODELS_ORDER_P_PERFO} ${FEELMODELS_ORDERGEO_PERFO} ${FLUID_PRESSURE_IS_CONTINUOUS}
#       ${FLUID_USE_CST_DENSITY_VISCOSITY} ${FEELMODELS_DENSITY_VISCOSITY_ORDER} ${FEELMODELS_USE_PERIODICITY} ${BCDESC_PERFO_FILE} ${MESH_PERFO_FILE})
#   endif()
#   # generate the solid0 library
#   if ( FEELMODELS_APP_SOLID_PERFO_BC_DESC )
#     SET(FLUIDMECHANICS0 0 )
#     SET(FLUIDMECHANICS1 0 )
#     SET(FLUIDMECHANICS2 0 )
#     SET(SOLIDMECHANICS0 1 )
#     set(FEELMODELS_ORDER_U_SOLID_PERFO ${FEELMODELS_APP_SOLID_PERFO_U_ORDER})
#     set(FEELMODELS_ORDERGEO_SOLID_PERFO ${FEELMODELS_APP_SOLID_PERFO_GEO_ORDER})
#     set(BCDESC_SOLID_PERFO_FILE ${FEELMODELS_APP_SOLID_PERFO_BC_DESC})
#     set(MESH_SOLID_PERFO_FILE ${FEELMODELS_APP_SOLID_PERFO_GEO_DESC})
#     genLibSolidMechanics(${FBM_APPLICATION_NAME}_SOLIDNUM0 ${FBM_APPLICATION_DIR}/solidmec0 ${BCMARKER_FILE}  ${FEELMODELS_DIM}
#       ${FEELMODELS_ORDER_U_SOLID_PERFO} ${FEELMODELS_ORDERGEO_SOLID_PERFO} ${SOLID_PERFO_USE_CST_DENSITY_COEFFLAME}  ${BCDESC_SOLID_PERFO_FILE} ${MESH_SOLID_PERFO_FILE})

#   endif()

#   # generate fbm library
#   SET(FLUIDMECHANICS0 1 )
#   SET(FLUIDMECHANICS1 1 )
#   if ( FEELMODELS_APP_PERFO_BC_DESC )
#     SET(FLUIDMECHANICS2 1 )
#   endif()
#   if ( FEELMODELS_APP_SOLID_PERFO_BC_DESC )
#     SET(SOLIDMECHANICS0 1 )
#   endif()
#   SET(FBMTOOL 1 )

#   #CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/${FBM_APPLICATION_DIR}/fsicoreconfig.h )

#   set( FBM_LIB_DEPENDS fsimesh fsialg fsicore ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} fm${FBM_APPLICATION_NAME}_NUM0 fm${FBM_APPLICATION_NAME}_NUM1 )
#   if ( FEELMODELS_APP_PERFO_BC_DESC )
#     set( FBM_LIB_DEPENDS ${FBM_LIB_DEPENDS} fm${FBM_APPLICATION_NAME}_NUM2 )
#   endif()
#   if ( FEELMODELS_APP_SOLID_PERFO_BC_DESC )
#     set( FBM_LIB_DEPENDS ${FBM_LIB_DEPENDS} sm${FBM_APPLICATION_NAME}_SOLIDNUM0 )
#   endif()

#   genLibFBMFluidMechanics(${FBM_APPLICATION_NAME}_FBM ${FBM_APPLICATION_DIR}/FBM
#     MARKERS ${BCMARKER_FILE}
#     LIB_DEPENDS ${FBM_LIB_DEPENDS} )

#   #set(APPLICATION_DIR codeGen_${FBM_APPLICATION_NAME})

#   set( LIB_DEPENDS ${FBM_LIB_DEPENDS} fbm${FBM_APPLICATION_NAME}_FBM )
#   genExecutableBase(
#     ${FBM_APPLICATION_NAME} ${FBM_APPLICATION_DIR}
#     SRC ${MAIN_FILE}
#     LIB_DEPENDS ${LIB_DEPENDS}
#     CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in
#     )

# endmacro( genExecutableFBMFluidMechanics )
# #############################################################################
# #############################################################################
# macro(fsi_add_fbm_fluid_application)
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "BC_MARKERS;DIM;GLOB_U_ORDER;GLOB_P_ORDER;GLOB_GEO_ORDER;GLOB_BC_DESC;GLOB_GEO_DESC;LOC_U_ORDER;LOC_P_ORDER;LOC_GEO_ORDER;LOC_BC_DESC;LOC_GEO_DESC;PERFO_U_ORDER;PERFO_P_ORDER;PERFO_GEO_ORDER;PERFO_BC_DESC;PERFO_GEO_DESC;SOLID_PERFO_U_ORDER;SOLID_PERFO_GEO_ORDER;SOLID_PERFO_BC_DESC;SOLID_PERFO_GEO_DESC;SRC"
#     ""
#     ${ARGN}
#     )

#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

#   if (NOT (FEELMODELS_APP_DIM OR FEELMODELS_APP_GLOB_U_ORDER OR FEELMODELS_APP_GLOB_P_ORDER OR FEELMODELS_APP_GLOB_GEO_ORDER ) )
#     message(FATAL_ERROR "miss argument! FEELMODELS_APP_DIM OR FEELMODELS_APP_GLOB_U_ORDER OR FEELMODELS_APP_GLOB_P_ORDER OR FEELMODELS_APP_GLOB_GEO_ORDER")
#   endif()

#   if (NOT (FEELMODELS_APP_LOC_U_ORDER OR FEELMODELS_APP_LOC_P_ORDER OR FEELMODELS_APP_LOC_GEO_ORDER ) )
#     message(FATAL_ERROR "miss argument! FEELMODELS_APP_LOC_U_ORDER OR FEELMODELS_APP_LOC_P_ORDER OR FEELMODELS_APP_LOC_GEO_ORDER")
#   endif()

#   set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
#   set(FEELMODELS_ORDER_U_GLOB ${FEELMODELS_APP_GLOB_U_ORDER})
#   set(FEELMODELS_ORDER_P_GLOB ${FEELMODELS_APP_GLOB_P_ORDER})
#   set(FEELMODELS_ORDERGEO_GLOB ${FEELMODELS_APP_GLOB_GEO_ORDER})
#   set(FEELMODELS_ORDER_U_LOC ${FEELMODELS_APP_LOC_U_ORDER})
#   set(FEELMODELS_ORDER_P_LOC ${FEELMODELS_APP_LOC_P_ORDER})
#   set(FEELMODELS_ORDERGEO_LOC ${FEELMODELS_APP_LOC_GEO_ORDER})

#   if ( FEELMODELS_APP_BC_MARKERS )
#     set(BCMARKER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_MARKERS})
#   else()
#     set(BCMARKER_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.bcmarkers )
#   endif()

#   if ( FEELMODELS_APP_GLOB_GEO_DESC )
#     set(GEO_GLOB_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GLOB_GEO_DESC})
#   else()
#     set(GEO_GLOB_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#   endif()

#   if ( FEELMODELS_APP_LOC_GEO_DESC )
#     set(GEO_LOC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_LOC_GEO_DESC})
#   else()
#     set(GEO_LOC_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#   endif()

#   set(BCDESC_GLOB_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GLOB_BC_DESC})
#   set(BCDESC_LOC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_LOC_BC_DESC})
#   set(MAIN_FILE ${FEELMODELS_APP_SRC})

#   if (0)
#     MESSAGE("*** Arguments for fbm application ${APPLICATION_NAME}")
#     MESSAGE("*** MARKERS ${BCMARKER_FILE}")
#     MESSAGE("*** DESC_GLOB ${BCDESC_GLOB_FILE}")
#     MESSAGE("*** GEO_GLOB ${GEO_GLOB_FILE}")
#     MESSAGE("*** DESC_LOC ${BCDESC_LOC_FILE}")
#     MESSAGE("*** GEO_LOC ${GEO_LOC_FILE}")
#     MESSAGE("*** SRC ${MAIN_FILE}")
#   endif()

#   unset(BCDESC_PERFO_FILE )
#   unset(GEO_PERFO_FILE )
#   if ( FEELMODELS_APP_PERFO_BC_DESC )
#     if (NOT (FEELMODELS_APP_PERFO_U_ORDER OR FEELMODELS_APP_PERFO_P_ORDER OR FEELMODELS_APP_PERFO_GEO_ORDER ) )
#       message(FATAL_ERROR "miss argument! FEELMODELS_APP_PERFO_U_ORDER OR FEELMODELS_APP_PERFO_P_ORDER OR FEELMODELS_APP_PERFO_GEO_ORDER")
#     endif()

#     set(BCDESC_PERFO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_PERFO_BC_DESC})
#     if ( FEELMODELS_APP_PERFO_GEO_DESC )
#       set(GEO_PERFO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_PERFO_GEO_DESC})
#     else()
#       set(GEO_PERFO_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#     endif()
#   endif()
#   unset(BCDESC_SOLID_PERFO_FILE )
#   unset(GEO_SOLID_PERFO_FILE )
#   if ( FEELMODELS_APP_SOLID_PERFO_BC_DESC )
#     set(BCDESC_SOLID_PERFO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_SOLID_PERFO_BC_DESC})
#     if ( FEELMODELS_APP_SOLID_PERFO_GEO_DESC )
#       set(GEO_SOLID_PERFO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_SOLID_PERFO_GEO_DESC})
#     else()
#       set(GEO_SOLID_PERFO_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#     endif()
#   endif()

#   genExecutableFBMFluidMechanics(
#     ${APPLICATION_NAME}
#     BC_MARKERS ${BCMARKER_FILE}
#     DIM ${FEELMODELS_DIM}
#     GLOB_U_ORDER ${FEELMODELS_ORDER_U_GLOB}
#     GLOB_P_ORDER ${FEELMODELS_ORDER_P_GLOB}
#     GLOB_GEO_ORDER ${FEELMODELS_ORDERGEO_GLOB}
#     GLOB_BC_DESC ${BCDESC_GLOB_FILE}
#     GLOB_GEO_DESC ${GEO_GLOB_FILE}
#     LOC_U_ORDER ${FEELMODELS_ORDER_U_LOC}
#     LOC_P_ORDER ${FEELMODELS_ORDER_P_LOC}
#     LOC_GEO_ORDER ${FEELMODELS_ORDERGEO_LOC}
#     LOC_BC_DESC ${BCDESC_LOC_FILE}
#     LOC_GEO_DESC ${GEO_LOC_FILE}
#     PERFO_U_ORDER ${FEELMODELS_APP_PERFO_U_ORDER}
#     PERFO_P_ORDER ${FEELMODELS_APP_PERFO_P_ORDER}
#     PERFO_GEO_ORDER ${FEELMODELS_APP_PERFO_GEO_ORDER}
#     PERFO_BC_DESC ${BCDESC_PERFO_FILE}
#     PERFO_GEO_DESC ${GEO_PERFO_FILE}
#     SOLID_PERFO_U_ORDER ${FEELMODELS_APP_SOLID_PERFO_U_ORDER}
#     SOLID_PERFO_GEO_ORDER ${FEELMODELS_APP_SOLID_PERFO_GEO_ORDER}
#     SOLID_PERFO_BC_DESC ${BCDESC_SOLID_PERFO_FILE}
#     SOLID_PERFO_GEO_DESC ${GEO_SOLID_PERFO_FILE}
#     SRC ${MAIN_FILE}
#     )

# endmacro(fsi_add_fbm_fluid_application)
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
# macro( genLibFBMThermoDynamics )
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "MARKERS;LIB_DEPENDS"
#     ""
#     ${ARGN}
#     )
#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})
#   CDR(APPLICATION_DIR  ${FEELMODELS_APP_DEFAULT_ARGS})
#   set(BCMARKER_FILE    ${FEELMODELS_APP_MARKERS})
#   set(LIB_DEPENDS      ${FEELMODELS_APP_LIB_DEPENDS})

#   set(PREFIX_FILES_TO_COPY fbmthermodyn)
#   set(PREFIX_LIB fbmthermodyn)
#   set(CODEGEN_FILES_TO_COPY
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmthermodyn.hpp
#     ${PROJECT_SOURCE_DIR}/fsi/fsicore/fbmthermodyn.cpp)
#   set(CODEGEN_SOURCES
#       ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}/fbmthermodyn.cpp
#     )

#   genLibBase( ${APPLICATION_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${APPLICATION_DIR}
#     MARKERS ${BCMARKER_FILE}
#     LIB_DEPENDS ${LIB_DEPENDS}
#     PREFIX_INCLUDE_USERCONFIG ${PREFIX_FILES_TO_COPY}
#     PREFIX_LIBNAME ${PREFIX_LIB}${APPLICATION_NAME}
#     FILES_TO_COPY ${CODEGEN_FILES_TO_COPY}
#     FILES_SOURCES ${CODEGEN_SOURCES}
#     CONFIG_PATH ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in
#     )
# endmacro( genLibFBMThermoDynamics )
# #############################################################################
# #############################################################################
# macro(genExecutableFBMThermoDynamics)
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "BC_MARKERS;DIM;GLOB_T_ORDER;GLOB_GEO_ORDER;GLOB_BC_DESC;GLOB_GEO_DESC;LOC_T_ORDER;LOC_GEO_ORDER;LOC_BC_DESC;LOC_GEO_DESC;PERFO_T_ORDER;PERFO_GEO_ORDER;PERFO_BC_DESC;PERFO_GEO_DESC;SRC"
#     ""
#     ${ARGN}
#     )

#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

#   set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
#   set(FEELMODELS_ORDERPOLY_GLOB ${FEELMODELS_APP_GLOB_T_ORDER})
#   set(FEELMODELS_ORDERGEO_GLOB ${FEELMODELS_APP_GLOB_GEO_ORDER})
#   set(FEELMODELS_ORDERPOLY_LOC ${FEELMODELS_APP_LOC_T_ORDER})
#   set(FEELMODELS_ORDERGEO_LOC ${FEELMODELS_APP_LOC_GEO_ORDER})

#   set(BCMARKER_FILE    ${FEELMODELS_APP_BC_MARKERS})
#   set(BCDESC_GLOB_FILE ${FEELMODELS_APP_GLOB_BC_DESC})
#   set(MESH_GLOB_FILE   ${FEELMODELS_APP_GLOB_GEO_DESC})
#   set(BCDESC_LOC_FILE  ${FEELMODELS_APP_LOC_BC_DESC})
#   set(MESH_LOC_FILE    ${FEELMODELS_APP_LOC_GEO_DESC})
#   set(MAIN_FILE        ${FEELMODELS_APP_SRC})

#   resetToZeroAllPhysicalVariables()
#   set(FBM_APPLICATION_DIR codeGen_${APPLICATION_NAME})
#   set(FBM_APPLICATION_NAME ${APPLICATION_NAME})

#   # generate the thermo0 library
#   SET(THERMODYNAMICS0 1 )
#   SET(THERMODYNAMICS1 0 )
#   SET(THERMODYNAMICS2 0 )
#   genLibThermoDynamics(${FBM_APPLICATION_NAME}_NUM0 ${FBM_APPLICATION_DIR}/thermodyn0
#     ${FEELMODELS_DIM} ${FEELMODELS_ORDERPOLY_GLOB} ${FEELMODELS_ORDERGEO_GLOB}
#     ${BCMARKER_FILE} ${BCDESC_GLOB_FILE} ${MESH_GLOB_FILE} )
#   # generate the thermo1 library
#   SET(THERMODYNAMICS0 0 )
#   SET(THERMODYNAMICS1 1 )
#   SET(THERMODYNAMICS2 0 )
#   genLibThermoDynamics(${FBM_APPLICATION_NAME}_NUM1 ${FBM_APPLICATION_DIR}/thermodyn1
#     ${FEELMODELS_DIM} ${FEELMODELS_ORDERPOLY_LOC} ${FEELMODELS_ORDERGEO_LOC}
#     ${BCMARKER_FILE} ${BCDESC_LOC_FILE} ${MESH_LOC_FILE} )
#   # generate the thermo2 library
#   if ( FEELMODELS_APP_PERFO_BC_DESC )
#     SET(THERMODYNAMICS0 0 )
#     SET(THERMODYNAMICS1 0 )
#     SET(THERMODYNAMICS2 1 )
#     set(FEELMODELS_ORDERPOLY_PERFO ${FEELMODELS_APP_PERFO_T_ORDER})
#     set(FEELMODELS_ORDERGEO_PERFO ${FEELMODELS_APP_PERFO_GEO_ORDER})
#     set(BCDESC_PERFO_FILE ${FEELMODELS_APP_PERFO_BC_DESC})
#     set(MESH_PERFO_FILE ${FEELMODELS_APP_PERFO_GEO_DESC})
#     genLibThermoDynamics(${FBM_APPLICATION_NAME}_NUM2 ${FBM_APPLICATION_DIR}/thermodyn2
#       ${FEELMODELS_DIM} ${FEELMODELS_ORDERPOLY_PERFO} ${FEELMODELS_ORDERGEO_PERFO}
#       ${BCMARKER_FILE} ${BCDESC_PERFO_FILE} ${MESH_PERFO_FILE})
#   endif()

#   # generate fbm library
#   SET(THERMODYNAMICS0 1 )
#   SET(THERMODYNAMICS1 1 )
#   if ( FEELMODELS_APP_PERFO_BC_DESC )
#     SET(THERMODYNAMICS2 1 )
#   endif()
#   SET(FBMTOOL 1 )

#   CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/fsi/fsicore/fsicoreconfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/${FBM_APPLICATION_DIR}/fsicoreconfig.h )

#   set( FBM_LIB_DEPENDS fsialg fsicore ${FEELMODELS_LIBRARIES_TO_LINK_FROM_FEELPP} thermo${FBM_APPLICATION_NAME}_NUM0 thermo${FBM_APPLICATION_NAME}_NUM1 )
#   if ( FEELMODELS_APP_PERFO_BC_DESC )
#     set( FBM_LIB_DEPENDS ${FBM_LIB_DEPENDS} thermo${FBM_APPLICATION_NAME}_NUM2 )
#   endif()

#   genLibFBMThermoDynamics(${FBM_APPLICATION_NAME}_FBM ${FBM_APPLICATION_DIR}/FBM
#     MARKERS ${BCMARKER_FILE}
#     LIB_DEPENDS ${FBM_LIB_DEPENDS} )

#   set(APPLICATION_DIR codeGen_${FBM_APPLICATION_NAME})

#   set( LIB_DEPENDS ${FBM_LIB_DEPENDS} fbmthermodyn${FBM_APPLICATION_NAME}_FBM )
#   genExecutableBase( ${FBM_APPLICATION_NAME} ${APPLICATION_DIR} SRC ${MAIN_FILE} LIB_DEPENDS ${LIB_DEPENDS} )

# endmacro( genExecutableFBMThermoDynamics )
# #############################################################################
# #############################################################################
# macro(fsi_add_fbm_thermo_application)
#   PARSE_ARGUMENTS(FEELMODELS_APP
#     "BC_MARKERS;DIM;GLOB_T_ORDER;GLOB_GEO_ORDER;GLOB_BC_DESC;GLOB_GEO_DESC;LOC_T_ORDER;LOC_GEO_ORDER;LOC_BC_DESC;LOC_GEO_DESC;PERFO_T_ORDER;PERFO_GEO_ORDER;PERFO_BC_DESC;PERFO_GEO_DESC;SRC"
#     ""
#     ${ARGN}
#     )

#   CAR(APPLICATION_NAME ${FEELMODELS_APP_DEFAULT_ARGS})

#   set(FEELMODELS_DIM ${FEELMODELS_APP_DIM})
#   set(FEELMODELS_ORDERPOLY_GLOB ${FEELMODELS_APP_GLOB_T_ORDER})
#   set(FEELMODELS_ORDERGEO_GLOB ${FEELMODELS_APP_GLOB_GEO_ORDER})
#   set(FEELMODELS_ORDERPOLY_LOC ${FEELMODELS_APP_LOC_T_ORDER})
#   set(FEELMODELS_ORDERGEO_LOC ${FEELMODELS_APP_LOC_GEO_ORDER})
#   set(FEELMODELS_ORDERPOLY_PERFO ${FEELMODELS_APP_PERFO_T_ORDER})
#   set(FEELMODELS_ORDERGEO_PERFO ${FEELMODELS_APP_PERFO_GEO_ORDER})


#   if ( FEELMODELS_APP_BC_MARKERS )
#     set(BCMARKER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_BC_MARKERS})
#   else()
#     set(BCMARKER_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.bcmarkers )
#   endif()

#   if ( FEELMODELS_APP_GLOB_GEO_DESC )
#     set(GEO_GLOB_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GLOB_GEO_DESC})
#   else()
#     set(GEO_GLOB_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#   endif()

#   if ( FEELMODELS_APP_LOC_GEO_DESC )
#     set(GEO_LOC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_LOC_GEO_DESC})
#   else()
#     set(GEO_LOC_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#   endif()

#   set(BCDESC_GLOB_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_GLOB_BC_DESC})
#   set(BCDESC_LOC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_LOC_BC_DESC})
#   set(MAIN_FILE ${FEELMODELS_APP_SRC})

#   if (0)
#     MESSAGE("*** Arguments for fbm application ${APPLICATION_NAME}")
#     MESSAGE("*** MARKERS ${BCMARKER_FILE}")
#     MESSAGE("*** ORDERPOLY_GLOB ${FEELMODELS_ORDERPOLY_GLOB}")
#     MESSAGE("*** DESC_GLOB ${BCDESC_GLOB_FILE}")
#     MESSAGE("*** GEO_GLOB ${GEO_GLOB_FILE}")
#     MESSAGE("*** DESC_LOC ${BCDESC_LOC_FILE}")
#     MESSAGE("*** GEO_LOC ${GEO_LOC_FILE}")
#     MESSAGE("*** SRC ${MAIN_FILE}")
#   endif()

#   unset(BCDESC_PERFO_FILE )
#   unset(GEO_PERFO_FILE )
#   if ( FEELMODELS_APP_PERFO_BC_DESC )
#     set(BCDESC_PERFO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_PERFO_BC_DESC})
#     if ( FEELMODELS_APP_PERFO_GEO_DESC )
#       set(GEO_PERFO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FEELMODELS_APP_PERFO_GEO_DESC})
#     else()
#       set(GEO_PERFO_FILE ${PROJECT_SOURCE_DIR}/fsi/fsicore/default.mesh )
#     endif()
#   endif()

#   genExecutableFBMThermoDynamics(
#     ${APPLICATION_NAME}
#     BC_MARKERS      ${BCMARKER_FILE}
#     DIM             ${FEELMODELS_DIM}
#     GLOB_T_ORDER    ${FEELMODELS_ORDERPOLY_GLOB}
#     GLOB_GEO_ORDER  ${FEELMODELS_ORDERGEO_GLOB}
#     GLOB_GEO_DESC   ${GEO_GLOB_FILE}
#     GLOB_BC_DESC    ${BCDESC_GLOB_FILE}
#     LOC_T_ORDER     ${FEELMODELS_ORDERPOLY_LOC}
#     LOC_GEO_ORDER   ${FEELMODELS_ORDERGEO_LOC}
#     LOC_GEO_DESC    ${GEO_LOC_FILE}
#     LOC_BC_DESC     ${BCDESC_LOC_FILE}
#     PERFO_T_ORDER   ${FEELMODELS_ORDERPOLY_PERFO}
#     PERFO_GEO_ORDER ${FEELMODELS_ORDERGEO_PERFO}
#     PERFO_GEO_DESC  ${GEO_PERFO_FILE}
#     PERFO_BC_DESC   ${BCDESC_PERFO_FILE}
#     SRC             ${MAIN_FILE}
#     )

# endmacro(fsi_add_fbm_thermo_application)
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################





