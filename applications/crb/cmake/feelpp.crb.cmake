
set(FEELPP_CRB_LIBRARIES_TEXT)
set(FEELPP_CRB_ALL_MODELS_CMAKE_VARIABLE)
foreach(mylib ${FEELPP_CRB_LIBRARIES})
  set(_crbmodelVarName ${mylib})
  set(FEELPP_CRB_ALL_MODELS_CMAKE_VARIABLE "${FEELPP_CRB_ALL_MODELS_CMAKE_VARIABLE} \${${_crbmodelVarName}}")
  set(FEELPP_CRB_MODEL_LIBRARIES_WITH_SPACE)
  foreach( _thecrblib ${${_crbmodelVarName}})
      if ( TARGET ${_thecrblib} )
        list(APPEND FEELPP_CRB_MODEL_LIBRARIES_WITH_SPACE  " $<TARGET_FILE:${_thecrblib}> ")
      else()
        list(APPEND FEELPP_CRB_MODEL_LIBRARIES_WITH_SPACE " ${_thecrblib} " )
        endif()
  endforeach()
  set(FEELPP_CRB_MODEL_LIBRARIES_TEXT "set(${_crbmodelVarName} ${FEELPP_CRB_MODEL_LIBRARIES_WITH_SPACE})")
  string(REPLACE ";" "" FEELPP_CRB_MODEL_LIBRARIES_TEXT ${FEELPP_CRB_MODEL_LIBRARIES_TEXT} )
  set(FEELPP_CRB_LIBRARIES_TEXT  " ${FEELPP_CRB_LIBRARIES_TEXT}\n${FEELPP_CRB_MODEL_LIBRARIES_TEXT}")
endforeach()

set(FEELPP_CRB_LIBRARIES_TEXT "${FEELPP_CRB_LIBRARIES_TEXT}\nset(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_ALL_MODELS_CMAKE_VARIABLE})")
#message("FEELPP_CRB_LIBRARIES_TEXT=${FEELPP_CRB_LIBRARIES_TEXT}")
file( GENERATE OUTPUT ${FEELPP_CRB_BINARY_DIR}/cmake/feelpp.crb.libraries.config.cmake
  CONTENT ${FEELPP_CRB_LIBRARIES_TEXT} )

set(FEELPP_CRB_INCLUDE_DIR_TEXT "set(FEELPP_CRB_INCLUDE_DIR ${FEELPP_CRB_SOURCE_DIR} )")
file( WRITE ${FEELPP_CRB_BINARY_DIR}/cmake/feelpp.crb.includes.config.cmake ${FEELPP_CRB_INCLUDE_DIR_TEXT} )


set( FEELPP_CRB_CONFIG_LIB_FILE  ${FEELPP_CRB_BINARY_DIR}/cmake/Feel++-CRBConfig.cmake)
configure_file(${FEELPP_CRB_SOURCE_DIR}/cmake/Feel++-CRBConfig.cmake.in ${FEELPP_CRB_CONFIG_LIB_FILE} @ONLY)
install(FILES  ${FEELPP_CRB_CONFIG_LIB_FILE} DESTINATION share/feelpp/crb/cmake/ COMPONENT Devel)

configure_file(${FEELPP_CRB_SOURCE_DIR}/cmake/feelpp.crb.install.config.cmake.in ${FEELPP_CRB_BINARY_DIR}/cmake/feelpp.crb.install.config.cmake @ONLY)
install(SCRIPT ${FEELPP_CRB_BINARY_DIR}/cmake/feelpp.crb.install.config.cmake COMPONENT Devel)
