
set(FEELPP_MOR_LIBRARIES_TEXT)
set(FEELPP_MOR_ALL_MODELS_CMAKE_VARIABLE)
foreach(mylib ${FEELPP_MOR_LIBRARIES})
  set(_mormodelVarName ${mylib})
  set(FEELPP_MOR_ALL_MODELS_CMAKE_VARIABLE "${FEELPP_MOR_ALL_MODELS_CMAKE_VARIABLE} \${${_mormodelVarName}}")
  set(FEELPP_MOR_MODEL_LIBRARIES_WITH_SPACE)
  foreach( _themorlib ${${_mormodelVarName}})
      if ( TARGET ${_themorlib} )
        list(APPEND FEELPP_MOR_MODEL_LIBRARIES_WITH_SPACE  " $<TARGET_FILE:${_themorlib}> ")
      else()
        list(APPEND FEELPP_MOR_MODEL_LIBRARIES_WITH_SPACE " ${_themorlib} " )
        endif()
  endforeach()
  set(FEELPP_MOR_MODEL_LIBRARIES_TEXT "set(${_mormodelVarName} ${FEELPP_MOR_MODEL_LIBRARIES_WITH_SPACE})")
  string(REPLACE ";" "" FEELPP_MOR_MODEL_LIBRARIES_TEXT ${FEELPP_MOR_MODEL_LIBRARIES_TEXT} )
  set(FEELPP_MOR_LIBRARIES_TEXT  " ${FEELPP_MOR_LIBRARIES_TEXT}\n${FEELPP_MOR_MODEL_LIBRARIES_TEXT}")
endforeach()

set(FEELPP_MOR_LIBRARIES_TEXT "${FEELPP_MOR_LIBRARIES_TEXT}\nset(FEELPP_MOR_LIBRARIES ${FEELPP_MOR_ALL_MODELS_CMAKE_VARIABLE})")
#message("FEELPP_MOR_LIBRARIES_TEXT=${FEELPP_MOR_LIBRARIES_TEXT}")
file( GENERATE OUTPUT ${FEELPP_MOR_BINARY_DIR}/cmake/feelpp.mor.libraries.config.cmake
  CONTENT ${FEELPP_MOR_LIBRARIES_TEXT} )

set(FEELPP_MOR_INCLUDE_DIR_TEXT "set(FEELPP_MOR_INCLUDE_DIR ${FEELPP_MOR_SOURCE_DIR} )")
file( WRITE ${FEELPP_MOR_BINARY_DIR}/cmake/feelpp.mor.includes.config.cmake ${FEELPP_MOR_INCLUDE_DIR_TEXT} )


set( FEELPP_MOR_CONFIG_LIB_FILE  ${FEELPP_MOR_BINARY_DIR}/cmake/Feel++-MORConfig.cmake)
configure_file(${FEELPP_MOR_SOURCE_DIR}/cmake/Feel++-MORConfig.cmake.in ${FEELPP_MOR_CONFIG_LIB_FILE} @ONLY)
install(FILES  ${FEELPP_MOR_CONFIG_LIB_FILE} DESTINATION share/feelpp/mor/cmake/ COMPONENT Devel)

configure_file(${FEELPP_MOR_SOURCE_DIR}/cmake/feelpp.mor.install.config.cmake.in ${FEELPP_MOR_BINARY_DIR}/cmake/feelpp.mor.install.config.cmake @ONLY)
install(SCRIPT ${FEELPP_MOR_BINARY_DIR}/cmake/feelpp.mor.install.config.cmake COMPONENT Devel)
