
find_package(Cling)
if(NOT Cling_FOUND)
  return()
endif()

set(FEELPP_CLING_BIN ${CLING_BIN})
set(FEELPP_HAS_CLING_INTERPRETER 1)
set(FEELPP_ENABLED_MODULES "${FEELPP_ENABLED_MODULES} Cling/Interpreter" )

add_custom_target(feelpp_cling_interpreter
  DEPENDS feelpp
  COMMAND ${CMAKE_COMMAND} -P "${CMAKE_SOURCE_DIR}/cmake/modules/feelpp.generate.cling.interpreter.bash.cmake"
  )