
message("Hola feelpp.generate.cling.interpreter.bash.cmake")

if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/cmake/modules/Feel++Config.cmake )
  message(ERROR "config file not found  ${CMAKE_CURRENT_BINARY_DIR}/cmake/modules/Feel++Config.cmake")
  return()
endif()


set(FEELPP_DONT_SETUP_CMAKE 1)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/modules/Feel++Config.cmake)
unset(FEELPP_DONT_SETUP_CMAKE)

if (0)
message("FEELPP_LIBRARY=${FEELPP_LIBRARY}")
message("FEELPP_LIBRARIES=${FEELPP_LIBRARIES}")
message("FEELPP_INCLUDE_DIR=${FEELPP_INCLUDE_DIR}")
message("FEELPP_DEPS_INCLUDE_DIR=${FEELPP_DEPS_INCLUDE_DIR}")
endif()

file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/feel++
"#!/usr/bin/env bash
printf \"\\\n\\e[1m\\e[1;31m\
_____ _____ _____ __      _     _
|   __|   __|   __|  |   _| |_ _| |
|   __|   __|   __|  |__|_   _|_   _|
|__|  |_____|_____|_____| |_|   |_|
\\e[0m \\e[1m\\e[1;36m
Feel++ interpreter (cling based)
For more infos, see www.feelpp.org
Type '.help' for help, '.q' to exit
\\e[0m
\"
"
)

file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/feel++
"${FEELPP_CLING_BIN} \\"
)


