# This file generate a bash script for the feel++ interpreter.
# It configure feel++ libraries for cling, but not a builtin cling executable!

# Create the script.
file( WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"#!/bin/sh
printf \"\
 _____ _____ _____ __      _     _
|   __|   __|   __|  |   _| |_ _| |
|   __|   __|   __|  |__|_   _|_   _|
|__|  |_____|_____|_____| |_|   |_|

Feel++ interpreter (cling based)
For more infos, see www.feelpp.org
Type '.help' for help, '.q' or 'exit' to quit

\"
${CLING_BIN} --nologo $1"
    )
# Set permissions as executable script.
file( INSTALL ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
    DESTINATION ${CLING_INSTALL_PREFIX}
    )
