# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
#       Date: 2013-10-21
#
#  Copyright (C) 2013 Université de Strasbourg
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 3.0 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#

# Set harts to not found by default
SET(HARTS_FOUND FALSE)

# Try to find harts installation in the system
FIND_PATH(HARTS_INCLUDE_DIR RunTimeSystem/Model/RunTimeSysEnv.h PATH_SUFFIXES HARTS)

FOREACH(_lib HARTSRuntimeSys HARTSUtils)
    FIND_LIBRARY(LIB_SUB_${_lib} ${_lib})

    if(LIB_SUB_${_lib})
        set(HARTS_LIBRARY ${HARTS_LIBRARY} ${LIB_SUB_${_lib}})
    else()
        #message(WARNING "A component of the HARTS library was not found in the system (${_lib})")
        set(HARTS_LIBRARY "")
        unset(HARTS_LIBRARY CACHE)
        break()
    endif()
endforeach()

IF (HARTS_INCLUDE_DIR AND HARTS_LIBRARY)
    SET(HARTS_FOUND TRUE)
ELSE()
    # Try to find a harts git checkout in contrib
    FIND_PATH(HARTS_SOURCE_DIR CMakeLists.txt
        PATH_SUFFIXES harts HARTS
        HINTS ${CMAKE_SOURCE_DIR}/contrib)

    # if we found harts in the contrib, we compile it
    IF(HARTS_SOURCE_DIR)
        message(STATUS "Building harts in ${CMAKE_BINARY_DIR}/contrib/harts-compile...")
        execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/contrib/harts-compile)
        execute_process(
            COMMAND cmake ${HARTS_SOURCE_DIR} -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/contrib/harts
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/harts-compile
            #      OUTPUT_QUIET
            )
        execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/contrib/harts)
        execute_process(
            COMMAND make -k -j${NProcs2} install
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/harts-compile
            #      OUTPUT_QUIET
            )

        FIND_PATH(HARTS_INCLUDE_DIR RunTimeSystem/Model/RunTimeSysEnv.h 
            HINTS ${CMAKE_BINARY_DIR}/contrib/harts/include/HARTS)

        FOREACH(_lib HARTSRuntimeSys HARTSUtils)
            FIND_LIBRARY(LIB_SUB_${_lib} ${_lib}
                  HINTS ${CMAKE_BINARY_DIR}/contrib/harts/lib)

            if(LIB_SUB_${_lib})
                set(HARTS_LIBRARY ${HARTS_LIBRARY} ${LIB_SUB_${_lib}})
            else()
                #message(WARNING "A component of the HARTS library was not found (${_lib})")
                set(HARTS_LIBRARY "")
                break()
            endif()
        endforeach()

        IF (HARTS_INCLUDE_DIR AND HARTS_LIBRARY)
            SET(HARTS_FOUND TRUE)
        ENDIF()
    ENDIF()
ENDIF()

# If we found harts, we create the variables holding the libraries and includes
IF (HARTS_FOUND)
    IF (NOT HARTS_FIND_QUIETLY)
        MESSAGE(STATUS "Found HARTS: ${HARTS_LIBRARY}")
    ENDIF (NOT HARTS_FIND_QUIETLY)
    set(FEELPP_HAS_HARTS 1)
    set(HARTS_INCLUDES ${HARTS_INCLUDE_DIR} CACHE STRING "HARTS include path")
    set(HARTS_LIBRARIES ${HARTS_LIBRARY} CACHE STRING "HARTS libraries")
ELSE ()
    MESSAGE(FATAL_ERROR "Could not find HARTS")
ENDIF ()

# Variables clean up
FOREACH(_lib HARTSRuntimeSys HARTSUtils)
    if(LIB_SUB_${_lib})
        unset(LIB_SUB_${_lib} CACHE)
    endif()
endforeach()
if(HARTS_LIBRARY)
    unset(HARTS_LIBRARY CACHE)
endif()
if(HARTS_INCLUDE_DIR)
    unset(HARTS_INCLUDE_DIR CACHE)
endif()

MARK_AS_ADVANCED( HARTS_INCLUDES HARTS_LIBRARIES )
