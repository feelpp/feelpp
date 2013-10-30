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

FIND_PATH(HARTS_INCLUDE_DIR RunTimeSystem/Model/RunTimeSysEnv.h PATH_SUFFIXES HARTS)

FOREACH(_lib HARTSRuntimeSys HARTSUtils)
    FIND_LIBRARY(LIB_SUB_${_lib} ${_lib})

    if(LIB_SUB_${_lib})
        set(HARTS_LIBRARY ${HARTS_LIBRARY} ${LIB_SUB_${_lib}})
    else()
        message(WARNING "A component of the HARTS library was not found (${_lib})")
        set(HARTS_LIBRARY "")
        break()
    endif()
endforeach()

MESSAGE(STATUS "HARTS: ${HARTS_LIBRARY}")

IF (HARTS_INCLUDE_DIR AND HARTS_LIBRARY)
    SET(HARTS_FOUND TRUE)
ENDIF (HARTS_INCLUDE_DIR AND HARTS_LIBRARY)

IF (HARTS_FOUND)
    IF (NOT HARTS_FIND_QUIETLY)
        MESSAGE(STATUS "Found HARTS: ${HARTS_LIBRARY}")
    ENDIF (NOT HARTS_FIND_QUIETLY)
    set(FEELPP_HAS_HARTS 1)
    set(HARTS_INCLUDES ${HARTS_INCLUDE_DIR} CACHE STRING "HARTS include path")
    set(HARTS_LIBRARIES ${HARTS_LIBRARY})
ELSE (HARTS_FOUND)
    IF (HARTS_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find HARTS")
    ENDIF (HARTS_FIND_REQUIRED)
ENDIF (HARTS_FOUND)

MARK_AS_ADVANCED( HARTS_INCLUDES HARTS_LIBRARIES )
