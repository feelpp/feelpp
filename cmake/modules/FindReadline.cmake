# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2010-04-10
#
#  Copyright (C) 2010 Université Joseph Fourier
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
FIND_PATH(READLINE_INCLUDE_DIR
  readline.h
  PATHS /opt/local/include
  PATH_SUFFIXES readline
  NO_DEFAULT_PATH
  )
if ( NOT READLINE_INCLUDE_DIR )
  FIND_PATH(READLINE_INCLUDE_DIR
    readline.h
    PATH_SUFFIXES readline
    )
endif()
FIND_LIBRARY(READLINE_LIBRARY NAMES readline PATHS /opt/local/lib/ NO_DEFAULT_PATH)
if ( NOT READLINE_LIBRARY )
  FIND_LIBRARY(READLINE_LIBRARY NAMES readline PATH_SUFFIXES x86_64-linux-gnu)
endif()

MESSAGE(STATUS "readline: ${READLINE_LIBRARY}")

IF (READLINE_INCLUDE_DIR AND READLINE_LIBRARY)
  SET(READLINE_FOUND TRUE)
  set(READLINE_INCLUDE_DIRS ${READLINE_INCLUDE_DIR})
  set(READLINE_LIBRARIES ${READLINE_LIBRARY} )
ENDIF (READLINE_INCLUDE_DIR AND READLINE_LIBRARY)

IF (READLINE_FOUND)
   IF (NOT Readline_FIND_QUIETLY)
      MESSAGE(STATUS "Found GNU readline: ${READLINE_LIBRARY}")
   ENDIF (NOT Readline_FIND_QUIETLY)
ELSE (READLINE_FOUND)
   IF (Readline_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find GNU readline")
   ENDIF (Readline_FIND_REQUIRED)
ENDIF (READLINE_FOUND)
