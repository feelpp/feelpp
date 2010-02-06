# -*- mode: cmake -*-
#
#  This file is part of the OPUS Project
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-02-05
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
FIND_PACKAGE(Threads)
CHECK_INCLUDE_FILE_CXX(pthread.h HAVE_PTHREAD_H )
if ( HAVE_PTHREAD_H )
  add_definitions( -DHAVE_PTHREAD_H )
  set( HAVE_PTHREAD_H 1 )
endif( HAVE_PTHREAD_H )

FIND_PACKAGE(LibXml2 REQUIRED)
SET(CMAKE_REQUIRED_FLAGS "${LIBXML2_INCLUDE_DIR};${CMAKE_REQUIRED_FLAGS}")

FIND_PATH(OT_INCLUDE_DIR
  OT.hxx
  PATHS /opt/openturns/include/ /usr/include/openturns
  DOC "Directory where OpenTURNS header files are stored" )
MARK_AS_ADVANCED( OT_INCLUDE_DIR )

if( OT_INCLUDE_DIR )

  FIND_PATH(OT_WRAPPERS_DIR
    wrapper.xml
    PATHS /opt/openturns/wrappers/ /usr/lib/openturns/wrappers
    DOC "Directory where OpenTURNS wrappers are stored" )
  MARK_AS_ADVANCED( WRAPPERS_DIR )
  FIND_PATH(OT_WRAPPER_DTD_PATH
    wrapper.dtd
    PATHS /opt/openturns/wrappers/ /usr/lib/openturns/wrappers
    DOC "Directory where OpenTURNS dtd for the wrappers  is stored" )
  MARK_AS_ADVANCED( OT_WRAPPER_DTD_PATH )

  FIND_LIBRARY(OT_LIB   NAMES OT     PATHS /usr/lib /opt/openturns/lib  /usr/lib/openturns)
  FIND_PATH(OT_LIBRARY_PATH   ${OT_LIB}     PATHS /usr/lib /opt/openturns/lib  /usr/lib/openturns)

  SET(OT_LIBRARIES ${OT_LIB} )
  MARK_AS_ADVANCED(OT_LIBRARIES OT_LIB)

  SET(CMAKE_REQUIRED_FLAGS "${OT_LIBRARIES};${CMAKE_REQUIRED_FLAGS}")
  SET(CMAKE_REQUIRED_INCLUDES "${OT_INCLUDE_DIR};${CMAKE_REQUIRED_INCLUDES}")
  MESSAGE( STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")


  CHECK_INCLUDE_FILE_CXX(OT.hxx HAVE_OT_HXX )
  CHECK_INCLUDE_FILE_CXX(WrapperCommon.h HAVE_WRAPPER_COMMON_H )
endif( OT_INCLUDE_DIR )

if (OT_INCLUDE_DIR AND OT_LIB)
  set(OT_FOUND TRUE)
endif (OT_INCLUDE_DIR AND OT_LIB)

if (OT_FOUND)
  if (NOT OT_FIND_QUIETLY)
    message(STATUS "Found OT: ${OT_LIBRARIES}")
  endif (NOT OT_FIND_QUIETLY)
else (OT_FOUND)
  if (OT_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find OT")
  endif (OT_FIND_REQUIRED)
endif (OT_FOUND)


