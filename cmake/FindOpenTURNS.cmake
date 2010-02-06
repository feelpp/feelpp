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

FIND_PACKAGE(LibXml2 REQUIRED)
SET(CMAKE_REQUIRED_FLAGS "${LIBXML2_INCLUDE_DIR};${CMAKE_REQUIRED_FLAGS}")

FIND_PATH(OT_INCLUDE_DIR
  OT.hxx
  PATHS /opt/openturns/include/ /usr/include/openturns
  DOC "Directory where OpenTURNS header files are stored" )
MARK_AS_ADVANCED( OT_INCLUDE_DIR )

if( OT_INCLUDE_DIR )

  FIND_PATH(OT_WRAPPERS_DIR
    wrapper.dtd
    PATHS /opt/openturns/wrappers/ /usr/lib/openturns/wrappers
    DOC "Directory where OpenTURNS wrappers are stored" )
  MARK_AS_ADVANCED( WRAPPERS_DIR )

  FIND_LIBRARY(OT_LIB   NAMES OT     PATHS /usr/lib /opt/openturns/lib  /usr/lib/openturns)
  FIND_PATH(OT_LIBRARY_PATH   ${OT_LIB}     PATHS /usr/lib /opt/openturns/lib  /usr/lib/openturns)

  SET(OT_LIBRARIES
    ${OT_LIB}
    ${LIBXML2_LIBRARIES}
    )
  MARK_AS_ADVANCED(OT_LIBRARIES OT_LIB)

  SET(CMAKE_REQUIRED_FLAGS "${OT_LIBRARIES};${CMAKE_REQUIRED_FLAGS}")
  SET(CMAKE_REQUIRED_INCLUDES "${OT_INCLUDE_DIR};${CMAKE_REQUIRED_INCLUDES}")
  MESSAGE( STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")


  CHECK_INCLUDE_FILE_CXX(OT.hxx HAVE_OT_HXX )
  CHECK_INCLUDE_FILE_CXX(WrapperCommon.h HAVE_WRAPPER_COMMON_H )


  if ( OT_LIB )
    SET( OT_FOUND "YES")
  endif (OT_LIB)

endif( OT_INCLUDE_DIR )
