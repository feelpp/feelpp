###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-04-12
#
#  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
execute_process(COMMAND uname -n OUTPUT_VARIABLE FEELPP_MACHINE_NAME )
STRING(STRIP "${FEELPP_MACHINE_NAME}" FEELPP_MACHINE_NAME )
message(STATUS "[Feel++] detected machine: ${FEELPP_MACHINE_NAME}")
#message(STATUS "[Feel++] script: ${FEELPP_SOURCE_DIR}/cmake/machines/feelpp.machines.${FEELPP_MACHINE_NAME}.cmake")
if ( EXISTS ${FEELPP_SOURCE_DIR}/cmake/machines/feelpp.machines.${FEELPP_MACHINE_NAME}.cmake )
  message( STATUS "[Feel++] Configuration found for : ${FEELPP_MACHINE_NAME}" )
  include( feelpp.machines.${FEELPP_MACHINE_NAME} )
endif()
# try harder by looking elsewhere to ensure we are on the proper machine
# we are more specific now
STRING(REGEX MATCH "login.*" FEELPP_NAME_LOGIN ${FEELPP_MACHINE_NAME} )
message(STATUS "FEELPP_MACHINE_NAME: ${FEELPP_NAME_LOGIN}")
if( FEELPP_NAME_LOGIN AND EXISTS /lrz )
  if ( EXISTS ${FEELPP_SOURCE_DIR}/cmake/machines/feelpp.machines.lrz.cmake )
    message( STATUS "[Feel++] Configuration found for : lrz(supermuc)" )
    include( feelpp.machines.lrz )
  endif()
endif()

STRING(REGEX MATCH "fen.*" FEELPP_NAME_LOGIN ${FEELPP_MACHINE_NAME} )
if( FEELPP_NAME_LOGIN AND EXISTS /cineca )
  message(STATUS "use cineca config")
  if ( EXISTS ${FEELPP_SOURCE_DIR}/cmake/machines/feelpp.machines.cineca.cmake )
    message( STATUS "[Feel++] Configuration found for : cineca(fermi)" )
    include( feelpp.machines.cineca )
  endif()
endif()

STRING(REGEX MATCH "curie.*" FEELPP_NAME_LOGIN ${FEELPP_MACHINE_NAME} )
if( FEELPP_NAME_LOGIN AND EXISTS /ccc )
  message(STATUS "use curie config")
  if ( EXISTS ${FEELPP_SOURCE_DIR}/cmake/machines/feelpp.machines.curie.cmake )
    message( STATUS "[Feel++] Configuration found for : ccc(curie)" )
    include( feelpp.machines.curie )
  endif()
endif()

STRING(REGEX MATCH "turing.*" FEELPP_NAME_LOGIN ${FEELPP_MACHINE_NAME} )
if( FEELPP_NAME_LOGIN AND EXISTS /bglocal )
  message(STATUS "use turing config")
  if ( EXISTS ${FEELPP_SOURCE_DIR}/cmake/machines/feelpp.machines.turing.cmake )
    message( STATUS "[Feel++] Configuration found for : idris(turing)" )
    include( feelpp.machines.turing )
  endif()
endif()


if( FEELPP_ENABLE_HOMEBREW AND EXISTS /usr/local/bin/brew )
  if ( EXISTS ${FEELPP_SOURCE_DIR}/cmake/machines/feelpp.machines.homebrew.cmake )
    message( STATUS "[Feel++] Configuration found for : Homebrew" )
    include( feelpp.machines.homebrew )
  endif()
endif()
