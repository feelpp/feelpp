# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
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
set(SRC
  info.cpp
  feel.cpp
  feelio.cpp
  bimap.cpp
  application.cpp
  environment.cpp
  repository.cpp
  simget.cpp
  xmlparser.hpp
  about.cpp
  typeinfo.cpp
  smartassert.cpp
  pslogger.cpp
  worldcomm.cpp
  debug.cpp
  utility.cpp
  ptreetools.cpp
  checker.cpp
  remotedata.cpp
  journalmanager.cpp
  journalwatcher.cpp
  table.cpp
  tabulateinformations.cpp
  terminalproperties.cpp
  )
if(FEELPP_HAS_MONGOCXX)
    set( SRC mongocxx.cpp ${SRC} )
endif()
if(FEELPP_HAS_HDF5)
    set(SRC hdf5.cpp ${SRC})
endif()
if (NOT FEELPP_HAS_STDIO_FMEMOPEN )
  set( SRC fmemopen.c ${SRC})
endif()
file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
feel_append_src(feelcore "${SRC};${HDR}")
