###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-19
#
#  Copyright (C) 2014-2015 Feel++ Consortium
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
option( FEELPP_ENABLE_NLOPT "Enable NLOPT (NonLinear Optimisation Library)" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )

if ( FEELPP_ENABLE_NLOPT )
  feelppContribPrepare( nlopt )

  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    set(NLOPT_INCLUDE_DIR ${FEELPP_SOURCE_DIR}/contrib/nlopt/api  ${FEELPP_BINARY_DIR}/contrib/nlopt/api)
    target_include_directories(feelpp_contrib INTERFACE
      $<BUILD_INTERFACE:${FEELPP_SOURCE_DIR}/contrib/nlopt/api>
      $<BUILD_INTERFACE:${FEELPP_BINARY_DIR}/contrib/nlopt/api>
      $<INSTALL_INTERFACE:include/feelpp>
      )
    target_link_libraries(feelpp_contrib INTERFACE feelpp_nlopt )
    SET( FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} NLOpt/Contrib")
    SET(FEELPP_HAS_NLOPT 1)
    list(INSERT FEELPP_LIBRARIES 0 feelpp_nlopt)
  endif()
endif()
