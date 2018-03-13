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
option( FEELPP_ENABLE_KWSYS "Enable KWSys (Kitware system Library)" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION})

if( FEELPP_ENABLE_KWSYS )
  feelppContribPrepare( kwsys )

  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    #add_subdirectory(contrib/ipopt)
    #    SET(KWSYS_INCLUDE_DIR
    #      ${FEELPP_SOURCE_DIR}/contrib/ipopt/Ipopt/src/Interfaces
    #      ${FEELPP_BINARY_DIR}/contrib/ipopt/include/
    #      )
    #
    #    # Compile/copy header in cmake binary dirs.
    #    include_directories(${KWSYS_INCLUDE_DIR})
    SET(KWSYS_NAMESPACE feelpp_kwsys)
    #  # Enable all components.
    #  set(KWSYS_ENABLE_C 1)
    #  set(KWSYS_USE_Base64 1)
    #  set(KWSYS_USE_Directory 1)
    #  set(KWSYS_USE_DynamicLoader 1)
    #  set(KWSYS_USE_Encoding 1)
    #  set(KWSYS_USE_Glob 1)
    #  set(KWSYS_USE_MD5 1)
    #  set(KWSYS_USE_Process 1)
    #  set(KWSYS_USE_RegularExpression 1)
    set(KWSYS_USE_System 1)
    #  set(KWSYS_USE_SystemTools 1)
    #  set(KWSYS_USE_CommandLineArguments 1)
    #  set(KWSYS_USE_Terminal 1)
    #  set(KWSYS_USE_IOStream 1)
    #  set(KWSYS_USE_FStream 1)
    #  set(KWSYS_USE_String 1)
    set(KWSYS_USE_SystemInformation 1)
    #  set(KWSYS_USE_ConsoleBuf 1)

    set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} KWSys/Contrib" )
    set(FEELPP_HAS_KWSYS 1)

    # install rules.
    set(KWSYS_INSTALL_BIN_DIR bin)
    set(KWSYS_INSTALL_LIB_DIR lib)
    set(KWSYS_INSTALL_INCLUDE_DIR include)
    set(KWSYS_INSTALL_COMPONENT_NAME_RUNTIME Runtime)
    set(KWSYS_INSTALL_COMPONENT_NAME_DEVELOPMENT Development)
    set(KWSYS_INSTALL_EXPORT_NAME feelpp-contrib-export-targets)

    list(INSERT FEELPP_LIBRARIES 0 feelpp_kwsys)

    endif() # Contrib prepare succeed
endif() # kwsys enable

