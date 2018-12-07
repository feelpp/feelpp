###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
#       Date: 2014-12-18
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

#
# JSONLAB
#
OPTION( FEELPP_ENABLE_JSONLAB "Enable JSONLAB" ON )

if ( FEELPP_ENABLE_JSONLAB )
  feelppContribPrepare( jsonlab )

  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    SET(FEELPP_HAS_JSONLAB 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} JSONLAB" )
    ADD_DEFINITIONS( -DFEELPP_HAS_JSONLAB )
  else()
      MESSAGE(WARNING "JSONLAB NOT FOUND!")
  endif()

endif()
