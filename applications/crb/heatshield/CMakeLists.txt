###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-04-23
#
#  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
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
foreach(O 1 2 3)

    configure_file(heatshield.xml.in ${CMAKE_CURRENT_SOURCE_DIR}/heatshield${O}.xml.in COPYONLY)
    configure_file(heatshield.cfg.in ${CMAKE_CURRENT_SOURCE_DIR}/heatshield${O}.cfg )
    crb_add_model(heatshield${O} HeatShield CLASS HeatShield<${O}> HDRS heatshield.hpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG heatshield${O}.cfg  XML heatshield${O}.xml )

    configure_file(heatshieldminimalversion.xml.in ${CMAKE_CURRENT_SOURCE_DIR}/heatshieldminimalversion${O}.xml.in COPYONLY)
    configure_file(heatshieldminimalversion.cfg.in ${CMAKE_CURRENT_SOURCE_DIR}/heatshieldminimalversion${O}.cfg )
    crb_add_model(heatshieldminimalversion${O} HeatShieldMinimalVersion CLASS HeatShieldMinimalVersion<${O}> HDRS heatshield-minimalversion.hpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG heatshieldminimalversion${O}.cfg  GEO heatshield.geo XML heatshieldminimalversion${O}.xml )

endforeach()
