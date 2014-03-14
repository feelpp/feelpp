###  TEMPLATE.txt.tpl; coding: utf-8 ---
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2014-01-19
#
#  Copyright (C) 2011-2014 Feel++ Consortium
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
    configure_file(benchmarkgrepl.xml.in ${CMAKE_CURRENT_SOURCE_DIR}/benchmarkgrepl${O}.xml.in COPYONLY)
    configure_file(benchmarkgrepl.cfg.in ${CMAKE_CURRENT_SOURCE_DIR}/benchmarkgrepl${O}.cfg )
    crb_add_model(benchmarkgrepl${O} BenchmarkGrepl CLASS BenchmarkGrepl<${O}> HDRS benchmarkgrepl.hpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG benchmarkgrepl${O}.cfg XML benchmarkgrepl${O}.xml )
endforeach()
