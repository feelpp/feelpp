###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-11
#
#  Copyright (C) 2009 Universit� Joseph Fourier (Grenoble I)
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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
feelpp_add_application( thermalfin SRCS thermalfin.cpp LABELS doc GEO thermalfin.geo CFG thermalfin.cfg )
feelpp_add_application( resistance SRCS resistance.cpp LABELS doc)
feelpp_add_application( heatsink SRCS heatsink.cpp LABELS doc CFG heatsink.cfg GEO heatsink2d.geo)

feelpp_add_application( ground2DP2 SRCS ground.cpp  LABELS doc DEFS FEELPP_DIM=2 FEELPP_ORDER=2
  CFG ground.cfg ground-realistic.cfg ground-realistic-2.cfg GEO ground2D.geo  )

add_custom_target(feelpp_doc_rht)
feelpp_add_application( rht_2dp1 SRCS rht.cpp LABELS doc DEFS FEELPP_DIM=2 FEELPP_ORDER=1 )
feelpp_add_application( rht_3dp1 SRCS rht.cpp LABELS doc DEFS FEELPP_DIM=3 FEELPP_ORDER=1)
add_dependencies(feelpp_doc_rht feelpp_doc_rht_2dp1 feelpp_doc_rht_3dp1)
feelpp_add_testcase(cases PREFIX doc_heat)