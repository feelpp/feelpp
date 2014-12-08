###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-03-15
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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
feelpp_add_application( dofpoints SRCS dofpoints.cpp LABELS doc)
feelpp_add_application( meshho SRCS meshho.cpp LABELS doc)
feelpp_add_application( bdfpod SRCS bdfpod.cpp LABELS doc)
feelpp_add_application( perf_int SRCS perf.cpp LABELS doc )
feelpp_add_application( pmesh SRCS pmesh.cpp CFG pmesh.cfg LABELS doc)
feelpp_add_application( dofneighbors SRCS dofneighbors.cpp LABELS doc)
feelpp_add_application( dofboundary SRCS dofboundary.cpp LABELS doc)
feelpp_add_application( forms SRCS forms.cpp LABELS doc)
