###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2013-10-14
#
#  Copyright (C) 2013 Feel++ Consortium
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
feelpp_add_application( doc_ge_sm SRCS ge_sm.cpp CFG ge_sm.cfg  )
feelpp_add_application( doc_ge_laplacian SRCS ge_laplacian.cpp CFG ge_laplacian.cfg  )
feelpp_add_application( doc_ge_laplacianv SRCS ge_laplacianv.cpp CFG ge_laplacianv.cfg  DEFS USE_LM)
feelpp_add_application( doc_ge_curlcurl SRCS ge_curlcurl.cpp CFG ge_curlcurl.cfg)
feelpp_add_application( doc_ge_stokes SRCS ge_stokes.cpp CFG ge_stokes.cfg)
feelpp_add_application( doc_ge_stokes_curl SRCS ge_stokes_curl.cpp CFG ge_stokes_curl.cfg)

file(COPY cylinder.geo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
