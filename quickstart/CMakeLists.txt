###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2012-12-16
#
#  Copyright (C) 2012 Feel++ Consortium
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
add_custom_target(quickstart)

feelpp_add_application( qs_laplacian SRCS qs_laplacian.cpp INCLUDE_IN_ALL CFG qs_laplacian.cfg)
install(TARGETS feelpp_qs_laplacian DESTINATION bin/ COMPONENT Bin )

feelpp_add_application( ls_laplacian SRCS ls_laplacian.cpp INCLUDE_IN_ALL CFG ls_laplacian.cfg TIMEOUT 60)
install(TARGETS feelpp_ls_laplacian DESTINATION bin/ COMPONENT Bin )

feelpp_add_application( qs_stokes SRCS qs_stokes.cpp INCLUDE_IN_ALL CFG qs_stokes.cfg )
install(TARGETS feelpp_qs_stokes DESTINATION bin/ COMPONENT Bin )

feelpp_add_application( qs_ns SRCS qs_ns.cpp INCLUDE_IN_ALL CFG qs_ns.cfg  )
install(TARGETS feelpp_qs_ns DESTINATION bin/ COMPONENT Bin )


add_dependencies(quickstart feelpp_qs_laplacian)
add_dependencies(quickstart feelpp_qs_stokes)
add_dependencies(quickstart feelpp_qs_ns)

add_custom_target(quickstart-check COMMAND "ctest" ARGS "-VV" "-R" "qs")
add_dependencies(quickstart-check quickstart)
add_dependencies(check quickstart-check)
