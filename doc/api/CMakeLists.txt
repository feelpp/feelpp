###  TEMPLATE.txt.tpl; coding: utf-8 ---
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-02-18
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

# Set to "ON" if you want the LaTeX documentation to be built.
set(DOXYFILE_LATEX 0)

find_program(DOXYGEN doxygen)
if ( DOXYGEN )
   # doxygen is optional
   include(UseDoxygen OPTIONAL)

   # make sure that doxygen is not run byt default
   set_target_properties(doxygen PROPERTIES EXCLUDE_FROM_ALL TRUE)

   configure_file(
     ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
     ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
     )

   configure_file(
     ${CMAKE_CURRENT_SOURCE_DIR}/latexheader.tex.in
     ${CMAKE_CURRENT_BINARY_DIR}/latexheader.tex
     )
endif()

