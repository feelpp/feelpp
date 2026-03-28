###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 07 Apr 2019
#
#  Copyright (C) 2019 Feel++ Consortium
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
add_custom_target( benchs )

feelpp_add_application( env SRCS 01-env.cpp LINK_LIBRARIES  benchmark EXEC B_ENV )
feelpp_add_application( mesh SRCS 02-mesh.cpp LINK_LIBRARIES  benchmark EXEC B_MESH )
feelpp_add_application( space SRCS 03-space.cpp LINK_LIBRARIES  benchmark EXEC B_SPACE )
feelpp_add_application( integrate SRCS 04-integrate.cpp LINK_LIBRARIES  benchmark EXEC B_INTEGRATE )
feelpp_add_application( basis SRCS 05-basis.cpp LINK_LIBRARIES  benchmark EXEC B_BASIS )
feelpp_add_application( tensor_basis SRCS 06-tensor-basis.cpp LINK_LIBRARIES  benchmark EXEC B_TENSOR_BASIS )
feelpp_add_application( voigt SRCS 07-voigt.cpp LINK_LIBRARIES  benchmark EXEC B_VOIGT )
feelpp_add_application( voigt_elasticity SRCS 08-voigt-elasticity.cpp LINK_LIBRARIES  benchmark EXEC B_VOIGT_ELASTICITY )
feelpp_add_application( matvec SRCS 09-matvec.cpp LINK_LIBRARIES  benchmark EXEC B_MATVEC )
if ( Python3_EXECUTABLE )
  add_custom_target( voigt_compile_bench
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/10-voigt-compile-bench.py ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running Voigt compile-time benchmark" )
endif()
feelpp_add_application( gmc SRCS gmc.cpp LINK_LIBRARIES  benchmark EXEC B_GMC )
if ( FEELPP_ENABLE_SPECX )
    feelpp_add_application( specx SRCS 99-specx.cpp LINK_LIBRARIES Feelpp::feelpp_specx )
endif()
add_dependencies( benchs ${B_ENV} ${B_MESH} ${B_SPACE} ${B_INTEGRATE} ${B_BASIS} ${B_TENSOR_BASIS} ${B_VOIGT} ${B_VOIGT_ELASTICITY} ${B_MATVEC} ${B_GMC})
