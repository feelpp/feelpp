###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 03 Mar 2016
#
#  Copyright (C) 2016 Feel++ Consortium
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
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/toolboxes)

if (FEELPP_TOOLBOXES_ENABLE_HDG)
   feelpp_add_application( electro-thermal SRCS et_main.cpp INCLUDE_IN_ALL CFG et.cfg DEFS FEELPP_DIM=2 LINK_LIBRARIES feelpp_modelcore )
   
   feelpp_add_application( electro-thermal-2 SRCS et2_main.cpp INCLUDE_IN_ALL CFG et.cfg DEFS FEELPP_DIM=2 LINK_LIBRARIES feelpp_modelcore)
   feelpp_add_application( electro-thermal-errchk SRCS et_main_errchk.cpp INCLUDE_IN_ALL CFG et.cfg DEFS FEELPP_DIM=2 LINK_LIBRARIES feelpp_modelcore)
   
   foreach(D 3)
      foreach(P RANGE 1 1 1)
   	feelpp_add_application( mixed-poisson-elasticity_${D}DP${P} SRCS mixedpoissonelasticity_main.cpp INCLUDE_IN_ALL LINK_LIBRARIES feelpp_modelcore DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P} )
      endforeach()
   endforeach()
   
   #foreach(D 3)
   #   foreach(P RANGE 1 1 1)
   	#configure_file(lc_main_model.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/lc${D}DP${P}_main.cpp )
   #	feelpp_add_application( lamina-cribrosa-model_${D}DP${P} SRCS lc_main_model.cpp INCLUDE_IN_ALL LINK_LIBRARIES feelpp_modelcore DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P} )
   #   endforeach()
   #endforeach()
   
   
      feelpp_add_application( lamina-cribrosa-model_3DP1 SRCS lc_main_model.cpp  LINK_LIBRARIES Feelpp::feelpp_modelcore Feelpp::feelpp_toolbox_hdg_poisson_lib DEFS FEELPP_DIM=3 FEELPP_ORDER=1 FEELPP_G_ORDER=1 )
endif()

feelpp_add_application( mixed-stokes SRCS mixedstokes.cpp DEFS FEELPP_DIM=2 FEELPP_ORDER=1)

feelpp_add_testcase(quarter-turn)
feelpp_add_testcase(lamina-cribrosa)
feelpp_add_testcase(kovasznay)
feelpp_add_testcase(eye2brain_doc)
