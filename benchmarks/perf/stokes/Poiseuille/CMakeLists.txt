###  TEMPLATE.txt.tpl; coding: utf-8 ---
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-02-03
#
#  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
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

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(gen_feel_perf_stokes_Poiseuille ALL )

  feelpp_add_application( perf_stokes_poiseuille_DD_2D SRCS stokes_poiseuille.cpp LABELS perf DEFS FEELPP_COND=DD STOKESPRESSMESHTYPE=1 NO_TEST)
  feelpp_add_application( perf_stokes_poiseuille_DD_3D SRCS stokes_poiseuille.cpp LABELS perf DEFS FEELPP_COND=DD STOKESPRESSMESHTYPE=2 NO_TEST)
  feelpp_add_application( perf_stokes_poiseuille_DD_2D_lm SRCS stokes_poiseuille.cpp LABELS perf DEFS FEELPP_COND=DD FEELPP_USE_LM STOKESPRESSMESHTYPE=1 NO_TEST)
  feelpp_add_application( perf_stokes_poiseuille_DD_3D_lm SRCS stokes_poiseuille.cpp LABELS perf DEFS FEELPP_COND=DD FEELPP_USE_LM STOKESPRESSMESHTYPE=2 NO_TEST)


  feelpp_add_application( perf_stokes_poiseuille_DN_2D SRCS stokes_poiseuille.cpp LABELS perf DEFS FEELPP_COND=DN STOKESPRESSMESHTYPE=1 NO_TEST)
  feelpp_add_application( perf_stokes_poiseuille_DN_3D SRCS stokes_poiseuille.cpp LABELS perf DEFS FEELPP_COND=DN STOKESPRESSMESHTYPE=2 NO_TEST)

  feelpp_add_application( perf_stokes_poiseuille_NN_2D SRCS stokes_poiseuille.cpp LABELS perf DEFS FEELPP_COND=NN STOKESPRESSMESHTYPE=1 NO_TEST)
  feelpp_add_application( perf_stokes_poiseuille_NN_3D SRCS stokes_poiseuille.cpp LABELS perf DEFS FEELPP_COND=NN STOKESPRESSMESHTYPE=2 NO_TEST)

 #feelpp_add_application( perf_stokes_poiseuille_test_2D SRCS stokes_poiseuille_test.cpp LABELS perf DEFS STOKESPRESSMESHTYPE=1 NO_TEST)
  feelpp_add_application( perf_stokes_poiseuille_test_3D SRCS stokes_poiseuille_test.cpp LABELS perf DEFS STOKESPRESSMESHTYPE=2 NO_TEST)


#  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#  foreach(PORDER 1 2 3)
#    foreach(GORDER 1 2 3)
#      file(WRITE stokesdd_P${PORDER}G${GORDER}.cpp "#include <stokes.hpp>
#int
#main( int argc, char** argv )
#{

#    using namespace Feel;

#    Environment env( _argc=argc, _argv=argv,
#                     _desc=makeOptions(),
#                     _about=about(_name=\"stokes_Dirichlet_Dirichlet\",
#                                  _author=\"Christophe Prud'homme\",
#                                  _email=\"christophe.prudhomme@feelpp.org\") );
#    Feel::Stokes_Dirichlet_Dirichlet<${PORDER},${GORDER}> stokesDD;
#    stokesDD.run();
#}
#" )
#      feelpp_add_application( stokesdd_P${PORDER}G${GORDER} SRCS stokesdd_P${PORDER}G${GORDER}.cpp LABELS doc DEFS STOKESPRESSMESHTYPE=2 NO_TEST )
#    endforeach()
#  endforeach()