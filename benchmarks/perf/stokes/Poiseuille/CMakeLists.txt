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

  feelpp_add_application( doc_stokes_Dirichlet_Dirichlet SRCS stokes_Dirichlet_Dirichlet.cpp LABEL benchmarks NO_TEST )
  feelpp_add_application( doc_stokes_Dirichlet_Dirichlet_lm SRCS stokes_Dirichlet_Dirichlet.cpp LABEL benchmarks DEFS FEELPP_USE_LM NO_TEST )
  feelpp_add_application( doc_stokes_Dirichlet_Dirichlet_cylinder SRCS stokes_Dirichlet_Dirichlet.cpp LABEL benchmarks DEFS STOKESPRESSMESHTYPE=2 NO_TEST )
  feelpp_add_application( doc_stokes_Dirichlet_Dirichlet_hypercube SRCS stokes_Dirichlet_Dirichlet.cpp LABEL benchmarks DEFS STOKESPRESSMESHTYPE=1 NO_TEST )

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  foreach(PORDER 1 2 3)
    foreach(GORDER 1 2 3)
      file(WRITE stokesdd_P${PORDER}G${GORDER}.cpp "#include <stokes.hpp>
int
main( int argc, char** argv )
{

    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name=\"stokes_Dirichlet_Dirichlet\",
                                  _author=\"Christophe Prud'homme\",
                                  _email=\"christophe.prudhomme@feelpp.org\") );
    Feel::Stokes_Dirichlet_Dirichlet<${PORDER},${GORDER}> stokesDD;
    stokesDD.run();
}
" )
      feelpp_add_application( stokesdd_P${PORDER}G${GORDER} SRCS stokesdd_P${PORDER}G${GORDER}.cpp LABEL benchmarks DEFS STOKESPRESSMESHTYPE=2 NO_TEST )
    endforeach()
  endforeach()

  feelpp_add_application( doc_stokes_Dirichlet_Neumann SRCS stokes_Dirichlet_Neumann.cpp LABEL benchmarks NO_TEST )
  feelpp_add_application( doc_stokes_Dirichlet_Neumann_lm SRCS stokes_Dirichlet_Neumann.cpp LABEL benchmarks DEFS FEELPP_USE_LM NO_TEST )
  feelpp_add_application( doc_stokes_Dirichlet_Neumann_cylinder SRCS stokes_Dirichlet_Neumann.cpp LABEL benchmarks DEFS STOKESPRESSMESHTYPE=2 NO_TEST)
  feelpp_add_application( doc_stokes_Dirichlet_Neumann_hypercube SRCS stokes_Dirichlet_Neumann.cpp LABEL benchmarks DEFS STOKESPRESSMESHTYPE=1 NO_TEST)

  feelpp_add_application( doc_stokes_Neumann_Neumann SRCS stokes_Neumann_Neumann.cpp LABEL benchmarks NO_TEST )
  #feelpp_add_application( doc_stokes_Neumann_Neumann_aorta SRCS stokes_Neumann_Neumann.cpp LABEL benchmarks DEFS STOKESPRESSMESHTYPE=3 NO_TEST )
  feelpp_add_application( doc_stokes_Neumann_Neumann_cylinder SRCS stokes_Neumann_Neumann.cpp LABEL benchmarks DEFS STOKESPRESSMESHTYPE=2 NO_TEST )
  feelpp_add_application( doc_stokes_Neumann_Neumann_hypercube SRCS stokes_Neumann_Neumann.cpp LABEL benchmarks DEFS STOKESPRESSMESHTYPE=1 NO_TEST)


