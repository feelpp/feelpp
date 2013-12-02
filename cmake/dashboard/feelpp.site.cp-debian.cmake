###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-05-03
#
#  Copyright (C) 2012 Universit� Joseph Fourier (Grenoble I)
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
set(OS_VERSION debian-wheezy)
set(ARCH x86_64)
set(WORK_DIR /home/vhuber)
set(GCC_MAKE_ARGS "-j1")
set(GCC_PARALLEL "1")
set(CLANG_MAKE_ARGS "-j4")
set(CLANG_PARALLEL "4")
set(BUILD_TYPE "release")
set(FEELPP_WORK_DIR ${WORK_DIR})

#Directories to update at run time
set(FEELPP_MODULES "research/hifimagnet" "research/fluid" )

#Options
set(FEELPP_ENABLE_ALTIVEC OFF)
set(FEELPP_ENABLE_BUILD_STATIC OFF)
set(FEELPP_ENABLE_DOXYGEN OFF)
set(FEELPP_ENABLE_NEON OFF)
set(FEELPP_ENABLE_OPENTURNS ON)
set(FEELPP_ENABLE_PCH_FOR_APPLICATIONS ON)
set(FEELPP_ENABLE_VERBOSE_CMAKE OFF)

#Directories
set(FEELPP_BENCHMARK_FLAG OFF) #
set(FEELPP_ENABLE_TESTS ON) #testsuite
set(FEELPP_ENABLE_DOCUMENTATION OFF) #doc
set(FEELPP_ENABLE_BENCHMARKS ON) #benchmarkes
set(FEELPP_ENABLE_RESEARCH OFF) #research
set(FEELPP_ENABLE_APPLICATIONS OFF) #applications
set(FEELPP_ENABLE_CRB_ALL ON) #Applications/CRB
set(FEELPP_ENABLE_APPLICATIONS_CRB ON)#Applications/CRB

