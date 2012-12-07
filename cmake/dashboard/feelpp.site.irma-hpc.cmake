###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-05-03
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
set(OS_VERSION debian-sid)
set(ARCH x86_64)
set(WORK_DIR /scratch/prudhomm/dashboard/)
set(MAKE_ARGS "-j6")
set(PARALLEL "6")
set(FEELPP_WORK_DIR ${WORK_DIR})
set(FEELPP_ENABLE_CRB_ALL ON)
set(FEELPP_ENABLE_BENCHMARKS ON)
set(FEELPP_MAKE_ARGS ${MAKE_ARGS})
set(CTEST_BUILD_FLAGS -j${PARALLEL})
set(CTEST_PARALLEL_LEVEL ${PARALLEL})


