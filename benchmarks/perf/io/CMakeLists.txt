###  CMakeLists.txt ---

#  Author(s): Yacine OUld Rouis <yacine.ouldrouis@cea.fr>
#       Date: 2015-06-22
#
#  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)
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

add_executable( feel_perf_io
bench_perf_io.cpp
  )

target_link_libraries( feel_perf_io feelpp ${FEELPP_LIBRARIES} )
feelpp_add_application( perf_io_vector SRCS bench_perf_io_vector.cpp )
