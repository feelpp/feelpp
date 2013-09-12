### FindGooglePerfTools.cmake; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-04-11
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
# -*- cmake -*-

# Find the Google perftools includes and libraries
# This module defines
#  GOOGLE_PERFTOOLS_INCLUDE_DIR, where to find heap-profiler.h, etc.
#  GOOGLE_PERFTOOLS_FOUND, If false, do not try to use Google perftools.
# also defined for general use are
#  TCMALLOC_LIBRARIES, where to find the tcmalloc library.
#  STACKTRACE_LIBRARIES, where to find the stacktrace library.
#  PROFILER_LIBRARIES, where to find the profiler library.

FIND_PATH(GOOGLE_PERFTOOLS_INCLUDE_DIR
  gperftools/heap-profiler.h
  PATHS
  /opt/local/include
  /usr/local/include
  /usr/include
  PATH_SUFFIXES gperftools google
  )

SET(TCMALLOC_NAMES ${TCMALLOC_NAMES} tcmalloc)
FIND_LIBRARY(TCMALLOC_LIBRARY
  NAMES ${TCMALLOC_NAMES}
  PATHS /usr/lib /usr/local/lib /opt/local/lib
  )

IF (TCMALLOC_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)
  SET(TCMALLOC_LIBRARIES ${TCMALLOC_LIBRARY})
ENDIF (TCMALLOC_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)

SET(STACKTRACE_NAMES ${STACKTRACE_NAMES} stacktrace)
FIND_LIBRARY(STACKTRACE_LIBRARY
  NAMES ${STACKTRACE_NAMES}
  PATHS /opt/local/lib /usr/lib /usr/local/lib
  )

IF (STACKTRACE_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)
    SET(STACKTRACE_LIBRARIES ${STACKTRACE_LIBRARY})
ENDIF (STACKTRACE_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)

SET(PROFILER_NAMES ${PROFILER_NAMES} profiler)
FIND_LIBRARY(PROFILER_LIBRARY
  NAMES ${PROFILER_NAMES}
  PATHS /opt/local/lib /usr/lib /usr/local/lib
  )

IF (PROFILER_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)
    SET(PROFILER_LIBRARIES ${PROFILER_LIBRARY})
ENDIF (PROFILER_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)


find_package_handle_standard_args(GOOGLE_PERFTOOLS "Could not find Google PerfTools " GOOGLE_PERFTOOLS_INCLUDE_DIR TCMALLOC_LIBRARIES  PROFILER_LIBRARIES )
MARK_AS_ADVANCED(
  TCMALLOC_LIBRARY
  STACKTRACE_LIBRARY
  PROFILER_LIBRARY
  GOOGLE_PERFTOOLS_INCLUDE_DIR
  )
