# -*- Mode : makefile -*-
#
# Copyright (C) 2004 EPFL, INRIA and Politechnico di Milano
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
include $(top_srcdir)/admin/common.am


lib_LTLIBRARIES		= liblifecore.la

liblifecore_la_SOURCES	= info.cpp \
application.cpp \
xmlparser.hpp \
applicationxml.cpp \
about.cpp  debug.cpp \
typeinfo.cpp\
smartassert.cpp $(MPFR_SOURCES) \
pslogger.cpp

liblifecore_la_LIBADD	= $(aztec_libs) $(BOOST_MPI_LIB) $(BOOST_SERIALIZATION_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB)
liblifecore_la_LDFLAGS  = -version-info 0:0:0

# if HAVE_PETSC
# liblifecore_la_SOURCES += applicationpetsc.cpp
# liblifecore_la_LIBADD += $(petsc_libs)
# endif

# if HAVE_TRILINOS
# liblifecore_la_SOURCES += applicationepetra.cpp
# endif


if HAVE_MPFR
liblifecore_la_SOURCES += mpfr.cpp
endif

#bin_PROGRAMS= t

#t_SOURCES   = t.cpp
#t_LDADD     = liblifecore.la

#BUILT_SOURCES = life.hpp.gch
#headers = life.hpp

## FIXME: need a better way to convince libtool to let us do this.
#life.hpp.gch: $(headers)
#	-rm -f life.hpp.gch
#	$(CXXCOMPILE) -fPIC -DPIC -x c++-header -o life.hpp.gch $(srcdir)/life.hpp

