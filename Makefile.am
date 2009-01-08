# -*- Mode : makefile; indent-tabs-mode: t -*-
#
#    AUTHOR: Christophe Prud'homme
#       ORG: EPFL
#    E-MAIL: christophe.prudhomme@ujf-grenoble.fr
#
# DESCRIPTION:
# ============
#  Distributed under the GPL(GNU Public License):
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# DESCRIP-END.
#
DO_NOT_COMPILE=polyvis
#SUBDIRS = admin $(TOPSUBDIRS)

if ENABLE_TESTS
MAYBE_TESTS=testsuite
endif

SUBDIRS = admin boost contrib life $(MAYBE_TESTS) .
DIST_SUBDIRS = admin boost contrib life $(MAYBE_TESTS)  .

if ENABLE_MODULE_DOC
  SUBDIRS += doc
  DIST_SUBDIRS += doc
endif
if ENABLE_MODULE_EXAMPLES
  SUBDIRS += examples
  DIST_SUBDIRS += examples
endif
if ENABLE_MODULE_BENCHMARKS
  SUBDIRS += benchmarks
  DIST_SUBDIRS += benchmarks
endif
if ENABLE_MODULE_APPLICATIONS
  SUBDIRS += applications
  DIST_SUBDIRS += applications
endif
if ENABLE_MODULE_RESEARCH
  SUBDIRS += research
endif





include $(top_srcdir)/Makefile.headers


EXTRA_DIST = Makefile.headers

if MAINTAINER_MODE
ACLOCAL_AMFLAGS=-I admin

DISTCHECK_CONFIGURE_FLAGS=--enable-debug=0  --enable-opt2 CXXFLAGS=""

$(nobase_include_HEADERS):
	cd $(top_srcdir) && $(SHELL) tools/update-headers $(top_srcdir)



Makefile.headers: $(top_srcdir)/tools/update-headers
	cd $(top_srcdir) && $(SHELL) tools/update-headers $(top_srcdir)

report:
	@if test -x tools/inspect/inspect; then\
        thedate=`date +%Y%m%d`;\
		thepwd=`pwd`;\
		cd $(srcdir) && $$thepwd/tools/inspect/inspect > $$thepwd/report-$$thedate.html;\
    else\
       echo "tools/inspect/inspect has not been built";\
       echo "consider typing: cd tools/inspect && make";\
	fi

endif


DISTCLEANFILES		= summary