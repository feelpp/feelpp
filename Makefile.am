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

if ENABLE_LIFESUITE
MAYBE_LIFESUITE=testsuite
endif

SUBDIRS = boost contrib life $(MAYBE_LIFESUITE) .
DIST_SUBDIRS = boost contrib life $(MAYBE_LIFESUITE)  .
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



EXTRA_DIST = Makefile.headers lifesuite.at

include $(top_srcdir)/Makefile.headers


if MAINTAINER_MODE
ACLOCAL_AMFLAGS=-I admin

DISTCHECK_CONFIGURE_FLAGS=--disable-testsuite --disable-examples --disable-benchmarks --enable-opt2 CXXFLAGS=""

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


if ENABLE_LIFESUITE

TESTSUITES=$(srcdir)/testsuite/testsuite.at\
$(srcdir)/examples/examplesuite.at\
$(srcdir)/benchmarks/benchmarksuite.at



package.m4: $(top_srcdir)/configure.ac
	{									   \
	echo '# Signature of the current package.'; \
	echo 'm4_define([AT_PACKAGE_NAME],		[@PACKAGE_NAME@])'; \
	echo 'm4_define([AT_PACKAGE_TARNAME],	[@PACKAGE_TARNAME@])'; \
	echo 'm4_define([AT_PACKAGE_VERSION],	[@PACKAGE_VERSION@])'; \
	echo 'm4_define([AT_PACKAGE_STRING],	[@PACKAGE_STRING@])'; \
	echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} > package.m4

EXTRA_DIST = package.m4 $(TESTSUITE) atlocal.in

TESTSUITE = lifesuite

check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) $(TESTSUITE)

check-testsuite: atconfig atlocal  $(TESTSUITE) $(LN_TESTSUITE)
	$(SHELL) $(TESTSUITE)

installcheck-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
              $(TESTSUITEFLAGS)

clean-local:
	test ! -f '$(TESTSUITE)' || \
     $(SHELL) '$(TESTSUITE)' --clean

AUTOM4TE    = autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/lifesuite.at package.m4 $(TESTSUITES)
	$(AUTOTEST) -I . -I $(srcdir) $@.at -o $@.tmp
	mv $@.tmp $@
#	@if test $(top_srcdir) != $(top_builddir); then\
#     ln -sf $(TESTSUITE);\
#    fi

TESTINCLUDES=-I . -I$(srcdir) -I$(srcdir)/testsuite -I$(srcdir)/examples -I$(srcdir)/benchmarks
lifesuite: $(srcdir)/lifesuite.at package.m4 $(TESTSUITES)
	$(AUTOTEST) $(TESTINCLUDES) $(srcdir)/$@.at -o $@.tmp
	mv $@.tmp $@
#	@if test $(top_srcdir) != $(top_builddir); then\
#     ln -sf $(TESTSUITE);\
#    fi


MOSTLYCLEANFILES	= lifesuite.log `find -type f -print | grep lifesuite.dir`

DISTCLEANFILES		= atconfig lifesuite

CLEANFILES 			= package.m4

endif # ENABLE_LIFESUITE