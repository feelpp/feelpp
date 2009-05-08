# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2007-02-23
#
#  Copyright (C) 2007-2008 Université de Grenoble 1 (Joseph Fourier)
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
include $(top_srcdir)/admin/common.am
ACLOCAL_AMFLAGS = -I admin

SUBDIRS= admin python convergence perf
#. $(ENABLE_EXTRAMODULE_VALIDATION) convergence perf kovasznay stokes turek structure
DIST_SUBDIRS = admin python convergence perf


$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	{									   \
	echo '# Signature of the current package.'; \
	echo 'm4_define([AT_PACKAGE_NAME],		[@PACKAGE_NAME@])'; \
	echo 'm4_define([AT_PACKAGE_TARNAME],	[@PACKAGE_TARNAME@])'; \
	echo 'm4_define([AT_PACKAGE_VERSION],	[@PACKAGE_VERSION@])'; \
	echo 'm4_define([AT_PACKAGE_STRING],	[@PACKAGE_STRING@])'; \
	echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} > $(srcdir)/package.m4

TESTSUITE = $(srcdir)/benchmarksuite

EXTRA_DIST = package.m4 benchmarksuite benchmarksuite.at atlocal.in

check-local: atconfig  atlocal $(TESTSUITE)
	$(SHELL) $(TESTSUITE)

installcheck-local: atlocal atconfig  $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
              $(TESTSUITEFLAGS)

clean-local:
	test ! -f '$(TESTSUITE)' || \
     $(SHELL) '$(TESTSUITE)' --clean

AUTOM4TE    = autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/benchmarksuite.at $(srcdir)/package.m4 $(srcdir)/convergence/convergence.at $(srcdir)/perf/perf.at
	$(AUTOTEST) -I . -I $(srcdir) -I $(top_srcdir) $@.at -o $@.tmp
	mv $@.tmp $@

# benchmarksuite: $(TESTSUITE)

atconfig: $(top_builddir)/config.status
	cd $(top_builddir) && \
      $(SHELL) ./config.status $@

atlocal: $(srcdir)/atlocal.in $(top_builddir)/config.status
	cd $(top_builddir) && \
      $(SHELL) ./config.status $@

MOSTLYCLEANFILES	= benchmarksuite.log `find -type f -print | grep benchmarksuite.dir`

DISTCLEANFILES		= atconfig atlocal summary


