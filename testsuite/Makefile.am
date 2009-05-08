# -*- makefile; indent-tabs-mode: t -*-
#
#
#
# Copyright (C) 2004 Christophe Prud'homme (christophe.prudhomme@epfl.ch)
#
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
#

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

check_HEADERS=testsuite.hpp

SUBDIRS		= lifecore lifealg lifemesh lifepoly lifediscr lifefilters lifematerial

EXTRA_DIST=testsuite.at


$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	{									   \
	echo '# Signature of the current package.'; \
	echo 'm4_define([AT_PACKAGE_NAME],		[@PACKAGE_NAME@])'; \
	echo 'm4_define([AT_PACKAGE_TARNAME],	[@PACKAGE_TARNAME@])'; \
	echo 'm4_define([AT_PACKAGE_VERSION],	[@PACKAGE_VERSION@])'; \
	echo 'm4_define([AT_PACKAGE_STRING],	[@PACKAGE_STRING@])'; \
	echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} > $(srcdir)/package.m4

EXTRA_DIST = package.m4 testsuite testsuite.at atlocal.in

TESTSUITE = $(srcdir)/testsuite

check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) $(TESTSUITE)

installcheck-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
              $(TESTSUITEFLAGS)

clean-local:
	test ! -f '$(TESTSUITE)' || \
     $(SHELL) '$(TESTSUITE)' --clean

TESTSUITES=testsuite.at \
           $(srcdir)/lifecore/testsuite.at\
           $(srcdir)/lifepoly/testsuite.at\
           $(srcdir)/lifealg/testsuite.at\
           $(srcdir)/lifemesh/testsuite.at\
           $(srcdir)/lifediscr/testsuite.at\
           $(srcdir)/lifefilters/testsuite.at


AUTOM4TE    = autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4 $(TESTSUITES)
	$(AUTOTEST) -I . -I $(srcdir) -I $(top_srcdir) $@.at -o $@.tmp
	mv $@.tmp $@


MOSTLYCLEANFILES	= testsuite.log `find -type f -print | grep testsuite.dir`

DISTCLEANFILES		= atconfig atlocal

subdir=testsuite
atconfig: $(top_builddir)/config.status
	cd $(top_builddir) && \
      $(SHELL) ./config.status $(subdir)/$@

atlocal: $(srcdir)/atlocal.in $(top_builddir)/config.status
	cd $(top_builddir) && \
      $(SHELL) ./config.status $(subdir)/$@



include $(top_srcdir)/admin/dist-files.mk