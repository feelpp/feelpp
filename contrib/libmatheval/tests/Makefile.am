# Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2011 Free
# Software Foundation, Inc.
#
# This file is part of GNU libmatheval
#
# GNU libmatheval is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# GNU libmatheval is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU libmatheval.  If not, see
# <http://www.gnu.org/licenses/>.

TESTSUITE_AT = testsuite.at basics.at bin_ops.at constants.at	\
functions.at numbers.at testsuite.at un_ops.at variables.at
TESTSUITE = $(srcdir)/testsuite

EXTRA_DIST = $(TESTSUITE_AT) testsuite

check_SCRIPTS = matheval.sh
DISTCLEANFILES = atconfig $(check_SCRIPTS)
MAINTAINERCLEANFILES = Makefile.in $(TESTSUITE)

$(srcdir)/package.m4: $(top_srcdir)/configure.in
	{					\
	  echo '# Signature of the current package.'; \
	  echo 'm4_define([AT_PACKAGE_NAME],      [@PACKAGE_NAME@])'; \
	  echo 'm4_define([AT_PACKAGE_TARNAME],   [@PACKAGE_TARNAME@])'; \
	  echo 'm4_define([AT_PACKAGE_VERSION],   [@PACKAGE_VERSION@])'; \
	  echo 'm4_define([AT_PACKAGE_STRING],    [@PACKAGE_STRING@])'; \
	  echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} >$(srcdir)/package.m4

EXTRA_DIST += package.m4

AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): package.m4 $(TESTSUITE_AT)
	$(AUTOTEST) -I $(srcdir) testsuite.at -o $@.tmp
	mv $@.tmp $@

atconfig: $(top_builddir)/config.status
	cd $(top_builddir) && ./config.status tests/$@

clean-local:
	if test -f $(TESTSUITE); then \
		$(SHELL) $(TESTSUITE) --clean; \
	fi

check-local: atconfig $(TESTSUITE)
	$(SHELL) $(TESTSUITE);

installcheck-local:
	if test -f $(TESTSUITE); then \
		$(SHELL) $(TESTSUITE) AUTOTEST_PATH=$(exec_prefix)/bin; \
	fi


INCLUDES = -I$(top_builddir) -I$(top_builddir)/lib -I$(top_srcdir)	\
-I$(top_srcdir)/lib

noinst_PROGRAMS = matheval

matheval_SOURCES = matheval.c

matheval_CFLAGS = $(GUILE_CFLAGS)
matheval_LDADD = $(top_builddir)/lib/libmatheval.la
matheval_LDFLAGS = $(GUILE_LDFLAGS)

.PHONY: beauty

beauty:
	-for file in $$(ls *.c); do indent -orig -i8 -ncdb $$file; done
	-rm -f *~
