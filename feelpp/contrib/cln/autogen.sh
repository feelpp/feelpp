#!/bin/sh
# Convenience script for regenerating all autogeneratable files that are
# omitted from the version control repository. In particular, this script
# also regenerates all aclocal.m4, config.h.in, Makefile.in, configure files
# with new versions of autoconf or automake.
#
# This script requires autoconf-2.63..2.69 and automake-1.11..1.16 in the PATH.

# Copyright (C) 2003-2019 Free Software Foundation, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Usage: ./autogen.sh

TEXINFO_VERSION=7.0.1

mkdir -p build-aux m4

# libtool
# Don't use libtoolize, as it may not be installed or may be outdated.
for f in build-aux/ltmain.sh m4/libtool.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 m4/lt~obsolete.m4; do
  { wget -nv --timeout=5 -O $f.tmp "https://git.savannah.gnu.org/gitweb/?p=gettext.git;a=blob_plain;f=${f}" \
      && mv $f.tmp $f; \
  } || rm -f $f.tmp
done

# config.{guess,sub}, gnulib module 'havelib', test-driver.diff
for f in build-aux/config.guess build-aux/config.sub \
         m4/lib-ld.m4 m4/lib-link.m4 m4/lib-prefix.m4 m4/host-cpu-c-abi.m4 build-aux/config.rpath \
         build-aux/test-driver.diff; do
  { wget -nv --timeout=5 -O $f.tmp "https://git.savannah.gnu.org/gitweb/?p=gnulib.git;a=blob_plain;f=${f}" \
      && mv $f.tmp $f; \
  } || rm -f $f.tmp
done

# autoconf archive
for f in m4/ax_cxx_compile_stdcxx.m4; do
  { wget -nv --timeout=5 -O $f.tmp "https://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=${f}" \
      && mv $f.tmp $f; \
  } || rm -f $f.tmp
done

# texinfo.tex
# The most recent snapshot of it is available in the gnulib repository.
# But this is a snapshot, with all possible dangers.
# A stable release of it is available through "automake --add-missing --copy",
# but that may be too old. So take the version which matches the latest stable
# texinfo release.
for f in texinfo.tex; do
  g="build-aux/$f"
  { wget -nv --timeout=5 -O $g.tmp "https://git.savannah.gnu.org/gitweb/?p=texinfo.git;a=blob_plain;f=doc/${f};hb=refs/tags/texinfo-${TEXINFO_VERSION}" \
      && mv $g.tmp $g; \
  } || rm -f $g.tmp
done

aclocal -I m4
autoconf
autoheader && touch autoconf/cl_config.h.in
# Make sure we get new versions of files brought in by automake.
(cd build-aux && rm -f ar-lib compile depcomp install-sh mdate-sh missing test-driver)
automake --add-missing --copy
patch build-aux/test-driver < build-aux/test-driver.diff
# Get rid of autom4te.cache directory.
rm -rf autom4te.cache
