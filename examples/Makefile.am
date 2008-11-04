# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2006-12-30
#
#  Copyright (C) 2006 Université Joseph Fourier
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
# When rebuilding, tell aclocal where to find additional macros
include $(top_srcdir)/admin/common.am

ACLOCAL_AMFLAGS = -I admin


# if ENABLE_FSI
# FSDIR=fsi
# endif

SUBDIRS      = admin generic $(FSIDIR) $(ENABLE_EXTRAMODULE_FBM) fluid solid heat levelset pbeq bloodflow # poly
DIST_SUBDIRS = admin generic $(FSIDIR) $(ENABLE_EXTRAMODULE_FBM) fluid solid heat levelset pbeq bloodflow # poly

