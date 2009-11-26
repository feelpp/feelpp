# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2007-08-07
#
#  Copyright (C) 2007 Université Joseph Fourier
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
check_PROGRAMS=solvepbeq test_atom

check_HEADERS=atom.hpp  heavysidefunction.hpp  molecule.hpp  mymesh.hpp  pbeqapplication.hpp  pbeq.hpp  pbeqspace.hpp

test_atom_SOURCES=atom.cpp test_atom.cpp
test_atom_LDADD=$(LIFE_LIBS)

solvepbeq_SOURCES=mymesh.cpp pbeqspace.cpp atom.cpp molecule.cpp heavysidefunction.cpp pbeq.cpp solvepbeq.cpp
solvepbeq_LDADD=$(LIFE_LIBS)

link:
	-ln -s $(srcdir)/molecule.pqr $(srcdir)/atom.pqr $(srcdir)/receptor.pqr $(srcdir)/ligand.pqr .


EXTRA_DIST=pbeq.at