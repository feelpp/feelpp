# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2008-11-11
#
#  Copyright (C) 2008 Université Joseph Fourier
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

SUBDIRS=styles pngs pdfs

ACLOCAL_AMFLAGS = -I admin

bin_PROGRAMS=life_myapp life_mymesh life_myintegrals life_laplacian life_stokes life_bratu life_nonlinearpow  life_convection

life_nonlinearpow_SOURCES=nonlinearpow.cpp
life_nonlinearpow_LDADD=$(LIFE_LIBS)

life_bratu_SOURCES=bratu.cpp
life_bratu_LDADD=$(LIFE_LIBS)

life_stokes_SOURCES=stokes.cpp
life_stokes_LDADD=$(LIFE_LIBS)

life_laplacian_SOURCES=laplacian.cpp
life_laplacian_LDADD=$(LIFE_LIBS)

life_myintegrals_SOURCES=myintegrals.cpp
life_myintegrals_LDADD=$(LIFE_LIBS)

life_mymesh_SOURCES=mymesh.cpp
life_mymesh_LDADD=$(LIFE_LIBS)

life_convection_SOURCES= convection.hpp convection_other.cpp convection_main.cpp convection_run.cpp convection_residual.cpp convection_jacobian.cpp convection_jacobian1.cpp convection_jacobian2.cpp convection_lo.cpp convection_lo2.cpp
life_convection_LDADD=$(LIFE_LIBS)

life_myapp_SOURCES=myapp.cpp
life_myapp_LDADD=$(LIFE_LIBS)

SUFFIXES=.cpp .hpp .tex .pdf

CLEANFILES = *.aux *.log *.lot *.lof *.toc *.bbl *.blg *.out *.idx *.ind *.ilg $(doc_DATA) version.tex

.tex.pdf:
	echo "*** INFO: first pdflatex compilation" ; \
	TEXINPUTS=:.//:$(top_builddir)/tutorial//:$(srcdir)//: $(PDFLATEX) $< ; \
	echo "*** INFO: index generation" ; \
	if test -f $*.idx; then $(MAKEINDEX) -s $(srcdir)/life-manual.ist $*.idx ; fi ; \
	echo "*** INFO: bibliography creation" ; \
	if test -f $(srcdir)/*.bib; then BIBINPUTS=:$(srcdir): $(BIBTEX) $*.aux ; fi ; \
	echo "*** INFO: table of content insertion" ; \
	TEXINPUTS_pdflatex=:.//:$(top_builddir)/tutorial//:$(srcdir)//: $(PDFLATEX) $< ; \
	echo "*** INFO: labels insertion" ; \
	TEXINPUTS_pdflatex=:.//:$(top_builddir)/tutorial//:$(srcdir)//: $(PDFLATEX) $< -o $(top_builddir)/$*.pdf



EXTRA_DIST = life-manual.tex version.tex.in feel.tex gfdl.tex life-manual.ist version.tex.in example-heatns.tex


life-manual.pdf : $(EXTRA_DIST) myapp.cpp mymesh.cpp myintegrals.cpp laplacian.cpp stokes.cpp bratu.cpp

pdf-local : life-manual.pdf
doc_DATA = life-manual.pdf
