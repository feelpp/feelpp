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

bin_PROGRAMS=myapp mymesh myintegrals laplacian stokes bratu nonlinearpow

nonlinearpow_SOURCES=nonlinearpow.cpp
nonlinearpow_LDADD=$(LIFE_LIBS)

bratu_SOURCES=bratu.cpp
bratu_LDADD=$(LIFE_LIBS)

stokes_SOURCES=stokes.cpp
stokes_LDADD=$(LIFE_LIBS)

laplacian_SOURCES=laplacian.cpp
laplacian_LDADD=$(LIFE_LIBS)

myintegrals_SOURCES=myintegrals.cpp
myintegrals_LDADD=$(LIFE_LIBS)

mymesh_SOURCES=mymesh.cpp
mymesh_LDADD=$(LIFE_LIBS)

myapp_SOURCES=myapp.cpp
myapp_LDADD=$(LIFE_LIBS)

SUFFIXES=.cpp .hpp .tex .pdf

CLEANFILES = *.aux *.log *.lot *.lof *.toc *.bbl *.blg *.out *.idx *.ind *.ilg $(doc_DATA) version.tex

.tex.pdf:
	echo "*** INFO: first pdflatex compilation" ; \
	TEXINPUTS=:.//:$(top_builddir)/tutorial//:$(srcdir)//: $(PDFLATEX) $< ; \
	echo "*** INFO: index generation" ; \
	if test -f $*.idx; then $(MAKEINDEX) -s $(srcdir)/tutorial.ist $*.idx ; fi ; \
	echo "*** INFO: bibliography creation" ; \
	if test -f $(srcdir)/*.bib; then BIBINPUTS=:$(srcdir): $(BIBTEX) $*.aux ; fi ; \
	echo "*** INFO: table of content insertion" ; \
	TEXINPUTS_pdflatex=:.//:$(top_builddir)/tutorial//:$(srcdir)//: $(PDFLATEX) $< ; \
	echo "*** INFO: labels insertion" ; \
	TEXINPUTS_pdflatex=:.//:$(top_builddir)/tutorial//:$(srcdir)//: $(PDFLATEX) $< -o $(top_builddir)/$*.pdf



EXTRA_DIST = tutorial.tex version.tex.in feel.tex gfdl.tex tutorial.ist version.tex.in


tutorial.pdf : $(EXTRA_DIST) myapp.cpp mymesh.cpp myintegrals.cpp laplacian.cpp stokes.cpp bratu.cpp

pdf-local : tutorial.pdf
doc_DATA = tutorial.pdf
