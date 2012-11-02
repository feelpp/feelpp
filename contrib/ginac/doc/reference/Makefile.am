## Process this file with automake to produce Makefile.in

# Targets
if CONFIG_DOXYGEN
html: html_files/index.html

if CONFIG_TEX
dvi: reference.dvi
ps: reference.ps
pdf: reference.pdf
endif


# This section produces HTML'ed and TeX'ed developer's reference from the
# sources with a JavaDoc-like tool (in this case doxygen):
html_files::
	@[ -d html_files ] || mkdir html_files >/dev/null 2>&1

html_files/index.html:
	  @echo "Running ${DOXYGEN} DoxyfileHTML..."; \
	  ${DOXYGEN} DoxyfileHTML

if CONFIG_TEX
latex::
	@[ -d latex ] || mkdir latex >/dev/null 2>&1

latex/reference.tex:
	@echo "Running ${DOXYGEN} DoxyfileTEX..."; \
	  ${DOXYGEN} DoxyfileTEX; \
	  mv latex/refman.tex latex/reference.tex

latex/reference.dvi: latex/reference.tex
	@set -e ;\
	  cd latex; \
	  ${LATEX} reference.tex ; \
	  ${MAKEINDEX} reference.idx ; \
	  ${LATEX} reference.tex

latex/reference.ps: latex/reference.dvi
	@echo "Running ${DVIPS} -o reference.ps reference.dvi..."; \
	  cd latex; ${DVIPS} -o reference.ps reference.dvi

pdflatex::
	@[ -d pdflatex ] || mkdir pdflatex >/dev/null 2>&1

pdflatex/reference.tex:
	@echo "Running ${DOXYGEN} DoxyfilePDF..."; \
	  ${DOXYGEN} DoxyfilePDF; \
	  mv pdflatex/refman.tex pdflatex/reference.tex

pdflatex/reference.pdf: pdflatex/reference.tex
	@set -e ; \
		cd pdflatex; \
	  ${PDFLATEX} reference.tex ;\
	  ${MAKEINDEX} reference.idx ;\
	  ${PYTHON} $(top_srcdir)/scripts/fixupind.py reference.ind; \
	  ${PDFLATEX} reference.tex

reference.dvi: latex latex/reference.dvi
	cp latex/reference.dvi reference.dvi

reference.ps: latex latex/reference.ps
	cp latex/reference.ps reference.ps

reference.pdf: pdflatex pdflatex/reference.pdf
	cp pdflatex/reference.pdf reference.pdf
endif 
## CONFIG_TEX

endif
## CONFIG_DOXYGEN

CLEANFILES = reference.dvi reference.ps reference.pdf html_files/* latex/* pdflatex/*
