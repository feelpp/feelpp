## Process this file with automake to produce Makefile.in

EXFILES = archive1.cpp compile1.cpp compile2.cpp compile3.cpp lanczos.cpp \
 mystring.cpp derivative.cpp

TEXI = ginac-examples.texi

info_TEXINFOS = $(TEXI)

TXT = ginac-examples.txt
HTML = ginac-examples.html

all: $(EXFILES)
	${MAKEINFO} --no-split --no-headers $(srcdir)/$(TEXI) > $(TXT)
	${MAKEINFO} --no-split --html $(srcdir)/$(TEXI)

CLEANFILES = $(TXT) $(HTML)

EXTRA_DIST = $(EXFILES)
