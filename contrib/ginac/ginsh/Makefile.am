## Process this file with automake to produce Makefile.in

bin_PROGRAMS = ginsh
ginsh_SOURCES = ginsh_parser.yy ginsh_lexer.ll ginsh.h ginsh_extensions.h
ginsh_LDADD = ../ginac/libginac.la $(GINSH_LIBS)

INCLUDES = -I$(srcdir)/../ginac -I../ginac -DIN_GINAC
AM_YFLAGS = -d

man_MANS = ginsh.1

CLEANFILES = ginsh_fcn_help.h ginsh_op_help.h
EXTRA_DIST = ginsh_parser.h ginsh_fcn_help.py ginsh_op_help.py

# files created by sed scripts
ginsh_fcn_help.h: ginsh.1.in $(srcdir)/ginsh_fcn_help.py
	$(PYTHON) $(srcdir)/ginsh_fcn_help.py -o $@ $<

ginsh_op_help.h: ginsh.1 $(srcdir)/ginsh_op_help.py
	$(PYTHON) $(srcdir)/ginsh_op_help.py -o $@ $<

ginsh_parser.o: ginsh_fcn_help.h ginsh_op_help.h ginsh_parser.h
