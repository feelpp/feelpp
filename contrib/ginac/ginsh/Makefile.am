## Process this file with automake to produce Makefile.in

bin_PROGRAMS = ginsh
BUILT_SOURCES = ginsh_parser.cpp ginsh_lexer.cpp ginsh_parser.hpp
ginsh_SOURCES = ginsh_parser.cpp ginsh_parser.hpp ginsh_lexer.lpp ginsh.h ginsh_extensions.h
ginsh_LDADD = ../ginac/libginac.la $(GINSH_LIBS)

AM_CPPFLAGS = -I$(srcdir)/../ginac -I../ginac -DIN_GINAC

man_MANS = ginsh.1

CLEANFILES = ginsh_fcn_help.h ginsh_op_help.h
EXTRA_DIST = ginsh_parser.ypp ginsh_fcn_help.py ginsh_op_help.py CMakeLists.txt

# files created by sed scripts
ginsh_fcn_help.h: ginsh.1.in $(srcdir)/ginsh_fcn_help.py
	$(PYTHON) $(srcdir)/ginsh_fcn_help.py -o $@ $<

ginsh_op_help.h: ginsh.1 $(srcdir)/ginsh_op_help.py
	$(PYTHON) $(srcdir)/ginsh_op_help.py -o $@ $<

# automake <=1.11 and automake >=1.12 have different conventions for naming C++ header files
# made by yacc.  To work with both, we write our own rule rather than using automake's.
ginsh_parser.cpp ginsh_parser.hpp: ginsh_parser.ypp
	$(YACC) $(YFLAGS) -o ginsh_parser.cpp --defines=ginsh_parser.hpp $<

ginsh_parser.o: ginsh_fcn_help.h ginsh_op_help.h ginsh_parser.hpp
