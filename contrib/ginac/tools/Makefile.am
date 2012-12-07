## Process this file with automake to produce Makefile.in

bin_PROGRAMS = viewgar
viewgar_SOURCES = viewgar.cpp
viewgar_LDADD = ../ginac/libginac.la

bin_SCRIPTS = 
if CONFIG_EXCOMPILER
bin_SCRIPTS += ginac-excompiler
endif

AM_CPPFLAGS = -I$(srcdir)/../ginac -I../ginac -DIN_GINAC

man_MANS = viewgar.1

EXTRA_DIST = ginac-excompiler.in
