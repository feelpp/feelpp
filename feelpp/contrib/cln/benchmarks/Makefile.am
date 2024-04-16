noinst_PROGRAMS = timebench1 timebench2a timebench2ap timebench2b

timebench1_SOURCES = timebench1.cc
timebench1_LDADD = ../src/libcln.la

timebench2a_SOURCES = timebench2a.cc
timebench2a_LDADD = ../src/libcln.la

timebench2ap_SOURCES = timebench2ap.cc
timebench2ap_LDADD = ../src/libcln.la

timebench2b_SOURCES = timebench2b.cc
timebench2b_LDADD = ../src/libcln.la

AM_CPPFLAGS = -I../include -I$(top_srcdir)/include
DEFAULT_INCLUDES = -I.@am__isrc@

EXTRA_DIST = README bench2.txt timebench2.sh


# Remove .deps directories created by 'configure'.
# I would say that this left-over is an Automake 1.16.1 bug.
distclean-local:
	rm -rf $(DEPDIR)
