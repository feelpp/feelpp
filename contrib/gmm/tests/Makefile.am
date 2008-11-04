
check_PROGRAMS = dummy


dummy_SOURCES = dummy.C

INCLUDES = -I$(top_srcdir)/include -I../include

LDADD    = -lm @SUPLDFLAGS@


EXTRA_DIST=\
          make_gmm_test.pl\
          dummy.C\
          gmm_torture01_lusolve.C\
          gmm_torture02_baseop.C\
          gmm_torture05_mult.C\
          gmm_torture06_mat_mult.C\
          gmm_torture10_qr.C\
          gmm_torture15_sub.C\
          gmm_torture20_iterative_solvers.C

#          auto_gmm_test_mat_mult.C

CLEANFILES = toto.mat ii_files/* auto_gmm*

# C.P. Commented this
# TESTS = $(top_srcdir)/tests/make_gmm_test.pl
# TESTS_ENVIRONMENT = perl

