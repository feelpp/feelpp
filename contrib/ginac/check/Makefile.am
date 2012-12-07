## Process this file with automake to produce Makefile.in

CHECKS = check_numeric \
	 check_inifcns \
	 check_matrices \
	 check_lsolve 

EXAMS = exam_paranoia \
	exam_heur_gcd \
	match_bug \
	parser_bugs \
	exam_numeric_archive \
	exam_numeric \
	exam_powerlaws  \
	exam_inifcns \
	exam_inifcns_nstdsums \
	exam_differentiation  \
	exam_polygcd  \
	exam_normalization  \
	exam_factor \
	exam_pseries  \
	exam_matrices  \
	exam_lsolve  \
	exam_indexed  \
	exam_color  \
	exam_clifford  \
	exam_archive  \
	exam_structure  \
	exam_hashmap  \
	exam_misc \
	exam_mod_gcd \
	bugme_chinrem_gcd \
	pgcd_relatively_prime_bug \
	pgcd_infinite_loop \
	exam_cra

TIMES = time_dennyfliegner \
	time_gammaseries \
	time_vandermonde \
	time_toeplitz \
	time_hashmap \
	time_lw_A \
	time_lw_B \
	time_lw_C \
	time_lw_D \
	time_lw_E \
	time_lw_F \
	time_lw_G \
	time_lw_H \
	time_lw_IJKL \
	time_lw_M1 \
	time_lw_M2 \
	time_lw_N \
	time_lw_O \
	time_lw_P \
	time_lw_Pprime \
	time_lw_Q \
	time_lw_Qprime \
	time_antipode \
	time_fateman_expand \
	time_uvar_gcd \
	time_parser

TESTS = $(CHECKS) $(EXAMS) $(TIMES)
check_PROGRAMS = $(CHECKS) $(EXAMS) $(TIMES)

check_numeric_SOURCES = check_numeric.cpp 
check_numeric_LDADD = ../ginac/libginac.la

check_inifcns_SOURCES = check_inifcns.cpp
check_inifcns_LDADD = ../ginac/libginac.la

check_matrices_SOURCES = check_matrices.cpp genex.cpp
check_matrices_LDADD = ../ginac/libginac.la

check_lsolve_SOURCES = check_lsolve.cpp genex.cpp
check_lsolve_LDADD = ../ginac/libginac.la

exam_paranoia_SOURCES = exam_paranoia.cpp
exam_paranoia_LDADD = ../ginac/libginac.la

exam_heur_gcd_SOURCES = heur_gcd_bug.cpp 
exam_heur_gcd_LDADD = ../ginac/libginac.la

match_bug_SOURCES = match_bug.cpp error_report.h
match_bug_LDADD = ../ginac/libginac.la

parser_bugs_SOURCES = parser_bugs.cpp
parser_bugs_LDADD = ../ginac/libginac.la

exam_numeric_archive_SOURCES = numeric_archive.cpp
exam_numeric_archive_LDADD = ../ginac/libginac.la

exam_numeric_SOURCES = exam_numeric.cpp
exam_numeric_LDADD = ../ginac/libginac.la

exam_powerlaws_SOURCES = exam_powerlaws.cpp
exam_powerlaws_LDADD = ../ginac/libginac.la

exam_inifcns_SOURCES = exam_inifcns.cpp
exam_inifcns_LDADD = ../ginac/libginac.la

exam_inifcns_nstdsums_SOURCES = exam_inifcns_nstdsums.cpp \
				exam_inifcns_nstdsums.h
exam_inifcns_nstdsums_LDADD = ../ginac/libginac.la

exam_differentiation_SOURCES = exam_differentiation.cpp
exam_differentiation_LDADD = ../ginac/libginac.la

exam_polygcd_SOURCES = exam_polygcd.cpp
exam_polygcd_LDADD = ../ginac/libginac.la

exam_normalization_SOURCES = exam_normalization.cpp
exam_normalization_LDADD = ../ginac/libginac.la

exam_factor_SOURCES = exam_factor.cpp
exam_factor_LDADD = ../ginac/libginac.la

exam_pseries_SOURCES = exam_pseries.cpp
exam_pseries_LDADD = ../ginac/libginac.la

exam_matrices_SOURCES = exam_matrices.cpp
exam_matrices_LDADD = ../ginac/libginac.la

exam_lsolve_SOURCES = exam_lsolve.cpp
exam_lsolve_LDADD = ../ginac/libginac.la

exam_indexed_SOURCES = exam_indexed.cpp
exam_indexed_LDADD = ../ginac/libginac.la

exam_color_SOURCES = exam_color.cpp
exam_color_LDADD = ../ginac/libginac.la

exam_clifford_SOURCES = exam_clifford.cpp
exam_clifford_LDADD = ../ginac/libginac.la

exam_archive_SOURCES = exam_archive.cpp
exam_archive_LDADD = ../ginac/libginac.la

exam_structure_SOURCES = exam_structure.cpp
exam_structure_LDADD = ../ginac/libginac.la

exam_hashmap_SOURCES = exam_hashmap.cpp
exam_hashmap_LDADD = ../ginac/libginac.la

exam_misc_SOURCES = exam_misc.cpp
exam_misc_LDADD = ../ginac/libginac.la

exam_mod_gcd_SOURCES = exam_mod_gcd.cpp
exam_mod_gcd_LDADD = ../ginac/libginac.la

exam_cra_SOURCES = exam_cra.cpp
exam_cra_LDADD = ../ginac/libginac.la

time_dennyfliegner_SOURCES = time_dennyfliegner.cpp \
			     randomize_serials.cpp timer.cpp timer.h
time_dennyfliegner_LDADD = ../ginac/libginac.la

time_gammaseries_SOURCES = time_gammaseries.cpp \
			   randomize_serials.cpp timer.cpp timer.h
time_gammaseries_LDADD = ../ginac/libginac.la

time_vandermonde_SOURCES = time_vandermonde.cpp \
			   randomize_serials.cpp timer.cpp timer.h
time_vandermonde_LDADD = ../ginac/libginac.la

time_toeplitz_SOURCES = time_toeplitz.cpp \
			randomize_serials.cpp timer.cpp timer.h
time_toeplitz_LDADD = ../ginac/libginac.la
time_hashmap_SOURCES = time_hashmap.cpp \
		       randomize_serials.cpp timer.cpp timer.h
time_hashmap_LDADD = ../ginac/libginac.la

time_lw_A_SOURCES = time_lw_A.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_A_LDADD = ../ginac/libginac.la

time_lw_B_SOURCES = time_lw_B.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_B_LDADD = ../ginac/libginac.la

time_lw_C_SOURCES = time_lw_C.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_C_LDADD = ../ginac/libginac.la

time_lw_D_SOURCES = time_lw_D.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_D_LDADD = ../ginac/libginac.la

time_lw_E_SOURCES = time_lw_E.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_E_LDADD = ../ginac/libginac.la

time_lw_F_SOURCES = time_lw_F.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_F_LDADD = ../ginac/libginac.la

time_lw_G_SOURCES = time_lw_G.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_G_LDADD = ../ginac/libginac.la

time_lw_H_SOURCES = time_lw_H.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_H_LDADD = ../ginac/libginac.la

time_lw_IJKL_SOURCES = time_lw_IJKL.cpp \
		       randomize_serials.cpp timer.cpp timer.h
time_lw_IJKL_LDADD = ../ginac/libginac.la

time_lw_M1_SOURCES = time_lw_M1.cpp \
		     randomize_serials.cpp timer.cpp timer.h
time_lw_M1_LDADD = ../ginac/libginac.la

time_lw_M2_SOURCES = time_lw_M2.cpp \
		     randomize_serials.cpp timer.cpp timer.h
time_lw_M2_LDADD = ../ginac/libginac.la

time_lw_N_SOURCES = time_lw_N.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_N_LDADD = ../ginac/libginac.la

time_lw_O_SOURCES = time_lw_O.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_O_LDADD = ../ginac/libginac.la

time_lw_P_SOURCES = time_lw_P.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_P_LDADD = ../ginac/libginac.la

time_lw_Pprime_SOURCES = time_lw_Pprime.cpp \
			 randomize_serials.cpp timer.cpp timer.h
time_lw_Pprime_LDADD = ../ginac/libginac.la

time_lw_Q_SOURCES = time_lw_Q.cpp \
		    randomize_serials.cpp timer.cpp timer.h
time_lw_Q_LDADD = ../ginac/libginac.la

time_lw_Qprime_SOURCES = time_lw_Qprime.cpp \
			 randomize_serials.cpp timer.cpp timer.h
time_lw_Qprime_LDADD = ../ginac/libginac.la

time_antipode_SOURCES = time_antipode.cpp \
			randomize_serials.cpp timer.cpp timer.h
time_antipode_LDADD = ../ginac/libginac.la

time_fateman_expand_SOURCES = time_fateman_expand.cpp \
			      randomize_serials.cpp timer.cpp timer.h
time_fateman_expand_LDADD = ../ginac/libginac.la

time_uvar_gcd_SOURCES = time_uvar_gcd.cpp test_runner.h timer.cpp timer.h
time_uvar_gcd_LDADD = ../ginac/libginac.la

time_parser_SOURCES = time_parser.cpp \
		      randomize_serials.cpp timer.cpp timer.h
time_parser_LDADD = ../ginac/libginac.la

bugme_chinrem_gcd_SOURCES = bugme_chinrem_gcd.cpp
bugme_chinrem_gcd_LDADD = ../ginac/libginac.la

pgcd_relatively_prime_bug_SOURCES = pgcd_relatively_prime_bug.cpp
pgcd_relatively_prime_bug_LDADD = ../ginac/libginac.la

pgcd_infinite_loop_SOURCES =  pgcd_infinite_loop.cpp
pgcd_infinite_loop_LDADD =  ../ginac/libginac.la

AM_CPPFLAGS = -I$(srcdir)/../ginac -I../ginac -DIN_GINAC

CLEANFILES = exam.gar
