/* 
 * Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2011 Free
 * Software Foundation, Inc.
 * 
 * This file is part of GNU libmatheval
 * 
 * GNU libmatheval is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * GNU libmatheval is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU libmatheval.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <libguile.h>
#include <matheval.h>
#include "config.h"

#ifndef HAVE_SCM_T_BITS
typedef long    scm_t_bits;
#endif

#ifndef HAVE_SCM_NUM2DBL
#ifdef SCM_NUM2DBL
#define scm_num2dbl(x,s) SCM_NUM2DBL(x)
#else
#error Neither scm_num2dbl() nor SCM_NUM2DBL available
#endif
#endif

#ifndef HAVE_SCM_C_DEFINE_GSUBR
#ifdef HAVE_SCM_MAKE_GSUBR
#define scm_c_define_gsubr scm_make_gsubr
#else
#error Neither scm_c_define_gsubr() nor scm_make_gsubr() available
#endif
#endif

static scm_t_bits evaluator_tag;	/* Unique identifier for Guile
					 * objects of evaluator type.  */

/* Guile interface for libmatheval library.  Procedures below are simple
 * wrappers for corresponding libmatheval procedures. */
static scm_sizet evaluator_destroy_scm(SCM evaluator_smob);
static SCM      evaluator_create_scm(SCM string);
static SCM      evaluator_evaluate_scm(SCM evaluator_smob, SCM count,
				       SCM names, SCM values);
static SCM      evaluator_get_string_scm(SCM evaluator_smob);
static SCM      evaluator_get_variables_scm(SCM evaluator_smob);
static SCM      evaluator_derivative_scm(SCM evaluator_smob, SCM name);
static SCM      evaluator_evaluate_x_scm(SCM evaluator_smob, SCM x);
static SCM      evaluator_evaluate_x_y_scm(SCM evaluator_smob, SCM x,
					   SCM y);
static SCM      evaluator_evaluate_x_y_z_scm(SCM evaluator_smob, SCM x,
					     SCM y, SCM z);
static SCM      evaluator_derivative_x_scm(SCM evaluator_smob);
static SCM      evaluator_derivative_y_scm(SCM evaluator_smob);
static SCM      evaluator_derivative_z_scm(SCM evaluator_smob);

static void
inner_main(void *closure, int argc, char **argv)
{
	/* Extend Guile with evaluator type and register procedure to free 
	 * objects of this type. */
	evaluator_tag = scm_make_smob_type("evaluator", sizeof(void *));
	scm_set_smob_free(evaluator_tag, evaluator_destroy_scm);

	/* Register other procedures working on evaluator type. */
	scm_c_define_gsubr("evaluator-create", 1, 0, 0,
			   (SCM(*)())evaluator_create_scm);
	scm_c_define_gsubr("evaluator-evaluate", 4, 0, 0,
			   (SCM(*)())evaluator_evaluate_scm);
	scm_c_define_gsubr("evaluator-get-string", 1, 0, 0,
			   (SCM(*)())evaluator_get_string_scm);
	scm_c_define_gsubr("evaluator-get-variables", 1, 0, 0,
			   (SCM(*)())evaluator_get_variables_scm);
	scm_c_define_gsubr("evaluator-derivative", 2, 0, 0,
			   (SCM(*)())evaluator_derivative_scm);
	scm_c_define_gsubr("evaluator-evaluate-x", 2, 0, 0,
			   (SCM(*)())evaluator_evaluate_x_scm);
	scm_c_define_gsubr("evaluator-evaluate-x-y", 3, 0, 0,
			   (SCM(*)())evaluator_evaluate_x_y_scm);
	scm_c_define_gsubr("evaluator-evaluate-x-y-z", 4, 0, 0,
			   (SCM(*)())evaluator_evaluate_x_y_z_scm);
	scm_c_define_gsubr("evaluator-derivative-x", 1, 0, 0,
			   (SCM(*)())evaluator_derivative_x_scm);
	scm_c_define_gsubr("evaluator-derivative-y", 1, 0, 0,
			   (SCM(*)())evaluator_derivative_y_scm);
	scm_c_define_gsubr("evaluator-derivative-z", 1, 0, 0,
			   (SCM(*)())evaluator_derivative_z_scm);

	/* Check is there exactly one argument left in command line. */
	assert(argc == 2);

	/* Interpret Guile code from file with name given through above
	 * argument.  */
	scm_primitive_load(scm_makfrom0str(argv[1]));
}

/* Program is demonstrating use of libmatheval library of procedures for
 * evaluating mathematical functions.  Program expects single argument
 * from command line and interpret Guile code (extended with procedures
 * from libmatheval Guile interface) from this file. */
int
main(int argc, char **argv)
{
	/* Initialize Guile library; inner_main() procedure gets called in 
	 * turn. */
	scm_boot_guile(argc, argv, inner_main, 0);

	exit(EXIT_SUCCESS);
}

/* Wrapper for evaluator_destroy() procedure from libmatheval library. */
static          scm_sizet
evaluator_destroy_scm(SCM evaluator_smob)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-destroy");

	evaluator_destroy((void *) SCM_CDR(evaluator_smob));

	return sizeof(void *);
}

/* Wrapper for evaluator_create() procedure from libmatheval library. */
static          SCM
evaluator_create_scm(SCM string)
{
	char           *stringz;
	void           *evaluator;

	SCM_ASSERT(SCM_NIMP(string)
		   && SCM_STRINGP(string), string, SCM_ARG1,
		   "evaluator-create");

	stringz = (char *) malloc((SCM_LENGTH(string) + 1) * sizeof(char));
	memcpy(stringz, SCM_CHARS(string), SCM_LENGTH(string));
	stringz[SCM_LENGTH(string)] = 0;

	evaluator = evaluator_create(stringz);

	free(stringz);

	SCM_RETURN_NEWSMOB(evaluator_tag, evaluator);
}

/* Wrapper for evaluator_evaluate() procedure from libmatheval library.
 * Variable names and values are passed as lists from Guile, so copies of
 * these argument have to be created in order to be able to call
 * evaluator_evaluate() procedure. */
static          SCM
evaluator_evaluate_scm(SCM evaluator_smob, SCM count, SCM names,
		       SCM values)
{
	SCM             name;
	char          **names_copy;
	SCM             value;
	double         *values_copy;
	double          result;
	int             i;

	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-evaluate");
	SCM_ASSERT(SCM_INUMP(count), count, SCM_ARG2,
		   "evaluator-evaluate");

	names_copy = (char **) malloc(SCM_INUM(count) * sizeof(char *));
	for (i = 0, name = names; i < SCM_INUM(count);
	     i++, name = SCM_CDR(name)) {
		SCM_ASSERT(SCM_NIMP(name) && SCM_CONSP(name)
			   && SCM_STRINGP(SCM_CAR(name)), names, SCM_ARG3,
			   "evaluator-evaluate");
		names_copy[i] =
		    (char *) malloc((SCM_LENGTH(SCM_CAR(name)) + 1) *
				    sizeof(char));
		memcpy(names_copy[i], SCM_CHARS(SCM_CAR(name)),
		       SCM_LENGTH(SCM_CAR(name)));
		names_copy[i][SCM_LENGTH(SCM_CAR(name))] = 0;
	}

	values_copy = (double *) malloc(SCM_INUM(count) * sizeof(double));
	for (i = 0, value = values; i < SCM_INUM(count);
	     i++, value = SCM_CDR(value)) {
		SCM_ASSERT(SCM_NIMP(value) && SCM_CONSP(value)
			   && SCM_NUMBERP(SCM_CAR(value)), values,
			   SCM_ARG4, "evaluator-evaluate");
		values_copy[i] =
		    scm_num2dbl(SCM_CAR(value), "evaluator-evaluate");
	}

	result =
	    evaluator_evaluate((void *) SCM_CDR(evaluator_smob),
			       SCM_INUM(count), names_copy, values_copy);

	for (i = 0; i < SCM_INUM(count); i++)
		free(names_copy[i]);
	free(names_copy);
	free(values_copy);

	return scm_make_real(result);
}

/* Wrapper for evaluator_get_string() procedure from libmatheval library. */
static          SCM
evaluator_get_string_scm(SCM evaluator_smob)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-get-string");

	return
	    scm_makfrom0str(evaluator_get_string
			    ((void *) SCM_CDR(evaluator_smob)));
}

/* Wrapper for evaluator_get_variables() procedure from libmatheval
 * library. */
static          SCM
evaluator_get_variables_scm(SCM evaluator_smob)
{
	char          **names;
	int             count;
	SCM             list;
	int             i;

	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-get-string");

	evaluator_get_variables((void *) SCM_CDR(evaluator_smob), &names,
				&count);
	list = SCM_EOL;
	for (i = 0; i < count; i++)
		list =
		    scm_append_x(scm_listify
				 (list,
				  scm_listify(scm_makfrom0str(names[i]),
					      SCM_UNDEFINED),
				  SCM_UNDEFINED));

	return list;
}

/* Wrapper for evaluator_derivative() procedure from libmatheval library. */
static          SCM
evaluator_derivative_scm(SCM evaluator_smob, SCM name)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-derivative");
	SCM_ASSERT(SCM_NIMP(name)
		   && SCM_STRINGP(name), name, SCM_ARG2,
		   "evaluator-derivative");
	SCM_RETURN_NEWSMOB(evaluator_tag,
			   evaluator_derivative((void *)
						SCM_CDR(evaluator_smob),
						SCM_CHARS(name)));
}

/* Wrapper for evaluator_evaluate_x() procedure from libmatheval library. */
static          SCM
evaluator_evaluate_x_scm(SCM evaluator_smob, SCM x)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-evaluate-x");
	SCM_ASSERT(SCM_NUMBERP(x), x, SCM_ARG2, "evaluator-evaluate-x");
	return
	    scm_make_real(evaluator_evaluate_x
			  ((void *) SCM_CDR(evaluator_smob),
			   scm_num2dbl(x, "evaluator-evaluate-x")));
}

/* Wrapper for evaluator_evaluate_x_y() procedure from libmatheval
 * library. */
static          SCM
evaluator_evaluate_x_y_scm(SCM evaluator_smob, SCM x, SCM y)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-evaluate-x-y");
	SCM_ASSERT(SCM_NUMBERP(x), x, SCM_ARG2, "evaluator-evaluate-x-y");
	SCM_ASSERT(SCM_NUMBERP(y), y, SCM_ARG3, "evaluator-evaluate-x-y");
	return
	    scm_make_real(evaluator_evaluate_x_y
			  ((void *) SCM_CDR(evaluator_smob),
			   scm_num2dbl(x, "evaluator-evaluate-x-y"),
			   scm_num2dbl(y, "evaluator-evaluate-x-y")));
}

/* Wrapper for evaluator_evaluate_x_y_z() procedure from libmatheval
 * library. */
static          SCM
evaluator_evaluate_x_y_z_scm(SCM evaluator_smob, SCM x, SCM y, SCM z)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-evaluate-x-y-z");
	SCM_ASSERT(SCM_NUMBERP(x), x, SCM_ARG2,
		   "evaluator-evaluate-x-y-z");
	SCM_ASSERT(SCM_NUMBERP(y), y, SCM_ARG3,
		   "evaluator-evaluate-x-y-z");
	SCM_ASSERT(SCM_NUMBERP(z), z, SCM_ARG4,
		   "evaluator-evaluate-x-y-z");
	return
	    scm_make_real(evaluator_evaluate_x_y_z
			  ((void *) SCM_CDR(evaluator_smob),
			   scm_num2dbl(x, "evaluator-evaluate-x-y-z"),
			   scm_num2dbl(y, "evaluator-evaluate-x-y-z"),
			   scm_num2dbl(z, "evaluator-evaluate-x-y-z")));
}

/* Wrapper for evaluator_derivative_x() procedure from libmatheval
 * library. */
static          SCM
evaluator_derivative_x_scm(SCM evaluator_smob)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-derivative-x");
	SCM_RETURN_NEWSMOB(evaluator_tag,
			   evaluator_derivative((void *)
						SCM_CDR(evaluator_smob),
						"x"));
}

/* Wrapper for evaluator_derivative_y() procedure from libmatheval
 * library. */
static          SCM
evaluator_derivative_y_scm(SCM evaluator_smob)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-derivative-y");
	SCM_RETURN_NEWSMOB(evaluator_tag,
			   evaluator_derivative((void *)
						SCM_CDR(evaluator_smob),
						"y"));
}

/* Wrapper for evaluator_derivative_z() procedure from libmatheval
 * library. */
static          SCM
evaluator_derivative_z_scm(SCM evaluator_smob)
{
	SCM_ASSERT((SCM_NIMP(evaluator_smob)
		    && SCM_SMOB_PREDICATE(evaluator_tag, evaluator_smob)),
		   evaluator_smob, SCM_ARG1, "evaluator-derivative-z");
	SCM_RETURN_NEWSMOB(evaluator_tag,
			   evaluator_derivative((void *)
						SCM_CDR(evaluator_smob),
						"z"));
}
