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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_STDINT_H
#include <stdint.h>
#else
#error no <stdint.h> avaialable
#endif

#include "common.h"
#include "matheval.h"

/* Wrapper for evaluator_create() function.  Evaluator objects will be
 * passed between Fortran and C as 64-bit integers instead of void
 * pointers, thus size of void* is assumed to be less or equal of 64 bits
 * throughout this interface. */
int64_t
evaluator_create__(char *string, int length)
{
	char           *stringz;	/* Zero-terminated string
					 * representing function.  */
	int64_t         evaluator;	/* Evaluator created for function. 
					 */

	/* Copy string passed from Fortran code and terminate it with
	 * zero. */
	stringz = XMALLOC(char, length + 1);
	memcpy(stringz, string, length * sizeof(char));
	stringz[length] = '\0';

	/* Call evaluator_create() function. */
	evaluator = (int64_t) evaluator_create(stringz);

	/* Free string used to create evaluator. */
	XFREE(stringz);

	return evaluator;
}

/* Wrapper for evaluator_destroy() function.  */
void
evaluator_destroy__(int64_t * evaluator)
{
	evaluator_destroy((void *) *evaluator);
}

/* Wrapper for evaluator_evaluate() function.  */
double
evaluator_evaluate__(int64_t * evaluator, int *count, char *names,
		     double *values, int length)
{
	char          **names_copy;	/* Copy of variable names.  Names
					 * are passed in single string
					 * from Fortran code, delimited by 
					 * blanks, while
					 * evaluator_evaluate() function
					 * expects array of strings.  */
	double          result;	/* Calculated value of function.  */
	int             i,
	                j,
	                n;	/* Loop counters.  */

	/* Parse string containing variable names and create array of
	 * strings with each string containing single name. */
	names_copy = XMALLOC(char *, *count);
	for (i = j = 0; i < *count && j < length; i++, j += n) {
		for (; names[j] == ' '; j++);
		for (n = 1; j + n < length && !(names[j + n] == ' '); n++);
		names_copy[i] = XMALLOC(char, n + 1);
		memcpy(names_copy[i], names + j, n * sizeof(char));
		names_copy[i][n] = '\0';
	}

	/* Call evaluator_evaluate() function. */
	result =
	    evaluator_evaluate((void *) *evaluator, *count, names_copy,
			       values);

	/* Free memory used. */
	for (i = 0; i < *count; i++)
		XFREE(names_copy[i]);
	XFREE(names_copy);

	return result;
}

/* First in pair of wrappers for evaluator_get_string() function.  */
int
evaluator_get_string_length__(int64_t * evaluator)
{
	/* Return length of evaluator textual respresentation. */
	return strlen(evaluator_get_string((void *) *evaluator));
}

/* Second in pair of wrappers for evaluator_get_string() function.  */
void
evaluator_get_string_chars__(int64_t * evaluator, char *string, int length)
{
	/* Copy evaluator textual respresentation to string passed from
	 * Fortran code. */
	memcpy(string, evaluator_get_string((void *) *evaluator),
	       length * sizeof(char));
}

/* First in pair of wrappers for evaluator_get_variables() function.  */
int
evaluator_get_variables_length__(int64_t * evaluator)
{
	char          **names;	/* Array with variable names. */
	int             count;	/* Number of elements in above array. */
	int             length;	/* Length of string with concatenated
				 * variable names */
	int             i;	/* Loop counter. */

	/* Get array of strings with variable names. */
	evaluator_get_variables((void *) *evaluator, &names, &count);

	/* Calculate length of string with concatenated names from above
	 * array. */
	length = 0;
	for (i = 0; i < count; i++) {
		if (i != 0)
			length++;
		length += strlen(names[i]);
	}

	return length;
}

/* Second in pair of wrappers for evaluator_get_variables() function.  */
void
evaluator_get_variables_chars__(int64_t * evaluator, char *string,
				int length)
{
	char          **names;	/* Array with variable names. */
	int             count;	/* Number of elements in above array. */
	int             n;	/* Number of characters to be copied from
				 * current variable name to string with
				 * concatenated variable names. */
	int             i;	/* Loop counter. */

	/* Get array of strings with variable names. */
	evaluator_get_variables((void *) *evaluator, &names, &count);

	/* Concatenate variable names from above array into string passed
	 * from Fortran code. */
	for (i = 0; i < count; i++) {
		if (i != 0 && length > 0) {
			*(string++) = ' ';
			length--;
		}
		n = strlen(names[i]);
		if (n > length)
			n = length;
		memcpy(string, names[i], n * sizeof(char));
		string += n;
		length -= n;
	}
}

/* Wrapper for evaluator_derivative() function.  */
int64_t
evaluator_derivative__(int64_t * evaluator, char *name, int length)
{
	char           *stringz;	/* Zero terminated string
					 * containing derivation variable
					 * name.  */
	int64_t         derivative;	/* Evaluator for function
					 * derivative.  */

	/* Copy variable name passed from Fortran code and terminate it
	 * with zero. */
	stringz = XMALLOC(char, length + 1);
	memcpy(stringz, name, length * sizeof(char));
	stringz[length] = '\0';

	/* Call evaluator_derivative() function. */
	derivative =
	    (int64_t) evaluator_derivative((void *) *evaluator, stringz);

	/* Free string containing derivation variable name. */
	XFREE(stringz);

	return derivative;
}

/* Wrapper for evaluator_evaluate_x() function.  */
double
evaluator_evaluate_x__(int64_t * evaluator, double *x)
{
	return evaluator_evaluate_x((void *) *evaluator, *x);
}

/* Wrapper for evaluator_evaluate_x_y() function.  */
double
evaluator_evaluate_x_y__(int64_t * evaluator, double *x, double *y)
{
	return evaluator_evaluate_x_y((void *) *evaluator, *x, *y);
}

/* Wrapper for evaluator_evaluate_x_y_z() function.  */
double
evaluator_evaluate_x_y_z__(int64_t * evaluator, double *x, double *y,
			   double *z)
{
	return evaluator_evaluate_x_y_z((void *) *evaluator, *x, *y, *z);
}

/* Wrapper for evaluator_derivative_x() function.  */
int64_t
evaluator_derivative_x__(int64_t * evaluator)
{
	return (int64_t) evaluator_derivative_x((void *) *evaluator);
}

/* Wrapper for evaluator_derivative_y() function.  */
int64_t
evaluator_derivative_y__(int64_t * evaluator)
{
	return (int64_t) evaluator_derivative_y((void *) *evaluator);
}

/* Wrapper for evaluator_derivative_z() function.  */
int64_t
evaluator_derivative_z__(int64_t * evaluator)
{
	return (int64_t) evaluator_derivative_z((void *) *evaluator);
}
