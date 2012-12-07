/* 
 * Copyright (C) 1999, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2011
 * Free Software Foundation, Inc.
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

#ifndef MATHEVAL_H
#define MATHEVAL_H 1

#ifdef __cplusplus
extern          "C" {
#endif

	/* Create evaluator from string representing function.  Function
	 * returns pointer that should be passed as first argument to all
	 * other library functions.  If an error occurs, function will
	 * return null pointer. */
	extern void    *evaluator_create(char *string);

	/* Destroy evaluator specified. */
	extern void     evaluator_destroy(void *evaluator);

	/* Evaluate function represented by evaluator given.  Variable
	 * names and respective values are represented by function third
	 * and fourth argument. Number of variables i.e. length of these
	 * two arrays is given by second argument.  Function returns
	 * evaluated function value.  In case that function contains
	 * variables with names not given through third function argument, 
	 * value of this variable is undeterminated. */
	extern double   evaluator_evaluate(void *evaluator, int count,
					   char **names, double *values);

	/* Return textual representation of function given by evaluator.
	 * Textual representation is built after evaluator simplification, 
	 * so it may differ from original string supplied when creating
	 * evaluator.  String representing function is allocated,
	 * remembered and later destroyed by evaluator object, thus caller 
	 * must not free returned pointer.  Returned information is valid
	 * until evaluator object destroyed. */
	extern char    *evaluator_get_string(void *evaluator);

	/* Get array of strings with names of variables appearing in
	 * function represented by given evaluator.  Only variables
	 * referenced by evaluator after simplification are returned.
	 * Address of first string in array is stored into location
	 * pointed by function second argument.  Number of array elements
	 * is stored into location pointed by third argument.  Array is
	 * allocated, remembered and later destroyed by evaluator object,
	 * thus caller must not free any of string nor array itself.
	 * Returned information is valid until evaluator object destroyed. 
	 */
	extern void     evaluator_get_variables(void *evaluator,
						char ***names, int *count);

	/* Create evaluator for first derivative of function represented
	 * by evaluator given as first argument using derivative variable
	 * given as second argument. */
	extern void    *evaluator_derivative(void *evaluator, char *name);

	/* Helper functions to simplify evaluation when variable names are 
	 * "x", "x" and "y" or "x" and "y" and "z" respectively. */
	extern double   evaluator_evaluate_x(void *evaluator, double x);
	extern double   evaluator_evaluate_x_y(void *evaluator, double x,
					       double y);
	extern double   evaluator_evaluate_x_y_z(void *evaluator, double x,
						 double y, double z);

	/* Helper functions to simplify differentiation when variable
	 * names are "x" or "y" or "z" respectively. */
	extern void    *evaluator_derivative_x(void *evaluator);
	extern void    *evaluator_derivative_y(void *evaluator);
	extern void    *evaluator_derivative_z(void *evaluator);

#ifdef __cplusplus
}
#endif
#endif
