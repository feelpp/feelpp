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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdarg.h>
#include "common.h"
#include "symbol_table.h"
#include "xmath.h"

/* Type definition for function accepting single argument of double type
 * and returning double value. */
typedef double  (*function_type) (double);

/* Calculate hash value for given name and hash table length.  */
static int      hash(char *name, int length);

SymbolTable    *
symbol_table_create(int length)
{
	SymbolTable    *symbol_table;	/* Pointer to symbol table.  */
	static char    *constants_names[] = {
		"e", "log2e", "log10e", "ln2", "ln10", "pi", "pi_2",
		"pi_4", "1_pi", "2_pi", "2_sqrtpi", "sqrt2", "sqrt1_2"
	};			/* Symbol table predefined constants
				 * names. */
	static double   constants[] = {
		2.7182818284590452354, 1.4426950408889634074,
		0.43429448190325182765, 0.69314718055994530942,
		2.30258509299404568402, 3.14159265358979323846,
		1.57079632679489661923, 0.78539816339744830962,
		0.31830988618379067154, 0.63661977236758134308,
		1.12837916709551257390, 1.41421356237309504880,
		0.70710678118654752440
	};			/* Symbol table predefined constants
				 * values. */
	static char    *functions_names[] = {
		"exp", "log", "sqrt", "sin", "cos", "tan", "cot", "sec",
		"csc", "asin", "acos", "atan", "acot", "asec", "acsc",
		"sinh", "cosh", "tanh", "coth", "sech", "csch",
		"asinh", "acosh", "atanh", "acoth", "asech", "acsch",
		"abs", "step", "delta", "nandelta", "erf", 
	};			/* Symbol table predefined functions
				 * names. */
	static double   (*functions[]) (double) = {
	exp, log, sqrt, sin, cos, tan, math_cot, math_sec, math_csc, asin, acos, atan, math_acot, math_asec, math_acsc, sinh, cosh, tanh, math_coth, math_sech, math_csch, math_asinh, math_acosh, math_atanh, math_acoth, math_asech, math_acsch, fabs, math_step, math_delta, math_nandelta, erf};	/* Symbol 
																																				 * table 
																																				 * predefined 
																																				 * functions 
																																				 * pointers 
																																				 * to 
																																				 * functions 
																																				 * to 
																																				 * calculate 
																																				 * them. 
																																				 */
	int             i;	/* Loop counter.  */

	/* Allocate memory for symbol table data structure as well as for
	 * corresponding hash table. */
	symbol_table = XMALLOC(SymbolTable, 1);
	symbol_table->length = length;
	symbol_table->records = XCALLOC(Record, symbol_table->length);

	/* Insert predefined constants into symbol table. */
	for (i = 0;
	     i < sizeof(constants_names) / sizeof(constants_names[0]); i++)
		symbol_table_insert(symbol_table, constants_names[i], 'c',
				    constants[i]);

	/* Insert predefined functions into symbol table. */
	for (i = 0;
	     i < sizeof(functions_names) / sizeof(functions_names[0]); i++)
		symbol_table_insert(symbol_table, functions_names[i], 'f',
				    functions[i]);

	/* Initialize symbol table reference count. */
	symbol_table->reference_count = 1;

	return symbol_table;
}

void
symbol_table_destroy(SymbolTable * symbol_table)
{
	Record         *curr,
	               *next;	/* Pointers to current and next record
				 * while traversing hash table bucket.  */
	int             i;	/* Loop counter.  */

	/* Decrement refernce count and return if symbol table still used
	 * elsewhere. */
	if (--symbol_table->reference_count > 0)
		return;

	/* Delete hash table as well as data structure representing symbol 
	 * table. */
	for (i = 0; i < symbol_table->length; i++)
		for (curr = symbol_table->records[i].next; curr;) {
			next = curr->next;
			XFREE(curr->name);
			XFREE(curr);
			curr = next;
		}
	XFREE(symbol_table->records);
	XFREE(symbol_table);
}

Record         *
symbol_table_insert(SymbolTable * symbol_table, char *name, char type, ...)
{
	Record         *record;	/* Pointer to symbol table record
				 * corresponding to name given.  */
	va_list         ap;	/* Function variable argument list.  */
	int             i;	/* Loop counter.  */

	/* Check if symbol already in table and, if affirmative and record 
	 * type same as type given, return corresponding record
	 * immediately. */
	if ((record = symbol_table_lookup(symbol_table, name))) {
		assert(record->type == type);
		return record;
	}

	/* Allocate memory for and initialize new record. */
	record = XMALLOC(Record, 1);
	record->name = XMALLOC(char, strlen(name) + 1);
	strcpy(record->name, name);
	record->type = type;
	record->flag = FALSE;

	/* Parse function variable argument list to complete record
	 * initialization. */
	va_start(ap, type);
	switch (record->type) {
	case 'c':
		record->data.value = va_arg(ap, double);
		break;

	case 'v':
		record->data.value = 0;
		break;

	case 'f':
		record->data.function = va_arg(ap, function_type);
		break;
	}
	va_end(ap);

	/* Calculate hash value and put record in corresponding hash table 
	 * bucket. */
	i = hash(name, symbol_table->length);
	record->next = symbol_table->records[i].next;
	symbol_table->records[i].next = record;

	return record;
}

Record         *
symbol_table_lookup(SymbolTable * symbol_table, char *name)
{
	int             i;	/* Hash value. */
	Record         *curr;	/* Pointer to current symbol table record. 
				 */

	/* 
	 * Calcuate hash value for name given.
	 */
	i = hash(name, symbol_table->length);

	/* Lookup for name in hash table bucket corresponding to above
	 * hash value. */
	for (curr = symbol_table->records[i].next; curr; curr = curr->next)
		if (!strcmp(curr->name, name))
			return curr;

	return NULL;
}

void
symbol_table_clear_flags(SymbolTable * symbol_table)
{
	Record         *curr;	/* Pointer to current symbol table record
				 * while traversing hash table bucket. */
	int             i;	/* Loop counter.  */

	/* Clear flag for all records in symbol table. */
	for (i = 0; i < symbol_table->length; i++)
		for (curr = symbol_table->records[i].next; curr;
		     curr = curr->next)
			curr->flag = FALSE;
}

int
symbol_table_get_flagged_count(SymbolTable * symbol_table)
{
	int             count;	/* Number of flagged symbol table records. 
				 */
	Record         *curr;	/* Pointer to current symbol table record
				 * while traversing hash table bucket. */
	int             i;	/* Loop counter.  */

	/* Calculate number of records in symbol table with flag set. */
	count = 0;
	for (i = 0; i < symbol_table->length; i++)
		for (curr = symbol_table->records[i].next; curr;
		     curr = curr->next)
			if (curr->flag)
				count++;
	return count;
}

int
symbol_table_get_flagged(SymbolTable * symbol_table, Record ** records,
			 int length)
{
	int             count;	/* Number of pointers to symbol table
				 * records put into given array. */
	Record         *curr;	/* Pointers to current symbol table record 
				 * while traversing hash table bucket.  */
	int             i;	/* Loop counter.  */

	/* Put pointers to records in symbol table with flag set into
	 * given array. */
	count = 0;
	for (i = 0; i < symbol_table->length; i++)
		for (curr = symbol_table->records[i].next; curr;
		     curr = curr->next)
			if (curr->flag) {
				records[count++] = curr;
				if (count == length)
					return count;
			}
	return count;
}

SymbolTable    *
symbol_table_assign(SymbolTable * symbol_table)
{
	/* Increase symbol table reference count and return pointer to
	 * data structure representing table. */
	symbol_table->reference_count++;
	return symbol_table;
}

/* Function below reused from A.V. Aho, R. Sethi, J.D. Ullman, "Compilers
 * - Principle, Techniques, and Tools", Addison-Wesley, 1986, pp 435-437,
 * and in turn from P.J. Weineberger's C compiler. */
static int
hash(char *s, int n)
{
	char           *p;
	unsigned        h,
	                g;

	h = 0;

	for (p = s; *p; p++) {
		h = (h << 4) + *p;
		if ((g = h & 0xf0000000)) {
			h = h ^ (g >> 24);
			h = h ^ g;
		}
	}

	return h % n;
}
