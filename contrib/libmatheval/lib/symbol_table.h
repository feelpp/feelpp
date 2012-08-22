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

#ifndef SYMBOL_TABLE_H
#define SYMBOL_TABLE_H 1

#if HAVE_CONFIG_H
#include "config.h"
#endif

/* Data structure representing symbol table record.  */
typedef struct _Record {
	struct _Record *next;	/* Pointer to next record.  */
	char           *name;	/* Symbol name.  */
	char            type;	/* Symbol type ('c' for constant, 'v' for
				 * variable, 'f' for function).  */
	union {
		double          value;	/* Constant or variable value.  */
		double          (*function) (double);	/* Pointer to
							 * function to
							 * calculate its
							 * value.  */
	} data;
	int             flag;	/* Record flag used for symbol table
				 * selective traversal.  */
} Record;

/* Data structure representing symbol table (hash table is used for this
 * purpose). */
typedef struct {
	int             length;	/* Hash table length.  */
	Record         *records;	/* Hash table buckets.  */
	int             reference_count;	/* Reference count for
						 * symbol table (evaluator 
						 * for derivative uses
						 * same symbol table as
						 * evaluator for
						 * corresponding
						 * function). */
} SymbolTable;

/* Create symbol table using specified length of hash table.  */
SymbolTable    *symbol_table_create(int length);

/* Destroy symbol table.  */
void            symbol_table_destroy(SymbolTable * symbol_table);

/* Insert symbol into given symbol table.  Further arguments are symbol
 * name and its type, as well as additional arguments according to symbol
 * type.  Return value is pointer to symbol table record created to
 * represent symbol.  If symbol already in symbol table, pointer to its
 * record is returned immediately. */
Record         *symbol_table_insert(SymbolTable * symbol_table, char *name,
				    char type, ...);

/* Lookup symbol by name from given symbol table.  Pointer to symbol
 * record is returned if symbol found, null pointer otherwise. */
Record         *symbol_table_lookup(SymbolTable * symbol_table,
				    char *name);

/* Clear flag for each symbol table record. */
void            symbol_table_clear_flags(SymbolTable * symbol_table);

/* Count number of flagged records in symbol table. */
int             symbol_table_get_flagged_count(SymbolTable * symbol_table);

/* Fill given array with pointers to records from given symbol table that
 * have flag set.  Further arguments are array to store pointers and array 
 * capacity.  Number of records that are actually put into array is
 * returned. */
int             symbol_table_get_flagged(SymbolTable * symbol_table,
					 Record ** records, int length);

/* Return symbol table pointer to be assigned to variable.  This function
 * should be used instead of simple pointer assignement for proper
 * reference counting.  Users willing to manage reference counts by
 * themselves are free to ignore this function. */
SymbolTable    *symbol_table_assign(SymbolTable * symbol_table);

#endif
