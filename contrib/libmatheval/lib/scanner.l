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

%{
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
#  include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include "node.h"
#include "parser.h"
#include "symbol_table.h"

/* Redefine macro to redirect scanner input from string instead of
 * standard input.  */
#define YY_INPUT( buffer, result, max_size ) \
{ result = input_from_string (buffer, max_size); }

/* Variables used to communicate with code using scanner.  */
extern SymbolTable *symbol_table; /* Evaluator symbol table.  */
extern char *input_string; /* String representing function.  */

/* Read next max_size character from string into buffer.  */
static int input_from_string (char *buffer, int max_size);
%}

/* Token definitions.  */
whitespace [ \t]+
digit [0-9]
number ({digit}+|{digit}+\.{digit}*|{digit}*\.{digit}+)([Ee][+\-]?{digit}+)?
constant "e"|"log2e"|"log10e"|"ln2"|"ln10"|"pi"|"pi_2"|"pi_4"|"1_pi"|"2_pi"|"2_sqrtpi"|"sqrt2"|"sqrt1_2"
function "exp"|"log"|"sqrt"|"sin"|"cos"|"tan"|"cot"|"sec"|"csc"|"asin"|"acos"|"atan"|"acot"|"asec"|"acsc"|"sinh"|"cosh"|"tanh"|"coth"|"sech"|"csch"|"asinh"|"acosh"|"atanh"|"acoth"|"asech"|"acsch"|"abs"|"step"|"delta"|"nandelta"|"erf"
identifier [_a-zA-Z][[_a-zA-Z0-9]*

%%

{whitespace}

{number} {
        /* Create node representing constant with appropriate value.  */
        yylval.node = node_create ('n', atof (yytext));
        return NUMBER;
}

{constant} {
        Record *record; /* Symbol table record.  */

        /* Find symbol table record corresponding to constant name.  */
        record = symbol_table_lookup (symbol_table, yytext);
        yylval.node = node_create ('c', record);
        return CONSTANT;
}

{function} {
        /* Find symbol table record corresponding to function name.  */
        yylval.record = symbol_table_lookup (symbol_table, yytext);
        return FUNCTION;
}

{identifier} {
        Record *record; /* Symbol table record.  */

        /* Insert variable into symbol table.  */
        record = symbol_table_insert (symbol_table, yytext, 'v');
        yylval.node = node_create ('v', record);
        return VARIABLE;
}

"+" {
        return '+';
}

"-" {
        return '-';
}

"*" {
        return '*';
}

"/" {
        return '/';
}

"^" {
        return '^';
}

"(" {
        return '(';
}

")" {
        return ')';
}

"\n" {
        return '\n';
}

%%

static int input_from_string (char *buffer, int max_size)
{
        int count; /* Count of characters to copy from input string to buffer.  */

        /* Calculate count of characters to copy.  */
        count = strlen (input_string);
        if (count > max_size)
                count = max_size;

        /* Perform copy operation and update input string.  */
        memcpy(buffer, input_string, count);
        input_string += count;

        return count;
}

void input_reset()
{
        /* Forget contents of input buffer.  */
        YY_FLUSH_BUFFER;
}
