/** @file ginsh.h
 *
 *  Global definitions for ginsh.
 *
 *  GiNaC Copyright (C) 1999-2011 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef GINSH_H
#define GINSH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <map>
#include <iostream>
#include <string>

using namespace std;

#ifdef HAVE_READLINE_READLINE_H
extern "C" {
#include <readline/readline.h>
}
#endif

#ifdef HAVE_READLINE_HISTORY_H
extern "C" {
#include <readline/history.h>
}
#endif

#ifdef IN_GINAC
#include "ginac.h"
#else
#include <ginac/ginac.h>
#endif

using namespace GiNaC;

// yacc stack type
#define YYSTYPE ex

// lex functions/variables
extern int yyerror(const char *s);
extern int yylex(void);
extern char *yytext;
extern FILE *yyin;

// List of input files to be processed
extern int num_files;
extern char **file_list;

// Table of all used symbols
typedef map<string, ex> sym_tab;
extern sym_tab syms;

// Type of symbols to generate (real or complex)
extern unsigned symboltype;

#endif
