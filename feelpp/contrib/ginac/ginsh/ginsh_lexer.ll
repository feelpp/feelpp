/** @file ginsh_lexer.ll
 *
 *  Lexical analyzer definition for ginsh.
 *  This file must be processed with flex. */

/*
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


/*
 *  Definitions
 */

%pointer

%{
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ginsh.h"
#include "ginsh_parser.h"

#define YY_INPUT(buf, result, max_size) (result = ginsh_input(buf, max_size))

// Table of all used symbols
sym_tab syms;

// Type of symbols to generate (real or complex)
unsigned symboltype = domain::complex;

// lex input function
static int ginsh_input(char *buf, int max_size);
%}

	/* Abbreviations */
D	[0-9]
E	[elEL][-+]?{D}+
A	[a-zA-Z_]
AN	[0-9a-zA-Z_]


/*
 *  Lexical rules
 */

%%
[ \t\n]+		/* skip whitespace */
\\$			/* skip line continuations */
"//".*			/* skip comments starting with "//" */
^"#".*			/* skip lines starting with "#" */
^"!".*			system(yytext + 1);	/* execute shell command */

			/* special values */
Pi			yylval = Pi; return T_LITERAL;
Euler			yylval = Euler; return T_LITERAL;
Catalan			yylval = Catalan; return T_LITERAL;
FAIL			yylval = *new fail(); return T_LITERAL;
I			yylval = I; return T_NUMBER;
Digits			yylval = (long)Digits; return T_DIGITS;

			/* keywords */
quit|exit		return T_QUIT;
warranty		return T_WARRANTY;
print			return T_PRINT;
iprint			return T_IPRINT;
print_latex		return T_PRINTLATEX;
print_csrc		return T_PRINTCSRC;
time			return T_TIME;
xyzzy			return T_XYZZY;
inventory		return T_INVENTORY;
look			return T_LOOK;
score			return T_SCORE;
complex_symbols return T_COMPLEX_SYMBOLS;
real_symbols    return T_REAL_SYMBOLS;

			/* comparison */
"=="			return T_EQUAL;
"!="			return T_NOTEQ;
"<="			return T_LESSEQ;
">="			return T_GREATEREQ;

			/* last 1..3 expressions */
\%			return T_QUOTE;
\%\%			return T_QUOTE2;
\%\%\%			return T_QUOTE3;

			/* numbers */
{D}+			|
"#"{D}+"R"{AN}+		|
"#b"([01])+		|
"#o"[0-7]+		|
"#x"[0-9a-fA-F]+	|
{D}+"."{D}*({E})?	|
{D}*"."{D}+({E})?	|
{D}+{E}			yylval = numeric(yytext); return T_NUMBER;

			/* symbols */
{A}{AN}*		{
				sym_tab::const_iterator i = syms.find(yytext);
				if (i == syms.end()) {
					if (symboltype == domain::complex) {
						symbol tmp(yytext);
						syms[yytext] = tmp;
						yylval = tmp;
					} else {
						realsymbol tmp(yytext);
						syms[yytext] = tmp;
						yylval = tmp;
					}
				} else
					yylval = i->second;
				return T_SYMBOL;
			}

			/* wildcards */
\${D}+			yylval = wild(atoi(yytext + 1)); return T_LITERAL;

			/* everything else */
.			return *yytext;

%%


/*
 *  Routines
 */

static int line_length = 0;
static char *line_read = NULL;
static char *line_ptr;

// Input function that uses libreadline for interactive input
static int ginsh_input(char *buf, int max_size)
{
	int result;
#if defined(YY_CURRENT_BUFFER)
	if (YY_CURRENT_BUFFER->yy_is_interactive) {
#else
	if (yy_current_buffer->yy_is_interactive) {
#endif
#ifdef HAVE_LIBREADLINE
		// Do we need to read a new line?
		int actual;
		if (line_length == 0) {

			// Free old line
			if (line_read)
				free(line_read);

			// Read new line, prompt "> "
			line_read = line_ptr = readline("> ");

			// EOF?
			if (!line_read) {
				line_length = 0;
				return YY_NULL;
			}

			// Add non-empty lines to history
			line_length = strlen(line_read) + 1;
			if (line_length > 1)
				add_history(line_read);

			// Reappend trailing '\n' which is stripped by readline()
			line_read[line_length - 1] = '\n';
		}

		// Copy data to lex buffer
		actual = line_length > max_size ? max_size : line_length;
		memcpy(buf, line_ptr, actual);
		line_length -= actual;
		line_ptr += actual;
		result = actual;
#else
		printf("> "); fflush(stdout);
		int c = '*', n;
		for (n = 0; n < max_size && (c = getc(yyin)) != EOF && c != '\n'; ++n)
			buf[n] = (char)c;
		if (c == '\n')
			buf[n++] = (char)c;
		if (c == EOF && ferror(yyin))
			YY_FATAL_ERROR("input in flex scanner failed");
		result = n;
#endif
	} else if (((result = fread(buf, 1, max_size, yyin)) == 0) && ferror(yyin))
		YY_FATAL_ERROR("input in flex scanner failed");

	return result;
}

// List of input files to be processed
int num_files = 0;
char **file_list = NULL;

// EOF encountered, connect to next file. If this was the last file,
// connect to stdin. If this was stdin, terminate the scanner.
int yywrap()
{
	if (yyin == stdin)
		return 1;

	fclose(yyin);
	if (num_files) {
		yyin = fopen(*file_list, "r");
		if (yyin == NULL) {
			cerr << "Can't open " << *file_list << endl;
			return 1;
		}
		num_files--;
		file_list++;
	} else
		yyin = stdin;
	return 0;
}
