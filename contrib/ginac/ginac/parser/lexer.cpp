/** @file lexer.cpp
 *
 *  Implementation of GiNaC's lexer. */

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

#include "lexer.h"
#include "compiler.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>

namespace GiNaC {

/// Skip to the end of line
static int skipline(std::istream* s);
/// Skip to the next non-whitespace character
static int skipspace(std::istream* s, int c, std::size_t& line);
/// Check if the identifier is predefined literal
static bool literal_p(const std::string& name);

/// gettok - Return the next token from standard input.
int lexer::gettok()
{
	// Skip any whitespace.
	c = skipspace(input, c, line_num);

	// identifier: [a-zA-Z][a-zA-Z0-9_]*
	if (isalpha(c)) { 
		str = c;
		do {
			c = input->get();
			if ( isalnum(c) || c=='_' )
				str += c;
			else
				break;
		} while (true);
		if (unlikely(literal_p(str)))
			return token_type::literal;
		else
			return token_type::identifier;
	}

	// Number: [0-9]+([.][0-9]*(eE[+-][0-9]+)*)*
	if (isdigit(c) || c == '.') {
		str = "";
		do {
			str += c;
			c = input->get();
		} while (isdigit(c) || c == '.');
		if (c == 'E' || c == 'e') {
			str += 'E';
			c = input->get();
			if (isdigit(c))
				str += '+';
			do {
				str += c;
				c = input->get();
			} while (isdigit(c));
		}
		return token_type::number;
	}

	// Comment until end of line.
	if (c == '#') {
		c = skipline(input);
		++line_num;
		if (c != EOF)
			return gettok();
	}

	// Check for end of file.  Don't eat the EOF.
	if (c == EOF)
		return token_type::eof;

	// Otherwise, just return the character as its ascii value.
	int current = c;
	c = input->get();
	return current;
}

static int skipline(std::istream* s)
{
	int c;
	do {
		c = s->get();
	} while (c != EOF && c != '\n' && c != '\r');
	return c;
}

static int skipspace(std::istream* s, int c, std::size_t& line)
{
	while (isspace(c)) {
		if (c == '\n')
			++line;
		c = s->get();
	}
	return c;
}

static bool literal_p(const std::string& name)
{
	if (name == "I")
		return true;
	else if (name == "Pi")
		return true;
	else if (name == "Euler")
		return true;
	else if (name == "Catalan")
		return true;
	else
		return false; 
}

lexer::lexer(std::istream* in, std::ostream* out, std::ostream* err)
{
	if (in)
		input = in;
	else
		in = &std::cin;

	if (out)
		output = out;
	else
		output = &std::cout;

	if (err)
		error = err;
	else
		error = &std::cerr;

	c = ' ';
	str = "";
	line_num = 0;
	column = 0;
}

lexer::~lexer() { }

void lexer::switch_input(std::istream* in)
{
	input = in;
	line_num = 0;
	column = 0;
	c = ' ';
}

/// Symbolic name of current token (for error reporting)
std::string lexer::tok2str(const int tok) const
{
	switch (tok) {
		case lexer::token_type::identifier:
		case lexer::token_type::number:
			return std::string("\"") + str + "\"";
		case lexer::token_type::eof:
			return std::string("EOF");
		default:
			return std::string("\"") + char(tok) + "\"";
	}
}

} // namespace GiNaC
