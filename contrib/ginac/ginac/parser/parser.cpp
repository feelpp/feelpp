/** @file parser.cpp
 *
 *  Implementation of GiNaC's parser. */

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

#include "parser.h"
#include "lst.h"
#include "lexer.h"
#include "debug.h"
#include "mul.h"
#include "constant.h"
#include "function.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h> // for uintptr_t
#endif
#include <sstream>
#include <stdexcept>

namespace GiNaC {

// <KLUDGE>
// Find out if ptr is a pointer to a function or a specially crafted integer.
// It's possible to distinguish between these because functions are aligned.
// Returns true if ptr is a pointer and false otherwise.
static bool decode_serial(unsigned& serial, const reader_func ptr)
{
	uintptr_t u = (uintptr_t)(void *)ptr;
	if (u & 1) {
		u >>= 1;
		serial = (unsigned)u;
		return false;
	}
	return true;
}

// Figures out if ptr is a pointer to function or a serial of GiNaC function.
// In the former case calls that function, in the latter case constructs
// GiNaC function with corresponding serial and arguments.
static ex dispatch_reader_fcn(const reader_func ptr, const exvector& args)
{
	unsigned serial = 0; // dear gcc, could you please shut up?
	bool is_ptr = decode_serial(serial, ptr);
	if (is_ptr)
		return ptr(args);
	else
		return function(serial, args);
}
// </KLUDGE>


/// identifier_expr:  identifier |  identifier '(' expression* ')'
ex parser::parse_identifier_expr()
{
	std::string name = scanner->str;
	get_next_tok();  // eat identifier.

	if (token != '(') // symbol
		return find_or_insert_symbol(name, syms, strict);

	// function/ctor call.
	get_next_tok();  // eat (
	exvector args;
	if (token != ')') {
		while (true) {
			ex e = parse_expression();
			args.push_back(e);

			if (token == ')')
				break;

			if (token != ',')
				Parse_error("expected ')' or ',' in argument list");

			get_next_tok();
		}
	}
	// Eat the ')'.
	get_next_tok();
	prototype the_prototype = make_pair(name, args.size());
	prototype_table::const_iterator reader = funcs.find(the_prototype);
	if (reader == funcs.end()) {
		Parse_error_("no function \"" << name << "\" with " <<
			     args.size() << " arguments");
	}
	// reader->second might be a pointer to a C++ function or a specially
	// crafted serial of a GiNaC::function.
	ex ret = dispatch_reader_fcn(reader->second, args);
	return ret;
}

/// paren_expr:  '(' expression ')'
ex parser::parse_paren_expr()
{
	get_next_tok();  // eat (.
	ex e = parse_expression();

	if (token != ')')
		Parse_error("expected ')'");
	get_next_tok();  // eat ).
	return e;
}

/// lst_expr:  '{' expression { ',' expression } '}'
ex parser::parse_lst_expr()
{
	get_next_tok();  // eat {.

	lst list;
	if (token != '}') {
		while (true) {
			ex e = parse_expression(); // expression();
			list.append(e);

			if (token == '}') {
				break;
			}

			if (token != ',') {
				Parse_error("expected '}'");
			}

			get_next_tok();  // eat ','.
		}
	}
	// Eat the '}'.
	get_next_tok();

	return list;
}

extern const ex _ex0;

/// unary_expr: [+-] expression
ex parser::parse_unary_expr()
{
	// Unlike most other parse_* method this one does NOT consume
	// current token so parse_binop_rhs() knows what kind of operator
	// is being parsed.
	
	// There are different kinds of expressions which need to be handled:
	// -a+b 
	// -(a) 
	// +a
	// +(a)
	// Delegete the work to parse_binop_rhs(), otherwise we end up
	// duplicating it here. 
	ex lhs = _ex0; // silly trick
	ex e = parse_binop_rhs(0, lhs);
	return e;
}

/// primary: identifier_expr | number_expr | paren_expr | unary_expr 
ex parser::parse_primary() 
{
	switch (token) {
		case lexer::token_type::identifier:
			 return parse_identifier_expr();
		case lexer::token_type::number:
			 return parse_number_expr();
		case '(': 
			 return parse_paren_expr();
		case '{': 
			 return parse_lst_expr();
		case '-':
		case '+':
			 return parse_unary_expr();
		case lexer::token_type::literal:
			 return parse_literal_expr();
		case lexer::token_type::eof:
		default:
			 Parse_error("unexpected token");
	}
}

/// expression ::= primary binoprhs
ex parser::parse_expression() 
{
	ex lhs = parse_primary();
	ex res = parse_binop_rhs(0, lhs);
	return res;
}

/// number_expr: number
ex parser::parse_number_expr()
{
	ex n = numeric(scanner->str.c_str());
	get_next_tok(); // consume the number
	return n;
}

/// literal_expr: 'I' | 'Pi' | 'Euler' | 'Catalan'
ex parser::parse_literal_expr()
{
	get_next_tok(); // consume the literal
	if (scanner->str == "I")
		return I;
	else if (scanner->str == "Pi")
		return Pi;
	else if (scanner->str == "Euler")
		return Euler;
	else if (scanner->str == "Catalan")
		return Catalan;
	bug("unknown literal: \"" << scanner->str << "\"");
}

ex parser::operator()(std::istream& input)
{
	scanner->switch_input(&input);
	get_next_tok();
	ex ret = parse_expression();
	// parse_expression() stops if it encounters an unknown token.
	// This is not a bug: since the parser is recursive checking
	// whether the next token is valid is responsibility of the caller.
	// Hence make sure nothing is left in the stream:
	if (token != lexer::token_type::eof)
		Parse_error("expected EOF");

	return ret;
}

ex parser::operator()(const std::string& input)
{
	std::istringstream is(input);
	ex ret = operator()(is);
	return ret;
}

int parser::get_next_tok()
{
	token = scanner->gettok();
	return token;
}

parser::parser(const symtab& syms_, const bool strict_,
	       const prototype_table& funcs_) : strict(strict_),
	funcs(funcs_), syms(syms_)
{
	scanner = new lexer();
}

parser::~parser()
{
	delete scanner;
}

} // namespace GiNaC
