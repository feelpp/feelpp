/** @file parse_binop_rhs.cpp
 *
 *  Code to deal with binary operators. */

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

#include "ex.h"
#include "symbol.h"
#include "mul.h"
#include "add.h"
#include "power.h"
#include "operators.h"
#include "parser.h"
#include "lexer.h"
#include "debug.h"

#include <sstream>
#include <stdexcept>

namespace GiNaC {

/// Make a sum or a product.
static ex make_binop_expr(const int binop, const exvector& args);
/// Check if the token is a binary operator. 
static inline bool is_binop(const int c);
/// Get the precedence of the pending binary operator.
static int get_tok_prec(const int c);

/// binoprhs: ([+*/^-] primary)*
ex parser::parse_binop_rhs(int expr_prec, ex& lhs)
{
	exvector args;
	args.push_back(lhs);
	int binop = -1, orig_binop = -1;
	bool need_sign_flip = false;
	while (1) {
		// check if this is a binop
		if (!is_binop(token)) {
			if (args.size() > 1)
				return make_binop_expr(orig_binop, args);
			else
				return lhs;
		}
		
		// Okay, we know this is a binop.
		if (args.size() == 1)
			orig_binop = token;

		binop = token;

		// If this is a binop that binds at least as tightly as
		// the current binop, consume it, otherwise we are done.
		int tok_prec = get_tok_prec(token);
		if (tok_prec < expr_prec) {
			if (args.size() > 1)
				return make_binop_expr(orig_binop, args);
			else 
				return lhs;
		}

		get_next_tok();  // eat binop

		// Parse the primary expression after the binary operator.
		ex rhs = parse_primary();

		// If binop binds less tightly with rhs than the operator after
		// rhs, let the pending operator take rhs as its lhs.
		int next_prec = get_tok_prec(token);
		if (tok_prec < next_prec)
			rhs = parse_binop_rhs(tok_prec + 1, rhs);

		// previous operator was '+', and current one is '-'
		// (or vice a versa).
		if (need_sign_flip)
			rhs = - rhs;

		args.push_back(rhs);

		// Minimize the number of eval() and ctor calls. This is
		// crucial for a reasonable performance. If the next operator
		// is compatible with the pending one (or the same) don't create
		// the expression and continue collecting operands instead.
		if (binop == token)
			continue;
		else if (binop == '+' && token == '-') {
			need_sign_flip = token != orig_binop;
			continue;
		} else if (binop == '-' && token == '+') {
			need_sign_flip = token != orig_binop;
			continue;
		} else { 
			if (args.size() <= 1)
				bug("binop has " << args.size() << " arguments, expected >= 2");
			lhs = make_binop_expr(orig_binop, args);
			args.clear();
			args.push_back(lhs);
		}
	}
}

extern const numeric* _num_1_p;

static ex make_minus_expr(const exvector& args)
{
	exvector rest_args;
	rest_args.reserve(args.size() - 1);
	std::copy(args.begin() + 1, args.end(), std::back_inserter(rest_args));
	ex rest_base = (new add(rest_args))->setflag(status_flags::dynallocated);
	ex rest = (new mul(rest_base, *_num_1_p))->setflag(status_flags::dynallocated);
	ex ret = (new add(args[0], rest))->setflag(status_flags::dynallocated);
	return ret;
}

static ex make_divide_expr(const exvector& args)
{
	exvector rest_args;
	rest_args.reserve(args.size() - 1);
	std::copy(args.begin() + 1, args.end(), std::back_inserter(rest_args));
	ex rest_base = (new mul(rest_args))->setflag(status_flags::dynallocated);
	ex rest = pow(rest_base, *_num_1_p);
	return (new mul(args[0], rest))->setflag(status_flags::dynallocated);
}

static ex make_binop_expr(const int binop, const exvector& args)
{
	switch (binop) {
		case '+':
			return (new add(args))->setflag(status_flags::dynallocated);
		case '-':
			return make_minus_expr(args);
		case '*':
			return (new mul(args))->setflag(status_flags::dynallocated);
		case '/':
			return make_divide_expr(args);
		case '^':
			if (args.size() != 2)
				throw std::invalid_argument(
						std::string(__func__) 
						+ ": power should have exactly 2 operands");
			return pow(args[0], args[1]);
		default:
			throw std::invalid_argument(
					std::string(__func__) 
					+ ": invalid binary operation: " 
					+ char(binop));
	}
}

static inline bool is_binop(const int c)
{
	switch (c) {
		case '+':
		case '-':
		case '*':
		case '/':
		case '^':
			return true;
		default:
			return false;
	}
}

/// Get the precedence of the pending binary operator.
static int get_tok_prec(const int c)
{
	switch (c) {
		case '+':
		case '-':
			return 20;
		case '*':
		case '/':
			return 40;
		case '^':
			return 60;
		default:
			return -1;
			// means 'this is not a binary operator'
	}
}

} // namespace GiNaC
