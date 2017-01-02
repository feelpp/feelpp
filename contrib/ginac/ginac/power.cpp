/** @file power.cpp
 *
 *  Implementation of GiNaC's symbolic exponentiation (basis^exponent). */

/*
 *  GiNaC Copyright (C) 1999-2016 Johannes Gutenberg University Mainz, Germany
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

#include "power.h"
#include "expairseq.h"
#include "add.h"
#include "mul.h"
#include "ncmul.h"
#include "numeric.h"
#include "constant.h"
#include "operators.h"
#include "inifcns.h" // for log() in power::derivative()
#include "matrix.h"
#include "indexed.h"
#include "symbol.h"
#include "lst.h"
#include "archive.h"
#include "utils.h"
#include "relational.h"
#include "compiler.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <algorithm>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(power, basic,
  print_func<print_dflt>(&power::do_print_dflt).
  print_func<print_latex>(&power::do_print_latex).
  print_func<print_csrc>(&power::do_print_csrc).
  print_func<print_python>(&power::do_print_python).
  print_func<print_python_repr>(&power::do_print_python_repr).
  print_func<print_csrc_cl_N>(&power::do_print_csrc_cl_N))

//////////
// default constructor
//////////

power::power() { }

//////////
// other constructors
//////////

// all inlined

//////////
// archiving
//////////

void power::read_archive(const archive_node &n, lst &sym_lst)
{
	inherited::read_archive(n, sym_lst);
	n.find_ex("basis", basis, sym_lst);
	n.find_ex("exponent", exponent, sym_lst);
}

void power::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_ex("basis", basis);
	n.add_ex("exponent", exponent);
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void power::print_power(const print_context & c, const char *powersymbol, const char *openbrace, const char *closebrace, unsigned level) const
{
	// Ordinary output of powers using '^' or '**'
	if (precedence() <= level)
		c.s << openbrace << '(';
	basis.print(c, precedence());
	c.s << powersymbol;
	c.s << openbrace;
	exponent.print(c, precedence());
	c.s << closebrace;
	if (precedence() <= level)
		c.s << ')' << closebrace;
}

void power::do_print_dflt(const print_dflt & c, unsigned level) const
{
	if (exponent.is_equal(_ex1_2)) {

		// Square roots are printed in a special way
		c.s << "sqrt(";
		basis.print(c);
		c.s << ')';

	} else
		print_power(c, "^", "", "", level);
}

void power::do_print_latex(const print_latex & c, unsigned level) const
{
	if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_negative()) {

		// Powers with negative numeric exponents are printed as fractions
		c.s << "\\frac{1}{";
		power(basis, -exponent).eval().print(c);
		c.s << '}';

	} else if (exponent.is_equal(_ex1_2)) {

		// Square roots are printed in a special way
		c.s << "\\sqrt{";
		basis.print(c);
		c.s << '}';

	} else
		print_power(c, "^", "{", "}", level);
}

static void print_sym_pow(const print_context & c, const symbol &x, int exp)
{
	// Optimal output of integer powers of symbols to aid compiler CSE.
	// C.f. ISO/IEC 14882:2011, section 1.9 [intro execution], paragraph 15
	// to learn why such a parenthesation is really necessary.
	if (exp == 1) {
		x.print(c);
	} else if (exp == 2) {
		x.print(c);
		c.s << "*";
		x.print(c);
	} else if (exp & 1) {
		x.print(c);
		c.s << "*";
		print_sym_pow(c, x, exp-1);
	} else {
		c.s << "(";
		print_sym_pow(c, x, exp >> 1);
		c.s << ")*(";
		print_sym_pow(c, x, exp >> 1);
		c.s << ")";
	}
}

void power::do_print_csrc_cl_N(const print_csrc_cl_N& c, unsigned level) const
{
	if (exponent.is_equal(_ex_1)) {
		c.s << "recip(";
		basis.print(c);
		c.s << ')';
		return;
	}
	c.s << "expt(";
	basis.print(c);
	c.s << ", ";
	exponent.print(c);
	c.s << ')';
}

void power::do_print_csrc(const print_csrc & c, unsigned level) const
{
	// Integer powers of symbols are printed in a special, optimized way
	if (exponent.info(info_flags::integer) &&
	    (is_a<symbol>(basis) || is_a<constant>(basis))) {
		int exp = ex_to<numeric>(exponent).to_int();
		if (exp > 0)
			c.s << '(';
		else {
			exp = -exp;
			c.s << "1.0/(";
		}
		print_sym_pow(c, ex_to<symbol>(basis), exp);
		c.s << ')';

	// <expr>^-1 is printed as "1.0/<expr>" or with the recip() function of CLN
	} else if (exponent.is_equal(_ex_1)) {
		c.s << "1.0/(";
		basis.print(c);
		c.s << ')';

	// Otherwise, use the pow() function
	} else {
		c.s << "pow(";
		basis.print(c);
		c.s << ',';
		exponent.print(c);
		c.s << ')';
	}
}

void power::do_print_python(const print_python & c, unsigned level) const
{
	print_power(c, "**", "", "", level);
}

void power::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << '(';
	basis.print(c);
	c.s << ',';
	exponent.print(c);
	c.s << ')';
}

bool power::info(unsigned inf) const
{
	switch (inf) {
		case info_flags::polynomial:
		case info_flags::integer_polynomial:
		case info_flags::cinteger_polynomial:
		case info_flags::rational_polynomial:
		case info_flags::crational_polynomial:
			return basis.info(inf) && exponent.info(info_flags::nonnegint);
		case info_flags::rational_function:
			return basis.info(inf) && exponent.info(info_flags::integer);
		case info_flags::real:
			return basis.info(inf) && exponent.info(info_flags::integer);
		case info_flags::expanded:
			return (flags & status_flags::expanded);
		case info_flags::positive:
			return basis.info(info_flags::positive) && exponent.info(info_flags::real);
		case info_flags::nonnegative:
			return (basis.info(info_flags::positive) && exponent.info(info_flags::real)) ||
			       (basis.info(info_flags::real) && exponent.info(info_flags::even));
		case info_flags::has_indices: {
			if (flags & status_flags::has_indices)
				return true;
			else if (flags & status_flags::has_no_indices)
				return false;
			else if (basis.info(info_flags::has_indices)) {
				setflag(status_flags::has_indices);
				clearflag(status_flags::has_no_indices);
				return true;
			} else {
				clearflag(status_flags::has_indices);
				setflag(status_flags::has_no_indices);
				return false;
			}
		}
	}
	return inherited::info(inf);
}

size_t power::nops() const
{
	return 2;
}

ex power::op(size_t i) const
{
	GINAC_ASSERT(i<2);

	return i==0 ? basis : exponent;
}

ex power::map(map_function & f) const
{
	const ex &mapped_basis = f(basis);
	const ex &mapped_exponent = f(exponent);

	if (!are_ex_trivially_equal(basis, mapped_basis)
	 || !are_ex_trivially_equal(exponent, mapped_exponent))
		return dynallocate<power>(mapped_basis, mapped_exponent);
	else
		return *this;
}

bool power::is_polynomial(const ex & var) const
{
	if (basis.is_polynomial(var)) {
		if (basis.has(var))
			// basis is non-constant polynomial in var
			return exponent.info(info_flags::nonnegint);
		else
			// basis is constant in var
			return !exponent.has(var);
	}
	// basis is a non-polynomial function of var
	return false;
}

int power::degree(const ex & s) const
{
	if (is_equal(ex_to<basic>(s)))
		return 1;
	else if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_integer()) {
		if (basis.is_equal(s))
			return ex_to<numeric>(exponent).to_int();
		else
			return basis.degree(s) * ex_to<numeric>(exponent).to_int();
	} else if (basis.has(s))
		throw(std::runtime_error("power::degree(): undefined degree because of non-integer exponent"));
	else
		return 0;
}

int power::ldegree(const ex & s) const 
{
	if (is_equal(ex_to<basic>(s)))
		return 1;
	else if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_integer()) {
		if (basis.is_equal(s))
			return ex_to<numeric>(exponent).to_int();
		else
			return basis.ldegree(s) * ex_to<numeric>(exponent).to_int();
	} else if (basis.has(s))
		throw(std::runtime_error("power::ldegree(): undefined degree because of non-integer exponent"));
	else
		return 0;
}

ex power::coeff(const ex & s, int n) const
{
	if (is_equal(ex_to<basic>(s)))
		return n==1 ? _ex1 : _ex0;
	else if (!basis.is_equal(s)) {
		// basis not equal to s
		if (n == 0)
			return *this;
		else
			return _ex0;
	} else {
		// basis equal to s
		if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_integer()) {
			// integer exponent
			int int_exp = ex_to<numeric>(exponent).to_int();
			if (n == int_exp)
				return _ex1;
			else
				return _ex0;
		} else {
			// non-integer exponents are treated as zero
			if (n == 0)
				return *this;
			else
				return _ex0;
		}
	}
}

/** Perform automatic term rewriting rules in this class.  In the following
 *  x, x1, x2,... stand for a symbolic variables of type ex and c, c1, c2...
 *  stand for such expressions that contain a plain number.
 *  - ^(x,0) -> 1  (also handles ^(0,0))
 *  - ^(x,1) -> x
 *  - ^(0,c) -> 0 or exception  (depending on the real part of c)
 *  - ^(1,x) -> 1
 *  - ^(c1,c2) -> *(c1^n,c1^(c2-n))  (so that 0<(c2-n)<1, try to evaluate roots, possibly in numerator and denominator of c1)
 *  - ^(^(x,c1),c2) -> ^(x,c1*c2)  if x is positive and c1 is real.
 *  - ^(^(x,c1),c2) -> ^(x,c1*c2)  (c2 integer or -1 < c1 <= 1 or (c1=-1 and c2>0), case c1=1 should not happen, see below!)
 *  - ^(*(x,y,z),c) -> *(x^c,y^c,z^c)  (if c integer)
 *  - ^(*(x,c1),c2) -> ^(x,c2)*c1^c2  (c1>0)
 *  - ^(*(x,c1),c2) -> ^(-x,c2)*c1^c2  (c1<0)
 */
ex power::eval() const
{
	if (flags & status_flags::evaluated)
		return *this;

	const numeric *num_basis = nullptr;
	const numeric *num_exponent = nullptr;

	if (is_exactly_a<numeric>(basis)) {
		num_basis = &ex_to<numeric>(basis);
	}
	if (is_exactly_a<numeric>(exponent)) {
		num_exponent = &ex_to<numeric>(exponent);
	}
	
	// ^(x,0) -> 1  (0^0 also handled here)
	if (exponent.is_zero()) {
		if (basis.is_zero())
			throw (std::domain_error("power::eval(): pow(0,0) is undefined"));
		else
			return _ex1;
	}
	
	// ^(x,1) -> x
	if (exponent.is_equal(_ex1))
		return basis;

	// ^(0,c1) -> 0 or exception  (depending on real value of c1)
	if (basis.is_zero() && num_exponent) {
		if ((num_exponent->real()).is_zero())
			throw (std::domain_error("power::eval(): pow(0,I) is undefined"));
		else if ((num_exponent->real()).is_negative())
			throw (pole_error("power::eval(): division by zero",1));
		else
			return _ex0;
	}

	// ^(1,x) -> 1
	if (basis.is_equal(_ex1))
		return _ex1;

	// power of a function calculated by separate rules defined for this function
	if (is_exactly_a<function>(basis))
		return ex_to<function>(basis).power(exponent);

	// Turn (x^c)^d into x^(c*d) in the case that x is positive and c is real.
	if (is_exactly_a<power>(basis) && basis.op(0).info(info_flags::positive) && basis.op(1).info(info_flags::real))
		return dynallocate<power>(basis.op(0), basis.op(1) * exponent);

	if ( num_exponent ) {

		// ^(c1,c2) -> c1^c2  (c1, c2 numeric(),
		// except if c1,c2 are rational, but c1^c2 is not)
		if ( num_basis ) {
			const bool basis_is_crational = num_basis->is_crational();
			const bool exponent_is_crational = num_exponent->is_crational();
			if (!basis_is_crational || !exponent_is_crational) {
				// return a plain float
				return dynallocate<numeric>(num_basis->power(*num_exponent));
			}

			const numeric res = num_basis->power(*num_exponent);
			if (res.is_crational()) {
				return res;
			}
			GINAC_ASSERT(!num_exponent->is_integer());  // has been handled by now

			// ^(c1,n/m) -> *(c1^q,c1^(n/m-q)), 0<(n/m-q)<1, q integer
			if (basis_is_crational && exponent_is_crational
			    && num_exponent->is_real()
			    && !num_exponent->is_integer()) {
				const numeric n = num_exponent->numer();
				const numeric m = num_exponent->denom();
				numeric r;
				numeric q = iquo(n, m, r);
				if (r.is_negative()) {
					r += m;
					--q;
				}
				if (q.is_zero()) {  // the exponent was in the allowed range 0<(n/m)<1
					if (num_basis->is_rational() && !num_basis->is_integer()) {
						// try it for numerator and denominator separately, in order to
						// partially simplify things like (5/8)^(1/3) -> 1/2*5^(1/3)
						const numeric bnum = num_basis->numer();
						const numeric bden = num_basis->denom();
						const numeric res_bnum = bnum.power(*num_exponent);
						const numeric res_bden = bden.power(*num_exponent);
						if (res_bnum.is_integer())
							return dynallocate<mul>(dynallocate<power>(bden,-*num_exponent),res_bnum).setflag(status_flags::evaluated);
						if (res_bden.is_integer())
							return dynallocate<mul>(dynallocate<power>(bnum,*num_exponent),res_bden.inverse()).setflag(status_flags::evaluated);
					}
					return this->hold();
				} else {
					// assemble resulting product, but allowing for a re-evaluation,
					// because otherwise we'll end up with something like
					//    (7/8)^(4/3)  ->  7/8*(1/2*7^(1/3))
					// instead of 7/16*7^(1/3).
					return pow(basis, r.div(m)) * pow(basis, q);
				}
			}
		}
	
		// ^(^(x,c1),c2) -> ^(x,c1*c2)
		// (c1, c2 numeric(), c2 integer or -1 < c1 <= 1 or (c1=-1 and c2>0),
		// case c1==1 should not happen, see below!)
		if (is_exactly_a<power>(basis)) {
			const power & sub_power = ex_to<power>(basis);
			const ex & sub_basis = sub_power.basis;
			const ex & sub_exponent = sub_power.exponent;
			if (is_exactly_a<numeric>(sub_exponent)) {
				const numeric & num_sub_exponent = ex_to<numeric>(sub_exponent);
				GINAC_ASSERT(num_sub_exponent!=numeric(1));
				if (num_exponent->is_integer() || (abs(num_sub_exponent) - (*_num1_p)).is_negative() ||
				    (num_sub_exponent == *_num_1_p && num_exponent->is_positive())) {
					return dynallocate<power>(sub_basis, num_sub_exponent.mul(*num_exponent));
				}
			}
		}
	
		// ^(*(x,y,z),c1) -> *(x^c1,y^c1,z^c1) (c1 integer)
		if (num_exponent->is_integer() && is_exactly_a<mul>(basis)) {
			return expand_mul(ex_to<mul>(basis), *num_exponent, false);
		}

		// (2*x + 6*y)^(-4) -> 1/16*(x + 3*y)^(-4)
		if (num_exponent->is_integer() && is_exactly_a<add>(basis)) {
			numeric icont = basis.integer_content();
			const numeric lead_coeff = 
				ex_to<numeric>(ex_to<add>(basis).seq.begin()->coeff).div(icont);

			const bool canonicalizable = lead_coeff.is_integer();
			const bool unit_normal = lead_coeff.is_pos_integer();
			if (canonicalizable && (! unit_normal))
				icont = icont.mul(*_num_1_p);
			
			if (canonicalizable && (icont != *_num1_p)) {
				const add& addref = ex_to<add>(basis);
				add & addp = dynallocate<add>(addref);
				addp.clearflag(status_flags::hash_calculated);
				addp.overall_coeff = ex_to<numeric>(addp.overall_coeff).div_dyn(icont);
				for (auto & i : addp.seq)
					i.coeff = ex_to<numeric>(i.coeff).div_dyn(icont);

				const numeric c = icont.power(*num_exponent);
				if (likely(c != *_num1_p))
					return dynallocate<mul>(dynallocate<power>(addp, *num_exponent), c);
				else
					return dynallocate<power>(addp, *num_exponent);
			}
		}

		// ^(*(...,x;c1),c2) -> *(^(*(...,x;1),c2),c1^c2)  (c1, c2 numeric(), c1>0)
		// ^(*(...,x;c1),c2) -> *(^(*(...,x;-1),c2),(-c1)^c2)  (c1, c2 numeric(), c1<0)
		if (is_exactly_a<mul>(basis)) {
			GINAC_ASSERT(!num_exponent->is_integer()); // should have been handled above
			const mul & mulref = ex_to<mul>(basis);
			if (!mulref.overall_coeff.is_equal(_ex1)) {
				const numeric & num_coeff = ex_to<numeric>(mulref.overall_coeff);
				if (num_coeff.is_real()) {
					if (num_coeff.is_positive()) {
						mul & mulp = dynallocate<mul>(mulref);
						mulp.overall_coeff = _ex1;
						mulp.clearflag(status_flags::evaluated | status_flags::hash_calculated);
						return dynallocate<mul>(dynallocate<power>(mulp, exponent),
						                        dynallocate<power>(num_coeff, *num_exponent));
					} else {
						GINAC_ASSERT(num_coeff.compare(*_num0_p)<0);
						if (!num_coeff.is_equal(*_num_1_p)) {
							mul & mulp = dynallocate<mul>(mulref);
							mulp.overall_coeff = _ex_1;
							mulp.clearflag(status_flags::evaluated | status_flags::hash_calculated);
							return dynallocate<mul>(dynallocate<power>(mulp, exponent),
							                        dynallocate<power>(abs(num_coeff), *num_exponent));
						}
					}
				}
			}
		}

		// ^(nc,c1) -> ncmul(nc,nc,...) (c1 positive integer, unless nc is a matrix)
		if (num_exponent->is_pos_integer() &&
		    basis.return_type() != return_types::commutative &&
		    !is_a<matrix>(basis)) {
			return ncmul(exvector(num_exponent->to_int(), basis));
		}
	}

	return this->hold();
}

ex power::evalf() const
{
	ex ebasis = basis.evalf();
	ex eexponent;
	
	if (!is_exactly_a<numeric>(exponent))
		eexponent = exponent.evalf();
	else
		eexponent = exponent;

	return dynallocate<power>(ebasis, eexponent);
}

ex power::evalm() const
{
	const ex ebasis = basis.evalm();
	const ex eexponent = exponent.evalm();
	if (is_a<matrix>(ebasis)) {
		if (is_exactly_a<numeric>(eexponent)) {
			return dynallocate<matrix>(ex_to<matrix>(ebasis).pow(eexponent));
		}
	}
	return dynallocate<power>(ebasis, eexponent);
}

bool power::has(const ex & other, unsigned options) const
{
	if (!(options & has_options::algebraic))
		return basic::has(other, options);
	if (!is_a<power>(other))
		return basic::has(other, options);
	if (!exponent.info(info_flags::integer) ||
	    !other.op(1).info(info_flags::integer))
		return basic::has(other, options);
	if (exponent.info(info_flags::posint) &&
	    other.op(1).info(info_flags::posint) &&
	    ex_to<numeric>(exponent) > ex_to<numeric>(other.op(1)) &&
	    basis.match(other.op(0)))
		return true;
	if (exponent.info(info_flags::negint) &&
	    other.op(1).info(info_flags::negint) &&
	    ex_to<numeric>(exponent) < ex_to<numeric>(other.op(1)) &&
	    basis.match(other.op(0)))
		return true;
	return basic::has(other, options);
}

// from mul.cpp
extern bool tryfactsubs(const ex &, const ex &, int &, exmap&);

ex power::subs(const exmap & m, unsigned options) const
{	
	const ex &subsed_basis = basis.subs(m, options);
	const ex &subsed_exponent = exponent.subs(m, options);

	if (!are_ex_trivially_equal(basis, subsed_basis)
	 || !are_ex_trivially_equal(exponent, subsed_exponent)) 
		return power(subsed_basis, subsed_exponent).subs_one_level(m, options);

	if (!(options & subs_options::algebraic))
		return subs_one_level(m, options);

	for (auto & it : m) {
		int nummatches = std::numeric_limits<int>::max();
		exmap repls;
		if (tryfactsubs(*this, it.first, nummatches, repls)) {
			ex anum = it.second.subs(repls, subs_options::no_pattern);
			ex aden = it.first.subs(repls, subs_options::no_pattern);
			ex result = (*this) * pow(anum/aden, nummatches);
			return (ex_to<basic>(result)).subs_one_level(m, options);
		}
	}

	return subs_one_level(m, options);
}

ex power::eval_ncmul(const exvector & v) const
{
	return inherited::eval_ncmul(v);
}

ex power::conjugate() const
{
	// conjugate(pow(x,y))==pow(conjugate(x),conjugate(y)) unless on the
	// branch cut which runs along the negative real axis.
	if (basis.info(info_flags::positive)) {
		ex newexponent = exponent.conjugate();
		if (are_ex_trivially_equal(exponent, newexponent)) {
			return *this;
		}
		return dynallocate<power>(basis, newexponent);
	}
	if (exponent.info(info_flags::integer)) {
		ex newbasis = basis.conjugate();
		if (are_ex_trivially_equal(basis, newbasis)) {
			return *this;
		}
		return dynallocate<power>(newbasis, exponent);
	}
	return conjugate_function(*this).hold();
}

ex power::real_part() const
{
	// basis == a+I*b, exponent == c+I*d
	const ex a = basis.real_part();
	const ex c = exponent.real_part();
	if (basis.is_equal(a) && exponent.is_equal(c)) {
		// Re(a^c)
		return *this;
	}

	const ex b = basis.imag_part();
	if (exponent.info(info_flags::integer)) {
		// Re((a+I*b)^c)  w/  c ∈ ℤ
		long N = ex_to<numeric>(c).to_long();
		// Use real terms in Binomial expansion to construct
		// Re(expand(pow(a+I*b, N))).
		long NN = N > 0 ? N : -N;
		ex numer = N > 0 ? _ex1 : pow(pow(a,2) + pow(b,2), NN);
		ex result = 0;
		for (long n = 0; n <= NN; n += 2) {
			ex term = binomial(NN, n) * pow(a, NN-n) * pow(b, n) / numer;
			if (n % 4 == 0) {
				result += term;  // sign: I^n w/ n == 4*m
			} else {
				result -= term;  // sign: I^n w/ n == 4*m+2
			}
		}
		return result;
	}

	// Re((a+I*b)^(c+I*d))
	const ex d = exponent.imag_part();
	return pow(abs(basis),c) * exp(-d*atan2(b,a)) * cos(c*atan2(b,a)+d*log(abs(basis)));
}

ex power::imag_part() const
{
	// basis == a+I*b, exponent == c+I*d
	const ex a = basis.real_part();
	const ex c = exponent.real_part();
	if (basis.is_equal(a) && exponent.is_equal(c)) {
		// Im(a^c)
		return 0;
	}

	const ex b = basis.imag_part();
	if (exponent.info(info_flags::integer)) {
		// Im((a+I*b)^c)  w/  c ∈ ℤ
		long N = ex_to<numeric>(c).to_long();
		// Use imaginary terms in Binomial expansion to construct
		// Im(expand(pow(a+I*b, N))).
		long p = N > 0 ? 1 : 3;  // modulus for positive sign
		long NN = N > 0 ? N : -N;
		ex numer = N > 0 ? _ex1 : pow(pow(a,2) + pow(b,2), NN);
		ex result = 0;
		for (long n = 1; n <= NN; n += 2) {
			ex term = binomial(NN, n) * pow(a, NN-n) * pow(b, n) / numer;
			if (n % 4 == p) {
				result += term;  // sign: I^n w/ n == 4*m+p
			} else {
				result -= term;  // sign: I^n w/ n == 4*m+2+p
			}
		}
		return result;
	}

	// Im((a+I*b)^(c+I*d))
	const ex d = exponent.imag_part();
	return pow(abs(basis),c) * exp(-d*atan2(b,a)) * sin(c*atan2(b,a)+d*log(abs(basis)));
}

// protected

/** Implementation of ex::diff() for a power.
 *  @see ex::diff */
ex power::derivative(const symbol & s) const
{
	if (is_a<numeric>(exponent)) {
		// D(b^r) = r * b^(r-1) * D(b) (faster than the formula below)
		const epvector newseq = {expair(basis, exponent - _ex1), expair(basis.diff(s), _ex1)};
		return dynallocate<mul>(std::move(newseq), exponent);
	} else {
		// D(b^e) = b^e * (D(e)*ln(b) + e*D(b)/b)
		return *this * (exponent.diff(s)*log(basis) + exponent*basis.diff(s)*pow(basis, _ex_1));
	}
}

int power::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<power>(other));
	const power &o = static_cast<const power &>(other);

	int cmpval = basis.compare(o.basis);
	if (cmpval)
		return cmpval;
	else
		return exponent.compare(o.exponent);
}

unsigned power::return_type() const
{
	return basis.return_type();
}

return_type_t power::return_type_tinfo() const
{
	return basis.return_type_tinfo();
}

ex power::expand(unsigned options) const
{
	if (is_a<symbol>(basis) && exponent.info(info_flags::integer)) {
		// A special case worth optimizing.
		setflag(status_flags::expanded);
		return *this;
	}

	// (x*p)^c -> x^c * p^c, if p>0
	// makes sense before expanding the basis
	if (is_exactly_a<mul>(basis) && !basis.info(info_flags::indefinite)) {
		const mul &m = ex_to<mul>(basis);
		exvector prodseq;
		epvector powseq;
		prodseq.reserve(m.seq.size() + 1);
		powseq.reserve(m.seq.size() + 1);
		bool possign = true;

		// search for positive/negative factors
		for (auto & cit : m.seq) {
			ex e=m.recombine_pair_to_ex(cit);
			if (e.info(info_flags::positive))
				prodseq.push_back(pow(e, exponent).expand(options));
			else if (e.info(info_flags::negative)) {
				prodseq.push_back(pow(-e, exponent).expand(options));
				possign = !possign;
			} else
				powseq.push_back(cit);
		}

		// take care on the numeric coefficient
		ex coeff=(possign? _ex1 : _ex_1);
		if (m.overall_coeff.info(info_flags::positive) && m.overall_coeff != _ex1)
			prodseq.push_back(pow(m.overall_coeff, exponent));
		else if (m.overall_coeff.info(info_flags::negative) && m.overall_coeff != _ex_1)
			prodseq.push_back(pow(-m.overall_coeff, exponent));
		else
			coeff *= m.overall_coeff;

		// If positive/negative factors are found, then extract them.
		// In either case we set a flag to avoid the second run on a part
		// which does not have positive/negative terms.
		if (prodseq.size() > 0) {
			ex newbasis = dynallocate<mul>(std::move(powseq), coeff);
			ex_to<basic>(newbasis).setflag(status_flags::purely_indefinite);
			return dynallocate<mul>(std::move(prodseq)) * pow(newbasis, exponent);
		} else
			ex_to<basic>(basis).setflag(status_flags::purely_indefinite);
	}

	const ex expanded_basis = basis.expand(options);
	const ex expanded_exponent = exponent.expand(options);
	
	// x^(a+b) -> x^a * x^b
	if (is_exactly_a<add>(expanded_exponent)) {
		const add &a = ex_to<add>(expanded_exponent);
		exvector distrseq;
		distrseq.reserve(a.seq.size() + 1);
		for (auto & cit : a.seq) {
			distrseq.push_back(pow(expanded_basis, a.recombine_pair_to_ex(cit)));
		}
		
		// Make sure that e.g. (x+y)^(2+a) expands the (x+y)^2 factor
		if (ex_to<numeric>(a.overall_coeff).is_integer()) {
			const numeric &num_exponent = ex_to<numeric>(a.overall_coeff);
			long int_exponent = num_exponent.to_int();
			if (int_exponent > 0 && is_exactly_a<add>(expanded_basis))
				distrseq.push_back(expand_add(ex_to<add>(expanded_basis), int_exponent, options));
			else
				distrseq.push_back(pow(expanded_basis, a.overall_coeff));
		} else
			distrseq.push_back(pow(expanded_basis, a.overall_coeff));
		
		// Make sure that e.g. (x+y)^(1+a) -> x*(x+y)^a + y*(x+y)^a
		ex r = dynallocate<mul>(distrseq);
		return r.expand(options);
	}
	
	if (!is_exactly_a<numeric>(expanded_exponent) ||
		!ex_to<numeric>(expanded_exponent).is_integer()) {
		if (are_ex_trivially_equal(basis,expanded_basis) && are_ex_trivially_equal(exponent,expanded_exponent)) {
			return this->hold();
		} else {
			return dynallocate<power>(expanded_basis, expanded_exponent).setflag(options == 0 ? status_flags::expanded : 0);
		}
	}
	
	// integer numeric exponent
	const numeric & num_exponent = ex_to<numeric>(expanded_exponent);
	long int_exponent = num_exponent.to_long();
	
	// (x+y)^n, n>0
	if (int_exponent > 0 && is_exactly_a<add>(expanded_basis))
		return expand_add(ex_to<add>(expanded_basis), int_exponent, options);
	
	// (x*y)^n -> x^n * y^n
	if (is_exactly_a<mul>(expanded_basis))
		return expand_mul(ex_to<mul>(expanded_basis), num_exponent, options, true);
	
	// cannot expand further
	if (are_ex_trivially_equal(basis,expanded_basis) && are_ex_trivially_equal(exponent,expanded_exponent))
		return this->hold();
	else
		return dynallocate<power>(expanded_basis, expanded_exponent).setflag(options == 0 ? status_flags::expanded : 0);
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

namespace {  // anonymous namespace for power::expand_add() helpers

/** Helper class to generate all bounded combinatorial partitions of an integer
 *  n with exactly m parts (including zero parts) in non-decreasing order.
 */
class partition_generator {
private:
	// Partitions n into m parts, not including zero parts.
	// (Cf. OEIS sequence A008284; implementation adapted from Jörg Arndt's
	// FXT library)
	struct mpartition2
	{
		// partition: x[1] + x[2] + ... + x[m] = n and sentinel x[0] == 0
		std::vector<int> x;
		int n;   // n>0
		int m;   // 0<m<=n
		mpartition2(unsigned n_, unsigned m_)
		  : x(m_+1), n(n_), m(m_)
		{
			for (int k=1; k<m; ++k)
				x[k] = 1;
			x[m] = n - m + 1;
		}
		bool next_partition()
		{
			int u = x[m];  // last element
			int k = m;
			int s = u;
			while (--k) {
				s += x[k];
				if (x[k] + 2 <= u)
					break;
			}
			if (k==0)
				return false;  // current is last
			int f = x[k] + 1;
			while (k < m) {
				x[k] = f;
				s -= f;
				++k;
			}
			x[m] = s;
			return true;
		}
	} mpgen;
	int m;  // number of parts 0<m<=n
	mutable std::vector<int> partition;  // current partition
public:
	partition_generator(unsigned n_, unsigned m_)
	  : mpgen(n_, 1), m(m_), partition(m_)
	{ }
	// returns current partition in non-decreasing order, padded with zeros
	const std::vector<int>& current() const
	{
		for (int i = 0; i < m - mpgen.m; ++i)
			partition[i] = 0;  // pad with zeros

		for (int i = m - mpgen.m; i < m; ++i)
			partition[i] = mpgen.x[i - m + mpgen.m + 1];

		return partition;
	}
	bool next()
	{
		if (!mpgen.next_partition()) {
			if (mpgen.m == m || mpgen.m == mpgen.n)
				return false;  // current is last
			// increment number of parts
			mpgen = mpartition2(mpgen.n, mpgen.m + 1);
		}
		return true;
	}
};

/** Helper class to generate all compositions of a partition of an integer n,
 *  starting with the compositions which has non-decreasing order.
 */
class composition_generator {
private:
	// Generates all distinct permutations of a multiset.
	// (Based on Aaron Williams' algorithm 1 from "Loopless Generation of
	// Multiset Permutations using a Constant Number of Variables by Prefix
	// Shifts." <http://webhome.csc.uvic.ca/~haron/CoolMulti.pdf>)
	struct coolmulti {
		// element of singly linked list
		struct element {
			int value;
			element* next;
			element(int val, element* n)
			  : value(val), next(n) {}
			~element()
			{   // recurses down to the end of the singly linked list
				delete next;
			}
		};
		element *head, *i, *after_i;
		// NB: Partition must be sorted in non-decreasing order.
		explicit coolmulti(const std::vector<int>& partition)
		  : head(nullptr), i(nullptr), after_i(nullptr)
		{
			for (unsigned n = 0; n < partition.size(); ++n) {
				head = new element(partition[n], head);
				if (n <= 1)
					i = head;
			}
			after_i = i->next;
		}
		~coolmulti()
		{   // deletes singly linked list
			delete head;
		}
		void next_permutation()
		{
			element *before_k;
			if (after_i->next != nullptr && i->value >= after_i->next->value)
				before_k = after_i;
			else
				before_k = i;
			element *k = before_k->next;
			before_k->next = k->next;
			k->next = head;
			if (k->value < head->value)
				i = k;
			after_i = i->next;
			head = k;
		}
		bool finished() const
		{
			return after_i->next == nullptr && after_i->value >= head->value;
		}
	} cmgen;
	bool atend;  // needed for simplifying iteration over permutations
	bool trivial;  // likewise, true if all elements are equal
	mutable std::vector<int> composition;  // current compositions
public:
	explicit composition_generator(const std::vector<int>& partition)
	  : cmgen(partition), atend(false), trivial(true), composition(partition.size())
	{
		for (unsigned i=1; i<partition.size(); ++i)
			trivial = trivial && (partition[0] == partition[i]);
	}
	const std::vector<int>& current() const
	{
		coolmulti::element* it = cmgen.head;
		size_t i = 0;
		while (it != nullptr) {
			composition[i] = it->value;
			it = it->next;
			++i;
		}
		return composition;
	}
	bool next()
	{
		// This ugly contortion is needed because the original coolmulti
		// algorithm requires code duplication of the payload procedure,
		// one before the loop and one inside it.
		if (trivial || atend)
			return false;
		cmgen.next_permutation();
		atend = cmgen.finished();
		return true;
	}
};

/** Helper function to compute the multinomial coefficient n!/(p1!*p2!*...*pk!)
 *  where n = p1+p2+...+pk, i.e. p is a partition of n.
 */
const numeric
multinomial_coefficient(const std::vector<int> & p)
{
	numeric n = 0, d = 1;
	for (auto & it : p) {
		n += numeric(it);
		d *= factorial(numeric(it));
	}
	return factorial(n) / d;
}

}  // anonymous namespace


/** expand a^n where a is an add and n is a positive integer.
 *  @see power::expand */
ex power::expand_add(const add & a, long n, unsigned options)
{
	// The special case power(+(x,...y;x),2) can be optimized better.
	if (n==2)
		return expand_add_2(a, options);

	// method:
	//
	// Consider base as the sum of all symbolic terms and the overall numeric
	// coefficient and apply the binomial theorem:
	// S = power(+(x,...,z;c),n)
	//   = power(+(+(x,...,z;0);c),n)
	//   = sum(binomial(n,k)*power(+(x,...,z;0),k)*c^(n-k), k=1..n) + c^n
	// Then, apply the multinomial theorem to expand all power(+(x,...,z;0),k):
	// The multinomial theorem is computed by an outer loop over all
	// partitions of the exponent and an inner loop over all compositions of
	// that partition. This method makes the expansion a combinatorial
	// problem and allows us to directly construct the expanded sum and also
	// to re-use the multinomial coefficients (since they depend only on the
	// partition, not on the composition).
	// 
	// multinomial power(+(x,y,z;0),3) example:
	// partition : compositions                : multinomial coefficient
	// [0,0,3]   : [3,0,0],[0,3,0],[0,0,3]     : 3!/(3!*0!*0!) = 1
	// [0,1,2]   : [2,1,0],[1,2,0],[2,0,1],... : 3!/(2!*1!*0!) = 3
	// [1,1,1]   : [1,1,1]                     : 3!/(1!*1!*1!) = 6
	//  =>  (x + y + z)^3 =
	//        x^3 + y^3 + z^3
	//      + 3*x^2*y + 3*x*y^2 + 3*y^2*z + 3*y*z^2 + 3*x*z^2 + 3*x^2*z
	//      + 6*x*y*z
	//
	// multinomial power(+(x,y,z;0),4) example:
	// partition : compositions                : multinomial coefficient
	// [0,0,4]   : [4,0,0],[0,4,0],[0,0,4]     : 4!/(4!*0!*0!) = 1
	// [0,1,3]   : [3,1,0],[1,3,0],[3,0,1],... : 4!/(3!*1!*0!) = 4
	// [0,2,2]   : [2,2,0],[2,0,2],[0,2,2]     : 4!/(2!*2!*0!) = 6
	// [1,1,2]   : [2,1,1],[1,2,1],[1,1,2]     : 4!/(2!*1!*1!) = 12
	// (no [1,1,1,1] partition since it has too many parts)
	//  =>  (x + y + z)^4 =
	//        x^4 + y^4 + z^4
	//      + 4*x^3*y + 4*x*y^3 + 4*y^3*z + 4*y*z^3 + 4*x*z^3 + 4*x^3*z
	//      + 6*x^2*y^2 + 6*y^2*z^2 + 6*x^2*z^2
	//      + 12*x^2*y*z + 12*x*y^2*z + 12*x*y*z^2
	//
	// Summary:
	// r = 0
	// for k from 0 to n:
	//     f = c^(n-k)*binomial(n,k)
	//     for p in all partitions of n with m parts (including zero parts):
	//         h = f * multinomial coefficient of p
	//         for c in all compositions of p:
	//             t = 1
	//             for e in all elements of c:
	//                 t = t * a[e]^e
	//             r = r + h*t
	// return r

	epvector result;
	// The number of terms will be the number of combinatorial compositions,
	// i.e. the number of unordered arrangements of m nonnegative integers
	// which sum up to n.  It is frequently written as C_n(m) and directly
	// related with binomial coefficients: binomial(n+m-1,m-1).
	size_t result_size = binomial(numeric(n+a.nops()-1), numeric(a.nops()-1)).to_long();
	if (!a.overall_coeff.is_zero()) {
		// the result's overall_coeff is one of the terms
		--result_size;
	}
	result.reserve(result_size);

	// Iterate over all terms in binomial expansion of
	// S = power(+(x,...,z;c),n)
	//   = sum(binomial(n,k)*power(+(x,...,z;0),k)*c^(n-k), k=1..n) + c^n
	for (int k = 1; k <= n; ++k) {
		numeric binomial_coefficient;  // binomial(n,k)*c^(n-k)
		if (a.overall_coeff.is_zero()) {
			// degenerate case with zero overall_coeff:
			// apply multinomial theorem directly to power(+(x,...z;0),n)
			binomial_coefficient = 1;
			if (k < n) {
				continue;
			}
		} else {
			binomial_coefficient = binomial(numeric(n), numeric(k)) * pow(ex_to<numeric>(a.overall_coeff), numeric(n-k));
		}

		// Multinomial expansion of power(+(x,...,z;0),k)*c^(n-k):
		// Iterate over all partitions of k with exactly as many parts as
		// there are symbolic terms in the basis (including zero parts).
		partition_generator partitions(k, a.seq.size());
		do {
			const std::vector<int>& partition = partitions.current();
			// All monomials of this partition have the same number of terms and the same coefficient.
			const unsigned msize = std::count_if(partition.begin(), partition.end(), [](int i) { return i > 0; });
			const numeric coeff = multinomial_coefficient(partition) * binomial_coefficient;

			// Iterate over all compositions of the current partition.
			composition_generator compositions(partition);
			do {
				const std::vector<int>& exponent = compositions.current();
				epvector monomial;
				monomial.reserve(msize);
				numeric factor = coeff;
				for (unsigned i = 0; i < exponent.size(); ++i) {
					const ex & r = a.seq[i].rest;
					GINAC_ASSERT(!is_exactly_a<add>(r));
					GINAC_ASSERT(!is_exactly_a<power>(r) ||
						     !is_exactly_a<numeric>(ex_to<power>(r).exponent) ||
						     !ex_to<numeric>(ex_to<power>(r).exponent).is_pos_integer() ||
						     !is_exactly_a<add>(ex_to<power>(r).basis) ||
						     !is_exactly_a<mul>(ex_to<power>(r).basis) ||
						     !is_exactly_a<power>(ex_to<power>(r).basis));
					GINAC_ASSERT(is_exactly_a<numeric>(a.seq[i].coeff));
					const numeric & c = ex_to<numeric>(a.seq[i].coeff);
					if (exponent[i] == 0) {
						// optimize away
					} else if (exponent[i] == 1) {
						// optimized
						monomial.push_back(expair(r, _ex1));
						if (c != *_num1_p)
							factor = factor.mul(c);
					} else { // general case exponent[i] > 1
						monomial.push_back(expair(r, exponent[i]));
						if (c != *_num1_p)
							factor = factor.mul(c.power(exponent[i]));
					}
				}
				result.push_back(expair(mul(std::move(monomial)).expand(options), factor));
			} while (compositions.next());
		} while (partitions.next());
	}

	GINAC_ASSERT(result.size() == result_size);
	if (a.overall_coeff.is_zero()) {
		return dynallocate<add>(std::move(result)).setflag(status_flags::expanded);
	} else {
		return dynallocate<add>(std::move(result), ex_to<numeric>(a.overall_coeff).power(n)).setflag(status_flags::expanded);
	}
}


/** Special case of power::expand_add. Expands a^2 where a is an add.
 *  @see power::expand_add */
ex power::expand_add_2(const add & a, unsigned options)
{
	epvector result;
	size_t result_size = (a.nops() * (a.nops()+1)) / 2;
	if (!a.overall_coeff.is_zero()) {
		// the result's overall_coeff is one of the terms
		--result_size;
	}
	result.reserve(result_size);

	auto last = a.seq.end();

	// power(+(x,...,z;c),2)=power(+(x,...,z;0),2)+2*c*+(x,...,z;0)+c*c
	// first part: ignore overall_coeff and expand other terms
	for (auto cit0=a.seq.begin(); cit0!=last; ++cit0) {
		const ex & r = cit0->rest;
		const ex & c = cit0->coeff;
		
		GINAC_ASSERT(!is_exactly_a<add>(r));
		GINAC_ASSERT(!is_exactly_a<power>(r) ||
		             !is_exactly_a<numeric>(ex_to<power>(r).exponent) ||
		             !ex_to<numeric>(ex_to<power>(r).exponent).is_pos_integer() ||
		             !is_exactly_a<add>(ex_to<power>(r).basis) ||
		             !is_exactly_a<mul>(ex_to<power>(r).basis) ||
		             !is_exactly_a<power>(ex_to<power>(r).basis));
		
		if (c.is_equal(_ex1)) {
			if (is_exactly_a<mul>(r)) {
				result.push_back(expair(expand_mul(ex_to<mul>(r), *_num2_p, options, true),
				                        _ex1));
			} else {
				result.push_back(expair(dynallocate<power>(r, _ex2),
				                        _ex1));
			}
		} else {
			if (is_exactly_a<mul>(r)) {
				result.push_back(expair(expand_mul(ex_to<mul>(r), *_num2_p, options, true),
				                        ex_to<numeric>(c).power_dyn(*_num2_p)));
			} else {
				result.push_back(expair(dynallocate<power>(r, _ex2),
				                        ex_to<numeric>(c).power_dyn(*_num2_p)));
			}
		}

		for (auto cit1=cit0+1; cit1!=last; ++cit1) {
			const ex & r1 = cit1->rest;
			const ex & c1 = cit1->coeff;
			result.push_back(expair(mul(r,r1).expand(options),
			                        _num2_p->mul(ex_to<numeric>(c)).mul_dyn(ex_to<numeric>(c1))));
		}
	}
	
	// second part: add terms coming from overall_coeff (if != 0)
	if (!a.overall_coeff.is_zero()) {
		for (auto & i : a.seq)
			result.push_back(a.combine_pair_with_coeff_to_pair(i, ex_to<numeric>(a.overall_coeff).mul_dyn(*_num2_p)));
	}

	GINAC_ASSERT(result.size() == result_size);

	if (a.overall_coeff.is_zero()) {
		return dynallocate<add>(std::move(result)).setflag(status_flags::expanded);
	} else {
		return dynallocate<add>(std::move(result), ex_to<numeric>(a.overall_coeff).power(2)).setflag(status_flags::expanded);
	}
}

/** Expand factors of m in m^n where m is a mul and n is an integer.
 *  @see power::expand */
ex power::expand_mul(const mul & m, const numeric & n, unsigned options, bool from_expand)
{
	GINAC_ASSERT(n.is_integer());

	if (n.is_zero()) {
		return _ex1;
	}

	// do not bother to rename indices if there are no any.
	if (!(options & expand_options::expand_rename_idx) &&
	    m.info(info_flags::has_indices))
		options |= expand_options::expand_rename_idx;
	// Leave it to multiplication since dummy indices have to be renamed
	if ((options & expand_options::expand_rename_idx) &&
	    (get_all_dummy_indices(m).size() > 0) && n.is_positive()) {
		ex result = m;
		exvector va = get_all_dummy_indices(m);
		sort(va.begin(), va.end(), ex_is_less());

		for (int i=1; i < n.to_int(); i++)
			result *= rename_dummy_indices_uniquely(va, m);
		return result;
	}

	epvector distrseq;
	distrseq.reserve(m.seq.size());
	bool need_reexpand = false;

	for (auto & cit : m.seq) {
		expair p = m.combine_pair_with_coeff_to_pair(cit, n);
		if (from_expand && is_exactly_a<add>(cit.rest) && ex_to<numeric>(p.coeff).is_pos_integer()) {
			// this happens when e.g. (a+b)^(1/2) gets squared and
			// the resulting product needs to be reexpanded
			need_reexpand = true;
		}
		distrseq.push_back(p);
	}

	const mul & result = dynallocate<mul>(std::move(distrseq), ex_to<numeric>(m.overall_coeff).power_dyn(n));
	if (need_reexpand)
		return ex(result).expand(options);
	if (from_expand)
		return result.setflag(status_flags::expanded);
	return result;
}

GINAC_BIND_UNARCHIVER(power);

} // namespace GiNaC
