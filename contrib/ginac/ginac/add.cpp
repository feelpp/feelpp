/** @file add.cpp
 *
 *  Implementation of GiNaC's sums of expressions. */

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

#include "add.h"
#include "mul.h"
#include "archive.h"
#include "operators.h"
#include "matrix.h"
#include "utils.h"
#include "clifford.h"
#include "ncmul.h"
#include "compiler.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(add, expairseq,
  print_func<print_context>(&add::do_print).
  print_func<print_latex>(&add::do_print_latex).
  print_func<print_csrc>(&add::do_print_csrc).
  print_func<print_tree>(&add::do_print_tree).
  print_func<print_python_repr>(&add::do_print_python_repr))

//////////
// default constructor
//////////

add::add()
{
}

//////////
// other constructors
//////////

// public

add::add(const ex & lh, const ex & rh)
{
	overall_coeff = _ex0;
	construct_from_2_ex(lh,rh);
	GINAC_ASSERT(is_canonical());
}

add::add(const exvector & v)
{
	overall_coeff = _ex0;
	construct_from_exvector(v);
	GINAC_ASSERT(is_canonical());
}

add::add(const epvector & v)
{
	overall_coeff = _ex0;
	construct_from_epvector(v);
	GINAC_ASSERT(is_canonical());
}

add::add(const epvector & v, const ex & oc)
{
	overall_coeff = oc;
	construct_from_epvector(v);
	GINAC_ASSERT(is_canonical());
}

add::add(epvector && vp)
{
	overall_coeff = _ex0;
	construct_from_epvector(std::move(vp));
	GINAC_ASSERT(is_canonical());
}

add::add(epvector && vp, const ex & oc)
{
	overall_coeff = oc;
	construct_from_epvector(std::move(vp));
	GINAC_ASSERT(is_canonical());
}

//////////
// archiving
//////////

GINAC_BIND_UNARCHIVER(add);

//////////
// functions overriding virtual functions from base classes
//////////

// public

void add::print_add(const print_context & c, const char *openbrace, const char *closebrace, const char *mul_sym, unsigned level) const
{
	if (precedence() <= level)
		c.s << openbrace << '(';

	numeric coeff;
	bool first = true;

	// First print the overall numeric coefficient, if present
	if (!overall_coeff.is_zero()) {
		overall_coeff.print(c, 0);
		first = false;
	}

	// Then proceed with the remaining factors
	for (auto & it : seq) {
		coeff = ex_to<numeric>(it.coeff);
		if (!first) {
			if (coeff.csgn() == -1) c.s << '-'; else c.s << '+';
		} else {
			if (coeff.csgn() == -1) c.s << '-';
			first = false;
		}
		if (!coeff.is_equal(*_num1_p) &&
		    !coeff.is_equal(*_num_1_p)) {
			if (coeff.is_rational()) {
				if (coeff.is_negative())
					(-coeff).print(c);
				else
					coeff.print(c);
			} else {
				if (coeff.csgn() == -1)
					(-coeff).print(c, precedence());
				else
					coeff.print(c, precedence());
			}
			c.s << mul_sym;
		}
		it.rest.print(c, precedence());
	}

	if (precedence() <= level)
		c.s << ')' << closebrace;
}

void add::do_print(const print_context & c, unsigned level) const
{
	print_add(c, "", "", "*", level);
}

void add::do_print_latex(const print_latex & c, unsigned level) const
{
	print_add(c, "{", "}", " ", level);
}

void add::do_print_csrc(const print_csrc & c, unsigned level) const
{
	if (precedence() <= level)
		c.s << "(";
	
	// Print arguments, separated by "+" or "-"
	char separator = ' ';
	for (auto & it : seq) {
		
		// If the coefficient is negative, separator is "-"
		if (it.coeff.is_equal(_ex_1) ||
			ex_to<numeric>(it.coeff).numer().is_equal(*_num_1_p))
			separator = '-';
		c.s << separator;
		if (it.coeff.is_equal(_ex1) || it.coeff.is_equal(_ex_1)) {
			it.rest.print(c, precedence());
		} else if (ex_to<numeric>(it.coeff).numer().is_equal(*_num1_p) ||
				 ex_to<numeric>(it.coeff).numer().is_equal(*_num_1_p))
		{
			it.rest.print(c, precedence());
			c.s << '/';
			ex_to<numeric>(it.coeff).denom().print(c, precedence());
		} else {
			it.coeff.print(c, precedence());
			c.s << '*';
			it.rest.print(c, precedence());
		}
		
		separator = '+';
	}
	
	if (!overall_coeff.is_zero()) {
		if (overall_coeff.info(info_flags::positive)
		 || is_a<print_csrc_cl_N>(c) || !overall_coeff.info(info_flags::real))  // sign inside ctor argument
			c.s << '+';
		overall_coeff.print(c, precedence());
	}
		
	if (precedence() <= level)
		c.s << ")";
}

void add::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << '(';
	op(0).print(c);
	for (size_t i=1; i<nops(); ++i) {
		c.s << ',';
		op(i).print(c);
	}
	c.s << ')';
}

bool add::info(unsigned inf) const
{
	switch (inf) {
		case info_flags::polynomial:
		case info_flags::integer_polynomial:
		case info_flags::cinteger_polynomial:
		case info_flags::rational_polynomial:
		case info_flags::real:
		case info_flags::rational:
		case info_flags::integer:
		case info_flags::crational:
		case info_flags::cinteger:
		case info_flags::positive:
		case info_flags::nonnegative:
		case info_flags::posint:
		case info_flags::nonnegint:
		case info_flags::even:
		case info_flags::crational_polynomial:
		case info_flags::rational_function: {
			for (auto & i : seq) {
				if (!(recombine_pair_to_ex(i).info(inf)))
					return false;
			}
			if (overall_coeff.is_zero() && (inf == info_flags::positive || inf == info_flags::posint))
				return true;
			return overall_coeff.info(inf);
		}
	}
	return inherited::info(inf);
}

bool add::is_polynomial(const ex & var) const
{
	for (auto & i : seq) {
		if (!i.rest.is_polynomial(var)) {
			return false;
		}
	}
	return true;
}

int add::degree(const ex & s) const
{
	int deg = std::numeric_limits<int>::min();
	if (!overall_coeff.is_zero())
		deg = 0;
	
	// Find maximum of degrees of individual terms
	for (auto & i : seq) {
		int cur_deg = i.rest.degree(s);
		if (cur_deg > deg)
			deg = cur_deg;
	}
	return deg;
}

int add::ldegree(const ex & s) const
{
	int deg = std::numeric_limits<int>::max();
	if (!overall_coeff.is_zero())
		deg = 0;
	
	// Find minimum of degrees of individual terms
	for (auto & i : seq) {
		int cur_deg = i.rest.ldegree(s);
		if (cur_deg < deg)
			deg = cur_deg;
	}
	return deg;
}

ex add::coeff(const ex & s, int n) const
{
	epvector coeffseq;
	epvector coeffseq_cliff;
	int rl = clifford_max_label(s);
	bool do_clifford = (rl != -1);
	bool nonscalar = false;

	// Calculate sum of coefficients in each term
	for (auto & i : seq) {
		ex restcoeff = i.rest.coeff(s, n);
		if (!restcoeff.is_zero()) {
			if (do_clifford) {
				if (clifford_max_label(restcoeff) == -1) {
					coeffseq_cliff.push_back(expair(ncmul(restcoeff, dirac_ONE(rl)), i.coeff));
				} else {
					coeffseq_cliff.push_back(expair(restcoeff, i.coeff));
					nonscalar = true;
				}
			}
			coeffseq.push_back(expair(restcoeff, i.coeff));
		}
	}

	return dynallocate<add>(nonscalar ? std::move(coeffseq_cliff) : std::move(coeffseq),
	                        n==0 ? overall_coeff : _ex0);
}

/** Perform automatic term rewriting rules in this class.  In the following
 *  x stands for a symbolic variables of type ex and c stands for such
 *  an expression that contain a plain number.
 *  - +(;c) -> c
 *  - +(x;0) -> x
 */
ex add::eval() const
{
	if (flags & status_flags::evaluated) {
		GINAC_ASSERT(seq.size()>0);
		GINAC_ASSERT(seq.size()>1 || !overall_coeff.is_zero());
		return *this;
	}

	const epvector evaled = evalchildren();
	if (unlikely(!evaled.empty())) {
		// start over evaluating a new object
		return dynallocate<add>(std::move(evaled), overall_coeff);
	}

#ifdef DO_GINAC_ASSERT
	for (auto & i : seq) {
		GINAC_ASSERT(!is_exactly_a<add>(i.rest));
	}
#endif // def DO_GINAC_ASSERT

	size_t seq_size = seq.size();
	if (seq_size == 0) {
		// +(;c) -> c
		return overall_coeff;
	} else if (seq_size == 1 && overall_coeff.is_zero()) {
		// +(x;0) -> x
		return recombine_pair_to_ex(*(seq.begin()));
	} else if (!overall_coeff.is_zero() && seq[0].rest.return_type() != return_types::commutative) {
		throw (std::logic_error("add::eval(): sum of non-commutative objects has non-zero numeric term"));
	}

	return this->hold();
}

ex add::evalm() const
{
	// Evaluate children first and add up all matrices. Stop if there's one
	// term that is not a matrix.
	epvector s;
	s.reserve(seq.size());

	bool all_matrices = true;
	bool first_term = true;
	matrix sum;

	for (auto & it : seq) {
		const ex &m = recombine_pair_to_ex(it).evalm();
		s.push_back(split_ex_to_pair(m));
		if (is_a<matrix>(m)) {
			if (first_term) {
				sum = ex_to<matrix>(m);
				first_term = false;
			} else
				sum = sum.add(ex_to<matrix>(m));
		} else
			all_matrices = false;
	}

	if (all_matrices)
		return sum + overall_coeff;
	else
		return dynallocate<add>(std::move(s), overall_coeff);
}

ex add::conjugate() const
{
	std::unique_ptr<exvector> v(nullptr);
	for (size_t i=0; i<nops(); ++i) {
		if (v) {
			v->push_back(op(i).conjugate());
			continue;
		}
		ex term = op(i);
		ex ccterm = term.conjugate();
		if (are_ex_trivially_equal(term, ccterm))
			continue;
		v.reset(new exvector);
		v->reserve(nops());
		for (size_t j=0; j<i; ++j)
			v->push_back(op(j));
		v->push_back(ccterm);
	}
	if (v) {
		return add(std::move(*v));
	}
	return *this;
}

ex add::real_part() const
{
	epvector v;
	v.reserve(seq.size());
	for (auto & it : seq)
		if (it.coeff.info(info_flags::real)) {
			ex rp = it.rest.real_part();
			if (!rp.is_zero())
				v.push_back(expair(rp, it.coeff));
		} else {
			ex rp = recombine_pair_to_ex(it).real_part();
			if (!rp.is_zero())
				v.push_back(split_ex_to_pair(rp));
		}
	return dynallocate<add>(std::move(v), overall_coeff.real_part());
}

ex add::imag_part() const
{
	epvector v;
	v.reserve(seq.size());
	for (auto & it : seq)
		if (it.coeff.info(info_flags::real)) {
			ex ip = it.rest.imag_part();
			if (!ip.is_zero())
				v.push_back(expair(ip, it.coeff));
		} else {
			ex ip = recombine_pair_to_ex(it).imag_part();
			if (!ip.is_zero())
				v.push_back(split_ex_to_pair(ip));
		}
	return dynallocate<add>(std::move(v), overall_coeff.imag_part());
}

ex add::eval_ncmul(const exvector & v) const
{
	if (seq.empty())
		return inherited::eval_ncmul(v);
	else
		return seq.begin()->rest.eval_ncmul(v);
}    

// protected

/** Implementation of ex::diff() for a sum. It differentiates each term.
 *  @see ex::diff */
ex add::derivative(const symbol & y) const
{
	epvector s;
	s.reserve(seq.size());
	
	// Only differentiate the "rest" parts of the expairs. This is faster
	// than the default implementation in basic::derivative() although
	// if performs the same function (differentiate each term).
	for (auto & it : seq)
		s.push_back(expair(it.rest.diff(y), it.coeff));

	return dynallocate<add>(std::move(s));
}

int add::compare_same_type(const basic & other) const
{
	return inherited::compare_same_type(other);
}

unsigned add::return_type() const
{
	if (seq.empty())
		return return_types::commutative;
	else
		return seq.begin()->rest.return_type();
}

return_type_t add::return_type_tinfo() const
{
	if (seq.empty())
		return make_return_type_t<add>();
	else
		return seq.begin()->rest.return_type_tinfo();
}

// Note: do_index_renaming is ignored because it makes no sense for an add.
ex add::thisexpairseq(const epvector & v, const ex & oc, bool do_index_renaming) const
{
	return dynallocate<add>(v, oc);
}

// Note: do_index_renaming is ignored because it makes no sense for an add.
ex add::thisexpairseq(epvector && vp, const ex & oc, bool do_index_renaming) const
{
	return dynallocate<add>(std::move(vp), oc);
}

expair add::split_ex_to_pair(const ex & e) const
{
	if (is_exactly_a<mul>(e)) {
		const mul &mulref(ex_to<mul>(e));
		const ex &numfactor = mulref.overall_coeff;
		if (numfactor.is_equal(_ex1))
			return expair(e, _ex1);
		mul & mulcopy = dynallocate<mul>(mulref);
		mulcopy.overall_coeff = _ex1;
		mulcopy.clearflag(status_flags::evaluated | status_flags::hash_calculated);
		return expair(mulcopy, numfactor);
	}
	return expair(e,_ex1);
}

expair add::combine_ex_with_coeff_to_pair(const ex & e,
                                          const ex & c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	if (is_exactly_a<mul>(e)) {
		const mul &mulref(ex_to<mul>(e));
		const ex &numfactor = mulref.overall_coeff;
		if (likely(numfactor.is_equal(_ex1)))
			return expair(e, c);
		mul & mulcopy = dynallocate<mul>(mulref);
		mulcopy.overall_coeff = _ex1;
		mulcopy.clearflag(status_flags::evaluated | status_flags::hash_calculated);
		if (c.is_equal(_ex1))
			return expair(mulcopy, numfactor);
		else
			return expair(mulcopy, ex_to<numeric>(numfactor).mul_dyn(ex_to<numeric>(c)));
	} else if (is_exactly_a<numeric>(e)) {
		if (c.is_equal(_ex1))
			return expair(e, _ex1);
		if (e.is_equal(_ex1))
			return expair(c, _ex1);
		return expair(ex_to<numeric>(e).mul_dyn(ex_to<numeric>(c)), _ex1);
	}
	return expair(e, c);
}

expair add::combine_pair_with_coeff_to_pair(const expair & p,
                                            const ex & c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c));

	if (is_exactly_a<numeric>(p.rest)) {
		GINAC_ASSERT(ex_to<numeric>(p.coeff).is_equal(*_num1_p)); // should be normalized
		return expair(ex_to<numeric>(p.rest).mul_dyn(ex_to<numeric>(c)),_ex1);
	}

	return expair(p.rest,ex_to<numeric>(p.coeff).mul_dyn(ex_to<numeric>(c)));
}

ex add::recombine_pair_to_ex(const expair & p) const
{
	if (ex_to<numeric>(p.coeff).is_equal(*_num1_p))
		return p.rest;
	else
		return dynallocate<mul>(p.rest, p.coeff);
}

ex add::expand(unsigned options) const
{
	epvector expanded = expandchildren(options);
	if (expanded.empty())
		return (options == 0) ? setflag(status_flags::expanded) : *this;

	return dynallocate<add>(std::move(expanded), overall_coeff).setflag(options == 0 ? status_flags::expanded : 0);
}

} // namespace GiNaC
extern template GiNaC::registered_class_info GiNaC::container<std::vector>::reg_info;
