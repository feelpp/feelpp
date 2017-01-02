/** @file pseries.cpp
 *
 *  Implementation of class for extended truncated power series and
 *  methods for series expansion. */

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

#include "pseries.h"
#include "add.h"
#include "inifcns.h" // for Order function
#include "lst.h"
#include "mul.h"
#include "power.h"
#include "relational.h"
#include "operators.h"
#include "symbol.h"
#include "integral.h"
#include "archive.h"
#include "utils.h"

#include <limits>
#include <numeric>
#include <stdexcept>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(pseries, basic,
  print_func<print_context>(&pseries::do_print).
  print_func<print_latex>(&pseries::do_print_latex).
  print_func<print_tree>(&pseries::do_print_tree).
  print_func<print_python>(&pseries::do_print_python).
  print_func<print_python_repr>(&pseries::do_print_python_repr))


/*
 *  Default constructor
 */

pseries::pseries() { }


/*
 *  Other ctors
 */

/** Construct pseries from a vector of coefficients and powers.
 *  expair.rest holds the coefficient, expair.coeff holds the power.
 *  The powers must be integers (positive or negative) and in ascending order;
 *  the last coefficient can be Order(_ex1) to represent a truncated,
 *  non-terminating series.
 *
 *  @param rel_  expansion variable and point (must hold a relational)
 *  @param ops_  vector of {coefficient, power} pairs (coefficient must not be zero)
 *  @return newly constructed pseries */
pseries::pseries(const ex &rel_, const epvector &ops_)
  : seq(ops_)
{
#ifdef DO_GINAC_ASSERT
	auto i = seq.begin();
	while (i != seq.end()) {
		auto ip1 = i+1;
		if (ip1 != seq.end())
			GINAC_ASSERT(!is_order_function(i->rest));
		else
			break;
		GINAC_ASSERT(is_a<numeric>(i->coeff));
		GINAC_ASSERT(ex_to<numeric>(i->coeff) < ex_to<numeric>(ip1->coeff));
		++i;
	}
#endif // def DO_GINAC_ASSERT
	GINAC_ASSERT(is_a<relational>(rel_));
	GINAC_ASSERT(is_a<symbol>(rel_.lhs()));
	point = rel_.rhs();
	var = rel_.lhs();
}
pseries::pseries(const ex &rel_, epvector &&ops_)
  : seq(std::move(ops_))
{
#ifdef DO_GINAC_ASSERT
	auto i = seq.begin();
	while (i != seq.end()) {
		auto ip1 = i+1;
		if (ip1 != seq.end())
			GINAC_ASSERT(!is_order_function(i->rest));
		else
			break;
		GINAC_ASSERT(is_a<numeric>(i->coeff));
		GINAC_ASSERT(ex_to<numeric>(i->coeff) < ex_to<numeric>(ip1->coeff));
		++i;
	}
#endif // def DO_GINAC_ASSERT
	GINAC_ASSERT(is_a<relational>(rel_));
	GINAC_ASSERT(is_a<symbol>(rel_.lhs()));
	point = rel_.rhs();
	var = rel_.lhs();
}


/*
 *  Archiving
 */

void pseries::read_archive(const archive_node &n, lst &sym_lst) 
{
	inherited::read_archive(n, sym_lst);
	auto first = n.find_first("coeff");
	auto last = n.find_last("power");
	++last;
	seq.reserve((last-first)/2);

	for (auto loc = first; loc < last;) {
		ex rest;
		ex coeff;
		n.find_ex_by_loc(loc++, rest, sym_lst);
		n.find_ex_by_loc(loc++, coeff, sym_lst);
		seq.push_back(expair(rest, coeff));
	}

	n.find_ex("var", var, sym_lst);
	n.find_ex("point", point, sym_lst);
}

void pseries::archive(archive_node &n) const
{
	inherited::archive(n);
	for (auto & it : seq) {
		n.add_ex("coeff", it.rest);
		n.add_ex("power", it.coeff);
	}
	n.add_ex("var", var);
	n.add_ex("point", point);
}


//////////
// functions overriding virtual functions from base classes
//////////

void pseries::print_series(const print_context & c, const char *openbrace, const char *closebrace, const char *mul_sym, const char *pow_sym, unsigned level) const
{
	if (precedence() <= level)
		c.s << '(';
		
	// objects of type pseries must not have any zero entries, so the
	// trivial (zero) pseries needs a special treatment here:
	if (seq.empty())
		c.s << '0';

	auto i = seq.begin(), end = seq.end();
	while (i != end) {

		// print a sign, if needed
		if (i != seq.begin())
			c.s << '+';

		if (!is_order_function(i->rest)) {

			// print 'rest', i.e. the expansion coefficient
			if (i->rest.info(info_flags::numeric) &&
				i->rest.info(info_flags::positive)) {
				i->rest.print(c);
			} else {
				c.s << openbrace << '(';
				i->rest.print(c);
				c.s << ')' << closebrace;
			}

			// print 'coeff', something like (x-1)^42
			if (!i->coeff.is_zero()) {
				c.s << mul_sym;
				if (!point.is_zero()) {
					c.s << openbrace << '(';
					(var-point).print(c);
					c.s << ')' << closebrace;
				} else
					var.print(c);
				if (i->coeff.compare(_ex1)) {
					c.s << pow_sym;
					c.s << openbrace;
					if (i->coeff.info(info_flags::negative)) {
						c.s << '(';
						i->coeff.print(c);
						c.s << ')';
					} else
						i->coeff.print(c);
					c.s << closebrace;
				}
			}
		} else
			Order(pow(var - point, i->coeff)).print(c);
		++i;
	}

	if (precedence() <= level)
		c.s << ')';
}

void pseries::do_print(const print_context & c, unsigned level) const
{
	print_series(c, "", "", "*", "^", level);
}

void pseries::do_print_latex(const print_latex & c, unsigned level) const
{
	print_series(c, "{", "}", " ", "^", level);
}

void pseries::do_print_python(const print_python & c, unsigned level) const
{
	print_series(c, "", "", "*", "**", level);
}

void pseries::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << std::endl;
	size_t num = seq.size();
	for (size_t i=0; i<num; ++i) {
		seq[i].rest.print(c, level + c.delta_indent);
		seq[i].coeff.print(c, level + c.delta_indent);
		c.s << std::string(level + c.delta_indent, ' ') << "-----" << std::endl;
	}
	var.print(c, level + c.delta_indent);
	point.print(c, level + c.delta_indent);
}

void pseries::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << "(relational(";
	var.print(c);
	c.s << ',';
	point.print(c);
	c.s << "),[";
	size_t num = seq.size();
	for (size_t i=0; i<num; ++i) {
		if (i)
			c.s << ',';
		c.s << '(';
		seq[i].rest.print(c);
		c.s << ',';
		seq[i].coeff.print(c);
		c.s << ')';
	}
	c.s << "])";
}

int pseries::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<pseries>(other));
	const pseries &o = static_cast<const pseries &>(other);
	
	// first compare the lengths of the series...
	if (seq.size()>o.seq.size())
		return 1;
	if (seq.size()<o.seq.size())
		return -1;
	
	// ...then the expansion point...
	int cmpval = var.compare(o.var);
	if (cmpval)
		return cmpval;
	cmpval = point.compare(o.point);
	if (cmpval)
		return cmpval;
	
	// ...and if that failed the individual elements
	auto it = seq.begin(), o_it = o.seq.begin();
	while (it!=seq.end() && o_it!=o.seq.end()) {
		cmpval = it->compare(*o_it);
		if (cmpval)
			return cmpval;
		++it;
		++o_it;
	}

	// so they are equal.
	return 0;
}

/** Return the number of operands including a possible order term. */
size_t pseries::nops() const
{
	return seq.size();
}

/** Return the ith term in the series when represented as a sum. */
ex pseries::op(size_t i) const
{
	if (i >= seq.size())
		throw (std::out_of_range("op() out of range"));

	if (is_order_function(seq[i].rest))
		return Order(pow(var-point, seq[i].coeff));
	return seq[i].rest * pow(var - point, seq[i].coeff);
}

/** Return degree of highest power of the series.  This is usually the exponent
 *  of the Order term.  If s is not the expansion variable of the series, the
 *  series is examined termwise. */
int pseries::degree(const ex &s) const
{
	if (seq.empty())
		return 0;

	if (var.is_equal(s))
		// Return last/greatest exponent
		return ex_to<numeric>((seq.end()-1)->coeff).to_int();

	int max_pow = std::numeric_limits<int>::min();
	for (auto & it : seq)
		max_pow = std::max(max_pow, it.rest.degree(s));
	return max_pow;
}

/** Return degree of lowest power of the series.  This is usually the exponent
 *  of the leading term.  If s is not the expansion variable of the series, the
 *  series is examined termwise.  If s is the expansion variable but the
 *  expansion point is not zero the series is not expanded to find the degree.
 *  I.e.: (1-x) + (1-x)^2 + Order((1-x)^3) has ldegree(x) 1, not 0. */
int pseries::ldegree(const ex &s) const
{
	if (seq.empty())
		return 0;

	if (var.is_equal(s))
		// Return first/smallest exponent
		return ex_to<numeric>((seq.begin())->coeff).to_int();

	int min_pow = std::numeric_limits<int>::max();
	for (auto & it : seq)
		min_pow = std::min(min_pow, it.rest.degree(s));
	return min_pow;
}

/** Return coefficient of degree n in power series if s is the expansion
 *  variable.  If the expansion point is nonzero, by definition the n=1
 *  coefficient in s of a+b*(s-z)+c*(s-z)^2+Order((s-z)^3) is b (assuming
 *  the expansion took place in the s in the first place).
 *  If s is not the expansion variable, an attempt is made to convert the
 *  series to a polynomial and return the corresponding coefficient from
 *  there. */
ex pseries::coeff(const ex &s, int n) const
{
	if (var.is_equal(s)) {
		if (seq.empty())
			return _ex0;
		
		// Binary search in sequence for given power
		numeric looking_for = numeric(n);
		int lo = 0, hi = seq.size() - 1;
		while (lo <= hi) {
			int mid = (lo + hi) / 2;
			GINAC_ASSERT(is_exactly_a<numeric>(seq[mid].coeff));
			int cmp = ex_to<numeric>(seq[mid].coeff).compare(looking_for);
			switch (cmp) {
				case -1:
					lo = mid + 1;
					break;
				case 0:
					return seq[mid].rest;
				case 1:
					hi = mid - 1;
					break;
				default:
					throw(std::logic_error("pseries::coeff: compare() didn't return -1, 0 or 1"));
			}
		}
		return _ex0;
	} else
		return convert_to_poly().coeff(s, n);
}

/** Does nothing. */
ex pseries::collect(const ex &s, bool distributed) const
{
	return *this;
}

/** Perform coefficient-wise automatic term rewriting rules in this class. */
ex pseries::eval() const
{
	if (flags & status_flags::evaluated) {
		return *this;
	}

	// Construct a new series with evaluated coefficients
	epvector new_seq;
	new_seq.reserve(seq.size());
	for (auto & it : seq)
		new_seq.push_back(expair(it.rest, it.coeff));

	return dynallocate<pseries>(relational(var,point), std::move(new_seq)).setflag(status_flags::evaluated);
}

/** Evaluate coefficients numerically. */
ex pseries::evalf() const
{
	// Construct a new series with evaluated coefficients
	epvector new_seq;
	new_seq.reserve(seq.size());
	for (auto & it : seq)
		new_seq.push_back(expair(it.rest, it.coeff));

	return dynallocate<pseries>(relational(var,point), std::move(new_seq)).setflag(status_flags::evaluated);
}

ex pseries::conjugate() const
{
	if(!var.info(info_flags::real))
		return conjugate_function(*this).hold();

	std::unique_ptr<epvector> newseq(conjugateepvector(seq));
	ex newpoint = point.conjugate();

	if (!newseq && are_ex_trivially_equal(point, newpoint)) {
		return *this;
	}

	return dynallocate<pseries>(var==newpoint, newseq ? std::move(*newseq) : seq);
}

ex pseries::real_part() const
{
	if(!var.info(info_flags::real))
		return real_part_function(*this).hold();
	ex newpoint = point.real_part();
	if(newpoint != point)
		return real_part_function(*this).hold();

	epvector v;
	v.reserve(seq.size());
	for (auto & it : seq)
		v.push_back(expair((it.rest).real_part(), it.coeff));
	return dynallocate<pseries>(var==point, std::move(v));
}

ex pseries::imag_part() const
{
	if(!var.info(info_flags::real))
		return imag_part_function(*this).hold();
	ex newpoint = point.real_part();
	if(newpoint != point)
		return imag_part_function(*this).hold();

	epvector v;
	v.reserve(seq.size());
	for (auto & it : seq)
		v.push_back(expair((it.rest).imag_part(), it.coeff));
	return dynallocate<pseries>(var==point, std::move(v));
}

ex pseries::eval_integ() const
{
	std::unique_ptr<epvector> newseq(nullptr);
	for (auto i=seq.begin(); i!=seq.end(); ++i) {
		if (newseq) {
			newseq->push_back(expair(i->rest.eval_integ(), i->coeff));
			continue;
		}
		ex newterm = i->rest.eval_integ();
		if (!are_ex_trivially_equal(newterm, i->rest)) {
			newseq.reset(new epvector);
			newseq->reserve(seq.size());
			for (auto j=seq.begin(); j!=i; ++j)
				newseq->push_back(*j);
			newseq->push_back(expair(newterm, i->coeff));
		}
	}

	ex newpoint = point.eval_integ();
	if (newseq || !are_ex_trivially_equal(newpoint, point))
		return dynallocate<pseries>(var==newpoint, std::move(*newseq));
	return *this;
}

ex pseries::evalm() const
{
	// evalm each coefficient
	epvector newseq;
	bool something_changed = false;
	for (auto i=seq.begin(); i!=seq.end(); ++i) {
		if (something_changed) {
			ex newcoeff = i->rest.evalm();
			if (!newcoeff.is_zero())
				newseq.push_back(expair(newcoeff, i->coeff));
		} else {
			ex newcoeff = i->rest.evalm();
			if (!are_ex_trivially_equal(newcoeff, i->rest)) {
				something_changed = true;
				newseq.reserve(seq.size());
				std::copy(seq.begin(), i, std::back_inserter<epvector>(newseq));
				if (!newcoeff.is_zero())
					newseq.push_back(expair(newcoeff, i->coeff));
			}
		}
	}
	if (something_changed)
		return dynallocate<pseries>(var==point, std::move(newseq));
	else
		return *this;
}

ex pseries::subs(const exmap & m, unsigned options) const
{
	// If expansion variable is being substituted, convert the series to a
	// polynomial and do the substitution there because the result might
	// no longer be a power series
	if (m.find(var) != m.end())
		return convert_to_poly(true).subs(m, options);
	
	// Otherwise construct a new series with substituted coefficients and
	// expansion point
	epvector newseq;
	newseq.reserve(seq.size());
	for (auto & it : seq)
		newseq.push_back(expair(it.rest.subs(m, options), it.coeff));
	return dynallocate<pseries>(relational(var,point.subs(m, options)), std::move(newseq));
}

/** Implementation of ex::expand() for a power series.  It expands all the
 *  terms individually and returns the resulting series as a new pseries. */
ex pseries::expand(unsigned options) const
{
	epvector newseq;
	for (auto & it : seq) {
		ex restexp = it.rest.expand();
		if (!restexp.is_zero())
			newseq.push_back(expair(restexp, it.coeff));
	}
	return dynallocate<pseries>(relational(var,point), std::move(newseq)).setflag(options == 0 ? status_flags::expanded : 0);
}

/** Implementation of ex::diff() for a power series.
 *  @see ex::diff */
ex pseries::derivative(const symbol & s) const
{
	epvector new_seq;

	if (s == var) {
		
		// FIXME: coeff might depend on var
		for (auto & it : seq) {
			if (is_order_function(it.rest)) {
				new_seq.push_back(expair(it.rest, it.coeff - 1));
			} else {
				ex c = it.rest * it.coeff;
				if (!c.is_zero())
					new_seq.push_back(expair(c, it.coeff - 1));
			}
		}

	} else {

		for (auto & it : seq) {
			if (is_order_function(it.rest)) {
				new_seq.push_back(it);
			} else {
				ex c = it.rest.diff(s);
				if (!c.is_zero())
					new_seq.push_back(expair(c, it.coeff));
			}
		}
	}

	return pseries(relational(var,point), std::move(new_seq));
}

ex pseries::convert_to_poly(bool no_order) const
{
	ex e;
	for (auto & it : seq) {
		if (is_order_function(it.rest)) {
			if (!no_order)
				e += Order(pow(var - point, it.coeff));
		} else
			e += it.rest * pow(var - point, it.coeff);
	}
	return e;
}

bool pseries::is_terminating() const
{
	return seq.empty() || !is_order_function((seq.end()-1)->rest);
}

ex pseries::coeffop(size_t i) const
{
	if (i >= nops())
		throw (std::out_of_range("coeffop() out of range"));
	return seq[i].rest;
}

ex pseries::exponop(size_t i) const
{
	if (i >= nops())
		throw (std::out_of_range("exponop() out of range"));
	return seq[i].coeff;
}


/*
 *  Implementations of series expansion
 */

/** Default implementation of ex::series(). This performs Taylor expansion.
 *  @see ex::series */
ex basic::series(const relational & r, int order, unsigned options) const
{
	epvector seq;
	const symbol &s = ex_to<symbol>(r.lhs());

	// default for order-values that make no sense for Taylor expansion
	if ((order <= 0) && this->has(s)) {
		seq.push_back(expair(Order(_ex1), order));
		return pseries(r, std::move(seq));
	}

	// do Taylor expansion
	numeric fac = 1;
	ex deriv = *this;
	ex coeff = deriv.subs(r, subs_options::no_pattern);

	if (!coeff.is_zero()) {
		seq.push_back(expair(coeff, _ex0));
	}

	int n;
	for (n=1; n<order; ++n) {
		fac = fac.div(n);
		// We need to test for zero in order to see if the series terminates.
		// The problem is that there is no such thing as a perfect test for
		// zero.  Expanding the term occasionally helps a little...
		deriv = deriv.diff(s).expand();
		if (deriv.is_zero())  // Series terminates
			return pseries(r, std::move(seq));

		coeff = deriv.subs(r, subs_options::no_pattern);
		if (!coeff.is_zero())
			seq.push_back(expair(fac * coeff, n));
	}
	
	// Higher-order terms, if present
	deriv = deriv.diff(s);
	if (!deriv.expand().is_zero())
		seq.push_back(expair(Order(_ex1), n));
	return pseries(r, std::move(seq));
}


/** Implementation of ex::series() for symbols.
 *  @see ex::series */
ex symbol::series(const relational & r, int order, unsigned options) const
{
	epvector seq;
	const ex point = r.rhs();
	GINAC_ASSERT(is_a<symbol>(r.lhs()));

	if (this->is_equal_same_type(ex_to<symbol>(r.lhs()))) {
		if (order > 0 && !point.is_zero())
			seq.push_back(expair(point, _ex0));
		if (order > 1)
			seq.push_back(expair(_ex1, _ex1));
		else
			seq.push_back(expair(Order(_ex1), numeric(order)));
	} else
		seq.push_back(expair(*this, _ex0));
	return pseries(r, std::move(seq));
}


/** Add one series object to another, producing a pseries object that
 *  represents the sum.
 *
 *  @param other  pseries object to add with
 *  @return the sum as a pseries */
ex pseries::add_series(const pseries &other) const
{
	// Adding two series with different variables or expansion points
	// results in an empty (constant) series 
	if (!is_compatible_to(other)) {
		epvector nul { expair(Order(_ex1), _ex0) };
		return pseries(relational(var,point), std::move(nul));
	}
	
	// Series addition
	epvector new_seq;
	auto a = seq.begin(), a_end = seq.end();
	auto b = other.seq.begin(), b_end = other.seq.end();
	int pow_a = std::numeric_limits<int>::max(), pow_b = std::numeric_limits<int>::max();
	for (;;) {
		// If a is empty, fill up with elements from b and stop
		if (a == a_end) {
			while (b != b_end) {
				new_seq.push_back(*b);
				++b;
			}
			break;
		} else
			pow_a = ex_to<numeric>((*a).coeff).to_int();
		
		// If b is empty, fill up with elements from a and stop
		if (b == b_end) {
			while (a != a_end) {
				new_seq.push_back(*a);
				++a;
			}
			break;
		} else
			pow_b = ex_to<numeric>((*b).coeff).to_int();
		
		// a and b are non-empty, compare powers
		if (pow_a < pow_b) {
			// a has lesser power, get coefficient from a
			new_seq.push_back(*a);
			if (is_order_function((*a).rest))
				break;
			++a;
		} else if (pow_b < pow_a) {
			// b has lesser power, get coefficient from b
			new_seq.push_back(*b);
			if (is_order_function((*b).rest))
				break;
			++b;
		} else {
			// Add coefficient of a and b
			if (is_order_function((*a).rest) || is_order_function((*b).rest)) {
				new_seq.push_back(expair(Order(_ex1), (*a).coeff));
				break;  // Order term ends the sequence
			} else {
				ex sum = (*a).rest + (*b).rest;
				if (!(sum.is_zero()))
					new_seq.push_back(expair(sum, numeric(pow_a)));
				++a;
				++b;
			}
		}
	}
	return pseries(relational(var,point), std::move(new_seq));
}


/** Implementation of ex::series() for sums. This performs series addition when
 *  adding pseries objects.
 *  @see ex::series */
ex add::series(const relational & r, int order, unsigned options) const
{
	ex acc; // Series accumulator
	
	// Get first term from overall_coeff
	acc = overall_coeff.series(r, order, options);
	
	// Add remaining terms
	for (auto & it : seq) {
		ex op;
		if (is_exactly_a<pseries>(it.rest))
			op = it.rest;
		else
			op = it.rest.series(r, order, options);
		if (!it.coeff.is_equal(_ex1))
			op = ex_to<pseries>(op).mul_const(ex_to<numeric>(it.coeff));
		
		// Series addition
		acc = ex_to<pseries>(acc).add_series(ex_to<pseries>(op));
	}
	return acc;
}


/** Multiply a pseries object with a numeric constant, producing a pseries
 *  object that represents the product.
 *
 *  @param other  constant to multiply with
 *  @return the product as a pseries */
ex pseries::mul_const(const numeric &other) const
{
	epvector new_seq;
	new_seq.reserve(seq.size());
	
	for (auto & it : seq) {
		if (!is_order_function(it.rest))
			new_seq.push_back(expair(it.rest * other, it.coeff));
		else
			new_seq.push_back(it);
	}
	return pseries(relational(var,point), std::move(new_seq));
}


/** Multiply one pseries object to another, producing a pseries object that
 *  represents the product.
 *
 *  @param other  pseries object to multiply with
 *  @return the product as a pseries */
ex pseries::mul_series(const pseries &other) const
{
	// Multiplying two series with different variables or expansion points
	// results in an empty (constant) series 
	if (!is_compatible_to(other)) {
		epvector nul { expair(Order(_ex1), _ex0) };
		return pseries(relational(var,point), std::move(nul));
	}

	if (seq.empty() || other.seq.empty()) {
		return dynallocate<pseries>(var==point, epvector());
	}
	
	// Series multiplication
	epvector new_seq;
	const int a_max = degree(var);
	const int b_max = other.degree(var);
	const int a_min = ldegree(var);
	const int b_min = other.ldegree(var);
	const int cdeg_min = a_min + b_min;
	int cdeg_max = a_max + b_max;
	
	int higher_order_a = std::numeric_limits<int>::max();
	int higher_order_b = std::numeric_limits<int>::max();
	if (is_order_function(coeff(var, a_max)))
		higher_order_a = a_max + b_min;
	if (is_order_function(other.coeff(var, b_max)))
		higher_order_b = b_max + a_min;
	const int higher_order_c = std::min(higher_order_a, higher_order_b);
	if (cdeg_max >= higher_order_c)
		cdeg_max = higher_order_c - 1;

	std::map<int, ex> rest_map_a, rest_map_b;
	for (const auto& it : seq)
		rest_map_a[ex_to<numeric>(it.coeff).to_int()] = it.rest;

	if (other.var.is_equal(var))
		for (const auto& it : other.seq)
			rest_map_b[ex_to<numeric>(it.coeff).to_int()] = it.rest;

	for (int cdeg=cdeg_min; cdeg<=cdeg_max; ++cdeg) {
		ex co = _ex0;
		// c(i)=a(0)b(i)+...+a(i)b(0)
		for (int i=a_min; cdeg-i>=b_min; ++i) {
			const auto& ita = rest_map_a.find(i);
			if (ita == rest_map_a.end())
				continue;
			const auto& itb = rest_map_b.find(cdeg-i);
			if (itb == rest_map_b.end())
				continue;
			if (!is_order_function(ita->second) && !is_order_function(itb->second))
				co += ita->second * itb->second;
		}
		if (!co.is_zero())
			new_seq.push_back(expair(co, numeric(cdeg)));
	}
	if (higher_order_c < std::numeric_limits<int>::max())
		new_seq.push_back(expair(Order(_ex1), numeric(higher_order_c)));
	return pseries(relational(var, point), std::move(new_seq));
}


/** Implementation of ex::series() for product. This performs series
 *  multiplication when multiplying series.
 *  @see ex::series */
ex mul::series(const relational & r, int order, unsigned options) const
{
	pseries acc; // Series accumulator

	GINAC_ASSERT(is_a<symbol>(r.lhs()));
	const ex& sym = r.lhs();
		
	// holds ldegrees of the series of individual factors
	std::vector<int> ldegrees;
	std::vector<bool> ldegree_redo;

	// find minimal degrees
	// first round: obtain a bound up to which minimal degrees have to be
	// considered
	for (auto & it : seq) {

		ex expon = it.coeff;
		int factor = 1;
		ex buf;
		if (expon.info(info_flags::integer)) {
			buf = it.rest;
			factor = ex_to<numeric>(expon).to_int();
		} else {
			buf = recombine_pair_to_ex(it);
		}

		int real_ldegree = 0;
		bool flag_redo = false;
		try {
			real_ldegree = buf.expand().ldegree(sym-r.rhs());
		} catch (std::runtime_error) {}

		if (real_ldegree == 0) {
			if ( factor < 0 ) {
				// This case must terminate, otherwise we would have division by
				// zero.
				int orderloop = 0;
				do {
					orderloop++;
					real_ldegree = buf.series(r, orderloop, options).ldegree(sym);
				} while (real_ldegree == orderloop);
			} else {
				// Here it is possible that buf does not have a ldegree, therefore
				// check only if ldegree is negative, otherwise reconsider the case
				// in the second round.
				real_ldegree = buf.series(r, 0, options).ldegree(sym);
				if (real_ldegree == 0)
					flag_redo = true;
			}
		}

		ldegrees.push_back(factor * real_ldegree);
		ldegree_redo.push_back(flag_redo);
	}

	int degbound = order-std::accumulate(ldegrees.begin(), ldegrees.end(), 0);
	// Second round: determine the remaining positive ldegrees by the series
	// method.
	// here we can ignore ldegrees larger than degbound
	size_t j = 0;
	for (auto & it : seq) {
		if ( ldegree_redo[j] ) {
			ex expon = it.coeff;
			int factor = 1;
			ex buf;
			if (expon.info(info_flags::integer)) {
				buf = it.rest;
				factor = ex_to<numeric>(expon).to_int();
			} else {
				buf = recombine_pair_to_ex(it);
			}
			int real_ldegree = 0;
			int orderloop = 0;
			do {
				orderloop++;
				real_ldegree = buf.series(r, orderloop, options).ldegree(sym);
			} while ((real_ldegree == orderloop)
			      && (factor*real_ldegree < degbound));
			ldegrees[j] = factor * real_ldegree;
			degbound -= factor * real_ldegree;
		}
		j++;
	}

	int degsum = std::accumulate(ldegrees.begin(), ldegrees.end(), 0);

	if (degsum >= order) {
		epvector epv { expair(Order(_ex1), order) };
		return dynallocate<pseries>(r, std::move(epv));
	}

	// Multiply with remaining terms
	auto itd = ldegrees.begin();
	for (auto it=seq.begin(), itend=seq.end(); it!=itend; ++it, ++itd) {

		// do series expansion with adjusted order
		ex op = recombine_pair_to_ex(*it).series(r, order-degsum+(*itd), options);

		// Series multiplication
		if (it == seq.begin())
			acc = ex_to<pseries>(op);
		else
			acc = ex_to<pseries>(acc.mul_series(ex_to<pseries>(op)));
	}

	return acc.mul_const(ex_to<numeric>(overall_coeff));
}


/** Compute the p-th power of a series.
 *
 *  @param p  power to compute
 *  @param deg  truncation order of series calculation */
ex pseries::power_const(const numeric &p, int deg) const
{
	// method:
	// (due to Leonhard Euler)
	// let A(x) be this series and for the time being let it start with a
	// constant (later we'll generalize):
	//     A(x) = a_0 + a_1*x + a_2*x^2 + ...
	// We want to compute
	//     C(x) = A(x)^p
	//     C(x) = c_0 + c_1*x + c_2*x^2 + ...
	// Taking the derivative on both sides and multiplying with A(x) one
	// immediately arrives at
	//     C'(x)*A(x) = p*C(x)*A'(x)
	// Multiplying this out and comparing coefficients we get the recurrence
	// formula
	//     c_i = (i*p*a_i*c_0 + ((i-1)*p-1)*a_{i-1}*c_1 + ...
	//                    ... + (p-(i-1))*a_1*c_{i-1})/(a_0*i)
	// which can easily be solved given the starting value c_0 = (a_0)^p.
	// For the more general case where the leading coefficient of A(x) is not
	// a constant, just consider A2(x) = A(x)*x^m, with some integer m and
	// repeat the above derivation.  The leading power of C2(x) = A2(x)^2 is
	// then of course x^(p*m) but the recurrence formula still holds.
	
	if (seq.empty()) {
		// as a special case, handle the empty (zero) series honoring the
		// usual power laws such as implemented in power::eval()
		if (p.real().is_zero())
			throw std::domain_error("pseries::power_const(): pow(0,I) is undefined");
		else if (p.real().is_negative())
			throw pole_error("pseries::power_const(): division by zero",1);
		else
			return *this;
	}
	
	const int ldeg = ldegree(var);
	if (!(p*ldeg).is_integer())
		throw std::runtime_error("pseries::power_const(): trying to assemble a Puiseux series");

	// adjust number of coefficients
	int numcoeff = deg - (p*ldeg).to_int();
	if (numcoeff <= 0) {
		epvector epv { expair(Order(_ex1), deg) };
		return dynallocate<pseries>(relational(var,point), std::move(epv));
	}
	
	// O(x^n)^(-m) is undefined
	if (seq.size() == 1 && is_order_function(seq[0].rest) && p.real().is_negative())
		throw pole_error("pseries::power_const(): division by zero",1);
	
	// Compute coefficients of the powered series
	exvector co;
	co.reserve(numcoeff);
	co.push_back(pow(coeff(var, ldeg), p));
	for (int i=1; i<numcoeff; ++i) {
		ex sum = _ex0;
		for (int j=1; j<=i; ++j) {
			ex c = coeff(var, j + ldeg);
			if (is_order_function(c)) {
				co.push_back(Order(_ex1));
				break;
			} else
				sum += (p * j - (i - j)) * co[i - j] * c;
		}
		co.push_back(sum / coeff(var, ldeg) / i);
	}
	
	// Construct new series (of non-zero coefficients)
	epvector new_seq;
	bool higher_order = false;
	for (int i=0; i<numcoeff; ++i) {
		if (!co[i].is_zero())
			new_seq.push_back(expair(co[i], p * ldeg + i));
		if (is_order_function(co[i])) {
			higher_order = true;
			break;
		}
	}
	if (!higher_order)
		new_seq.push_back(expair(Order(_ex1), p * ldeg + numcoeff));

	return pseries(relational(var,point), std::move(new_seq));
}


/** Return a new pseries object with the powers shifted by deg. */
pseries pseries::shift_exponents(int deg) const
{
	epvector newseq = seq;
	for (auto & it : newseq)
		it.coeff += deg;
	return pseries(relational(var, point), std::move(newseq));
}


/** Implementation of ex::series() for powers. This performs Laurent expansion
 *  of reciprocals of series at singularities.
 *  @see ex::series */
ex power::series(const relational & r, int order, unsigned options) const
{
	// If basis is already a series, just power it
	if (is_exactly_a<pseries>(basis))
		return ex_to<pseries>(basis).power_const(ex_to<numeric>(exponent), order);

	// Basis is not a series, may there be a singularity?
	bool must_expand_basis = false;
	try {
		basis.subs(r, subs_options::no_pattern);
	} catch (pole_error) {
		must_expand_basis = true;
	}

	bool exponent_is_regular = true;
	try {
		exponent.subs(r, subs_options::no_pattern);
	} catch (pole_error) {
		exponent_is_regular = false;
	}

	if (!exponent_is_regular) {
		ex l = exponent*log(basis);
		// this == exp(l);
		ex le = l.series(r, order, options);
		// Note: expanding exp(l) won't help, since that will attempt
		// Taylor expansion, and fail (because exponent is "singular")
		// Still l itself might be expanded in Taylor series.
		// Examples:
		// sin(x)/x*log(cos(x))
		// 1/x*log(1 + x)
		return exp(le).series(r, order, options);
		// Note: if l happens to have a Laurent expansion (with
		// negative powers of (var - point)), expanding exp(le)
		// will barf (which is The Right Thing).
	}

	// Is the expression of type something^(-int)?
	if (!must_expand_basis && !exponent.info(info_flags::negint)
	 && (!is_a<add>(basis) || !is_a<numeric>(exponent)))
		return basic::series(r, order, options);

	// Is the expression of type 0^something?
	if (!must_expand_basis && !basis.subs(r, subs_options::no_pattern).is_zero()
	 && (!is_a<add>(basis) || !is_a<numeric>(exponent)))
		return basic::series(r, order, options);

	// Singularity encountered, is the basis equal to (var - point)?
	if (basis.is_equal(r.lhs() - r.rhs())) {
		epvector new_seq;
		if (ex_to<numeric>(exponent).to_int() < order)
			new_seq.push_back(expair(_ex1, exponent));
		else
			new_seq.push_back(expair(Order(_ex1), exponent));
		return pseries(r, std::move(new_seq));
	}

	// No, expand basis into series

	numeric numexp;
	if (is_a<numeric>(exponent)) {
		numexp = ex_to<numeric>(exponent);
	} else {
		numexp = 0;
	}
	const ex& sym = r.lhs();
	// find existing minimal degree
	ex eb = basis.expand();
	int real_ldegree = 0;
	if (eb.info(info_flags::rational_function))
		real_ldegree = eb.ldegree(sym-r.rhs());
	if (real_ldegree == 0) {
		int orderloop = 0;
		do {
			orderloop++;
			real_ldegree = basis.series(r, orderloop, options).ldegree(sym);
		} while (real_ldegree == orderloop);
	}

	if (!(real_ldegree*numexp).is_integer())
		throw std::runtime_error("pseries::power_const(): trying to assemble a Puiseux series");
	ex e = basis.series(r, (order + real_ldegree*(1-numexp)).to_int(), options);
	
	ex result;
	try {
		result = ex_to<pseries>(e).power_const(numexp, order);
	} catch (pole_error) {
		epvector ser { expair(Order(_ex1), order) };
		result = pseries(r, std::move(ser));
	}

	return result;
}


/** Re-expansion of a pseries object. */
ex pseries::series(const relational & r, int order, unsigned options) const
{
	const ex p = r.rhs();
	GINAC_ASSERT(is_a<symbol>(r.lhs()));
	const symbol &s = ex_to<symbol>(r.lhs());
	
	if (var.is_equal(s) && point.is_equal(p)) {
		if (order > degree(s))
			return *this;
		else {
			epvector new_seq;
			for (auto & it : seq) {
				int o = ex_to<numeric>(it.coeff).to_int();
				if (o >= order) {
					new_seq.push_back(expair(Order(_ex1), o));
					break;
				}
				new_seq.push_back(it);
			}
			return pseries(r, std::move(new_seq));
		}
	} else
		return convert_to_poly().series(r, order, options);
}

ex integral::series(const relational & r, int order, unsigned options) const
{
	if (x.subs(r) != x)
		throw std::logic_error("Cannot series expand wrt dummy variable");
	
	// Expanding integrand with r substituted taken in boundaries.
	ex fseries = f.series(r, order, options);
	epvector fexpansion;
	fexpansion.reserve(fseries.nops());
	for (size_t i=0; i<fseries.nops(); ++i) {
		ex currcoeff = ex_to<pseries>(fseries).coeffop(i);
		currcoeff = (currcoeff == Order(_ex1))
			? currcoeff
			: integral(x, a.subs(r), b.subs(r), currcoeff);
		if (currcoeff != 0)
			fexpansion.push_back(
				expair(currcoeff, ex_to<pseries>(fseries).exponop(i)));
	}

	// Expanding lower boundary
	ex result = dynallocate<pseries>(r, std::move(fexpansion));
	ex aseries = (a-a.subs(r)).series(r, order, options);
	fseries = f.series(x == (a.subs(r)), order, options);
	for (size_t i=0; i<fseries.nops(); ++i) {
		ex currcoeff = ex_to<pseries>(fseries).coeffop(i);
		if (is_order_function(currcoeff))
			break;
		ex currexpon = ex_to<pseries>(fseries).exponop(i);
		int orderforf = order-ex_to<numeric>(currexpon).to_int()-1;
		currcoeff = currcoeff.series(r, orderforf);
		ex term = ex_to<pseries>(aseries).power_const(ex_to<numeric>(currexpon+1),order);
		term = ex_to<pseries>(term).mul_const(ex_to<numeric>(-1/(currexpon+1)));
		term = ex_to<pseries>(term).mul_series(ex_to<pseries>(currcoeff));
		result = ex_to<pseries>(result).add_series(ex_to<pseries>(term));
	}

	// Expanding upper boundary
	ex bseries = (b-b.subs(r)).series(r, order, options);
	fseries = f.series(x == (b.subs(r)), order, options);
	for (size_t i=0; i<fseries.nops(); ++i) {
		ex currcoeff = ex_to<pseries>(fseries).coeffop(i);
		if (is_order_function(currcoeff))
			break;
		ex currexpon = ex_to<pseries>(fseries).exponop(i);
		int orderforf = order-ex_to<numeric>(currexpon).to_int()-1;
		currcoeff = currcoeff.series(r, orderforf);
		ex term = ex_to<pseries>(bseries).power_const(ex_to<numeric>(currexpon+1),order);
		term = ex_to<pseries>(term).mul_const(ex_to<numeric>(1/(currexpon+1)));
		term = ex_to<pseries>(term).mul_series(ex_to<pseries>(currcoeff));
		result = ex_to<pseries>(result).add_series(ex_to<pseries>(term));
	}

	return result;
}


/** Compute the truncated series expansion of an expression.
 *  This function returns an expression containing an object of class pseries 
 *  to represent the series. If the series does not terminate within the given
 *  truncation order, the last term of the series will be an order term.
 *
 *  @param r  expansion relation, lhs holds variable and rhs holds point
 *  @param order  truncation order of series calculations
 *  @param options  of class series_options
 *  @return an expression holding a pseries object */
ex ex::series(const ex & r, int order, unsigned options) const
{
	ex e;
	relational rel_;
	
	if (is_a<relational>(r))
		rel_ = ex_to<relational>(r);
	else if (is_a<symbol>(r))
		rel_ = relational(r,_ex0);
	else
		throw (std::logic_error("ex::series(): expansion point has unknown type"));
	
	e = bp->series(rel_, order, options);
	return e;
}

GINAC_BIND_UNARCHIVER(pseries);

} // namespace GiNaC
