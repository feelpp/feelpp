/** @file ncmul.cpp
 *
 *  Implementation of GiNaC's non-commutative products of expressions. */

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

#include "ncmul.h"
#include "ex.h"
#include "add.h"
#include "mul.h"
#include "clifford.h"
#include "matrix.h"
#include "archive.h"
#include "indexed.h"
#include "utils.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(ncmul, exprseq,
  print_func<print_context>(&ncmul::do_print).
  print_func<print_tree>(&ncmul::do_print_tree).
  print_func<print_csrc>(&ncmul::do_print_csrc).
  print_func<print_python_repr>(&ncmul::do_print_csrc))


//////////
// default constructor
//////////

ncmul::ncmul()
{
}

//////////
// other constructors
//////////

// public

ncmul::ncmul(const ex & lh, const ex & rh) : inherited{lh,rh}
{
}

ncmul::ncmul(const ex & f1, const ex & f2, const ex & f3) : inherited{f1,f2,f3}
{
}

ncmul::ncmul(const ex & f1, const ex & f2, const ex & f3,
             const ex & f4) : inherited{f1,f2,f3,f4}
{
}

ncmul::ncmul(const ex & f1, const ex & f2, const ex & f3,
             const ex & f4, const ex & f5) : inherited{f1,f2,f3,f4,f5}
{
}

ncmul::ncmul(const ex & f1, const ex & f2, const ex & f3,
             const ex & f4, const ex & f5, const ex & f6) : inherited{f1,f2,f3,f4,f5,f6}
{
}

ncmul::ncmul(const exvector & v) : inherited(v)
{
}

ncmul::ncmul(exvector && v) : inherited(std::move(v))
{
}

//////////
// archiving
//////////


//////////
// functions overriding virtual functions from base classes
//////////

// public

void ncmul::do_print(const print_context & c, unsigned level) const
{
	printseq(c, '(', '*', ')', precedence(), level);
}

void ncmul::do_print_csrc(const print_context & c, unsigned level) const
{
	c.s << class_name();
	printseq(c, '(', ',', ')', precedence(), precedence());
}

bool ncmul::info(unsigned inf) const
{
	return inherited::info(inf);
}

typedef std::vector<std::size_t> uintvector;

ex ncmul::expand(unsigned options) const
{
	// First, expand the children
	exvector v = expandchildren(options);
	const exvector &expanded_seq = v.empty() ? this->seq : v;
	
	// Now, look for all the factors that are sums and remember their
	// position and number of terms.
	uintvector positions_of_adds(expanded_seq.size());
	uintvector number_of_add_operands(expanded_seq.size());

	size_t number_of_adds = 0;
	size_t number_of_expanded_terms = 1;

	size_t current_position = 0;
	for (auto & it : expanded_seq) {
		if (is_exactly_a<add>(it)) {
			positions_of_adds[number_of_adds] = current_position;
			size_t num_ops = it.nops();
			number_of_add_operands[number_of_adds] = num_ops;
			number_of_expanded_terms *= num_ops;
			number_of_adds++;
		}
		++current_position;
	}

	// If there are no sums, we are done
	if (number_of_adds == 0) {
		if (!v.empty())
			return dynallocate<ncmul>(std::move(v)).setflag(options == 0 ? status_flags::expanded : 0);
		else
			return *this;
	}

	// Now, form all possible products of the terms of the sums with the
	// remaining factors, and add them together
	exvector distrseq;
	distrseq.reserve(number_of_expanded_terms);

	uintvector k(number_of_adds);

	/* Rename indices in the static members of the product */
	exvector expanded_seq_mod;
	size_t j = 0;
	exvector va;

	for (size_t i=0; i<expanded_seq.size(); i++) {
		if (i == positions_of_adds[j]) {
			expanded_seq_mod.push_back(_ex1);
			j++;
		} else {
			expanded_seq_mod.push_back(rename_dummy_indices_uniquely(va, expanded_seq[i], true));
		}
	}

	while (true) {
		exvector term = expanded_seq_mod;
		for (size_t i=0; i<number_of_adds; i++) {
			term[positions_of_adds[i]] = rename_dummy_indices_uniquely(va, expanded_seq[positions_of_adds[i]].op(k[i]), true);
		}

		distrseq.push_back(dynallocate<ncmul>(std::move(term)).setflag(options == 0 ? status_flags::expanded : 0));

		// increment k[]
		int l = number_of_adds-1;
		while ((l>=0) && ((++k[l]) >= number_of_add_operands[l])) {
			k[l] = 0;
			l--;
		}
		if (l<0)
			break;
	}

	return dynallocate<add>(distrseq).setflag(options == 0 ? status_flags::expanded : 0);
}

int ncmul::degree(const ex & s) const
{
	if (is_equal(ex_to<basic>(s)))
		return 1;

	// Sum up degrees of factors
	int deg_sum = 0;
	for (auto & i : seq)
		deg_sum += i.degree(s);
	return deg_sum;
}

int ncmul::ldegree(const ex & s) const
{
	if (is_equal(ex_to<basic>(s)))
		return 1;

	// Sum up degrees of factors
	int deg_sum = 0;
	for (auto & i : seq)
		deg_sum += i.degree(s);
	return deg_sum;
}

ex ncmul::coeff(const ex & s, int n) const
{
	if (is_equal(ex_to<basic>(s)))
		return n==1 ? _ex1 : _ex0;

	exvector coeffseq;
	coeffseq.reserve(seq.size());

	if (n == 0) {
		// product of individual coeffs
		// if a non-zero power of s is found, the resulting product will be 0
		for (auto & it : seq)
			coeffseq.push_back(it.coeff(s,n));
		return dynallocate<ncmul>(std::move(coeffseq));
	}
		 
	bool coeff_found = false;
	for (auto & i : seq) {
		ex c = i.coeff(s,n);
		if (c.is_zero()) {
			coeffseq.push_back(i);
		} else {
			coeffseq.push_back(c);
			coeff_found = true;
		}
	}

	if (coeff_found)
		return dynallocate<ncmul>(std::move(coeffseq));
	
	return _ex0;
}

size_t ncmul::count_factors(const ex & e) const
{
	if ((is_exactly_a<mul>(e)&&(e.return_type()!=return_types::commutative))||
		(is_exactly_a<ncmul>(e))) {
		size_t factors=0;
		for (size_t i=0; i<e.nops(); i++)
			factors += count_factors(e.op(i));
		
		return factors;
	}
	return 1;
}
		
void ncmul::append_factors(exvector & v, const ex & e) const
{
	if ((is_exactly_a<mul>(e)&&(e.return_type()!=return_types::commutative))||
		(is_exactly_a<ncmul>(e))) {
		for (size_t i=0; i<e.nops(); i++)
			append_factors(v, e.op(i));
	} else 
		v.push_back(e);
}

typedef std::vector<unsigned> unsignedvector;
typedef std::vector<exvector> exvectorvector;

/** Perform automatic term rewriting rules in this class.  In the following
 *  x, x1, x2,... stand for a symbolic variables of type ex and c, c1, c2...
 *  stand for such expressions that contain a plain number.
 *  - ncmul(...,*(x1,x2),...,ncmul(x3,x4),...) -> ncmul(...,x1,x2,...,x3,x4,...)  (associativity)
 *  - ncmul(x) -> x
 *  - ncmul() -> 1
 *  - ncmul(...,c1,...,c2,...) -> *(c1,c2,ncmul(...))  (pull out commutative elements)
 *  - ncmul(x1,y1,x2,y2) -> *(ncmul(x1,x2),ncmul(y1,y2))  (collect elements of same type)
 *  - ncmul(x1,x2,x3,...) -> x::eval_ncmul(x1,x2,x3,...)
 */
ex ncmul::eval() const
{
	// The following additional rule would be nice, but produces a recursion,
	// which must be trapped by introducing a flag that the sub-ncmuls()
	// are already evaluated (maybe later...)
	//                  ncmul(x1,x2,...,X,y1,y2,...) ->
	//                      ncmul(ncmul(x1,x2,...),X,ncmul(y1,y2,...)
	//                      (X noncommutative_composite)

	if (flags & status_flags::evaluated) {
		return *this;
	}

	// ncmul(...,*(x1,x2),...,ncmul(x3,x4),...) ->
	//     ncmul(...,x1,x2,...,x3,x4,...)  (associativity)
	size_t factors = 0;
	for (auto & it : seq)
		factors += count_factors(it);
	
	exvector assocseq;
	assocseq.reserve(factors);
	make_flat_inserter mf(seq, true);
	for (auto & it : seq) {
		ex factor = mf.handle_factor(it, 1);
		append_factors(assocseq, factor);
	}
	
	// ncmul(x) -> x
	if (assocseq.size()==1) return *(seq.begin());

	// ncmul() -> 1
	if (assocseq.empty()) return _ex1;

	// determine return types
	unsignedvector rettypes(assocseq.size());
	size_t i = 0;
	size_t count_commutative=0;
	size_t count_noncommutative=0;
	size_t count_noncommutative_composite=0;
	for (auto & it : assocseq) {
		rettypes[i] = it.return_type();
		switch (rettypes[i]) {
		case return_types::commutative:
			count_commutative++;
			break;
		case return_types::noncommutative:
			count_noncommutative++;
			break;
		case return_types::noncommutative_composite:
			count_noncommutative_composite++;
			break;
		default:
			throw(std::logic_error("ncmul::eval(): invalid return type"));
		}
		++i;
	}
	GINAC_ASSERT(count_commutative+count_noncommutative+count_noncommutative_composite==assocseq.size());

	// ncmul(...,c1,...,c2,...) ->
	//     *(c1,c2,ncmul(...)) (pull out commutative elements)
	if (count_commutative!=0) {
		exvector commutativeseq;
		commutativeseq.reserve(count_commutative+1);
		exvector noncommutativeseq;
		noncommutativeseq.reserve(assocseq.size()-count_commutative);
		size_t num = assocseq.size();
		for (size_t i=0; i<num; ++i) {
			if (rettypes[i]==return_types::commutative)
				commutativeseq.push_back(assocseq[i]);
			else
				noncommutativeseq.push_back(assocseq[i]);
		}
		commutativeseq.push_back(dynallocate<ncmul>(std::move(noncommutativeseq)));
		return dynallocate<mul>(std::move(commutativeseq));
	}
		
	// ncmul(x1,y1,x2,y2) -> *(ncmul(x1,x2),ncmul(y1,y2))
	//     (collect elements of same type)

	if (count_noncommutative_composite==0) {
		// there are neither commutative nor noncommutative_composite
		// elements in assocseq
		GINAC_ASSERT(count_commutative==0);

		size_t assoc_num = assocseq.size();
		exvectorvector evv;
		std::vector<return_type_t> rttinfos;
		evv.reserve(assoc_num);
		rttinfos.reserve(assoc_num);

		for (auto & it : assocseq) {
			return_type_t ti = it.return_type_tinfo();
			size_t rtt_num = rttinfos.size();
			// search type in vector of known types
			for (i=0; i<rtt_num; ++i) {
				if(ti == rttinfos[i]) {
					evv[i].push_back(it);
					break;
				}
			}
			if (i >= rtt_num) {
				// new type
				rttinfos.push_back(ti);
				evv.push_back(exvector());
				(evv.end()-1)->reserve(assoc_num);
				(evv.end()-1)->push_back(it);
			}
		}

		size_t evv_num = evv.size();
#ifdef DO_GINAC_ASSERT
		GINAC_ASSERT(evv_num == rttinfos.size());
		GINAC_ASSERT(evv_num > 0);
		size_t s=0;
		for (i=0; i<evv_num; ++i)
			s += evv[i].size();
		GINAC_ASSERT(s == assoc_num);
#endif // def DO_GINAC_ASSERT
		
		// if all elements are of same type, simplify the string
		if (evv_num == 1) {
			return evv[0][0].eval_ncmul(evv[0]);
		}
		
		exvector splitseq;
		splitseq.reserve(evv_num);
		for (i=0; i<evv_num; ++i)
			splitseq.push_back(dynallocate<ncmul>(evv[i]));
		
		return dynallocate<mul>(splitseq);
	}
	
	return dynallocate<ncmul>(assocseq).setflag(status_flags::evaluated);
}

ex ncmul::evalm() const
{
	// Evaluate children first
	exvector s;
	s.reserve(seq.size());
	for (auto & it : seq)
		s.push_back(it.evalm());

	// If there are only matrices, simply multiply them
	auto it = s.begin(), itend = s.end();
	if (is_a<matrix>(*it)) {
		matrix prod(ex_to<matrix>(*it));
		it++;
		while (it != itend) {
			if (!is_a<matrix>(*it))
				goto no_matrix;
			prod = prod.mul(ex_to<matrix>(*it));
			it++;
		}
		return prod;
	}

no_matrix:
	return dynallocate<ncmul>(std::move(s));
}

ex ncmul::thiscontainer(const exvector & v) const
{
	return dynallocate<ncmul>(v);
}

ex ncmul::thiscontainer(exvector && v) const
{
	return dynallocate<ncmul>(std::move(v));
}

ex ncmul::conjugate() const
{
	if (return_type() != return_types::noncommutative) {
		return exprseq::conjugate();
	}

	if (!is_clifford_tinfo(return_type_tinfo())) {
		return exprseq::conjugate();
	}

	exvector ev;
	ev.reserve(nops());
	for (auto i=end(); i!=begin();) {
		--i;
		ev.push_back(i->conjugate());
	}
	return dynallocate<ncmul>(std::move(ev));
}

ex ncmul::real_part() const
{
	return basic::real_part();
}

ex ncmul::imag_part() const
{
	return basic::imag_part();
}

// protected

/** Implementation of ex::diff() for a non-commutative product. It applies
 *  the product rule.
 *  @see ex::diff */
ex ncmul::derivative(const symbol & s) const
{
	size_t num = seq.size();
	exvector addseq;
	addseq.reserve(num);
	
	// D(a*b*c) = D(a)*b*c + a*D(b)*c + a*b*D(c)
	exvector ncmulseq = seq;
	for (size_t i=0; i<num; ++i) {
		ex e = seq[i].diff(s);
		e.swap(ncmulseq[i]);
		addseq.push_back(dynallocate<ncmul>(ncmulseq));
		e.swap(ncmulseq[i]);
	}
	return dynallocate<add>(addseq);
}

int ncmul::compare_same_type(const basic & other) const
{
	return inherited::compare_same_type(other);
}

unsigned ncmul::return_type() const
{
	if (seq.empty())
		return return_types::commutative;

	bool all_commutative = true;
	exvector::const_iterator noncommutative_element; // point to first found nc element

	auto i = seq.begin(), end = seq.end();
	while (i != end) {
		unsigned rt = i->return_type();
		if (rt == return_types::noncommutative_composite)
			return rt; // one ncc -> mul also ncc
		if ((rt == return_types::noncommutative) && (all_commutative)) {
			// first nc element found, remember position
			noncommutative_element = i;
			all_commutative = false;
		}
		if ((rt == return_types::noncommutative) && (!all_commutative)) {
			// another nc element found, compare type_infos
			if(noncommutative_element->return_type_tinfo() != i->return_type_tinfo())
					return return_types::noncommutative_composite;
		}
		++i;
	}
	// all factors checked
	GINAC_ASSERT(!all_commutative); // not all factors should commutate, because this is a ncmul();
	return all_commutative ? return_types::commutative : return_types::noncommutative;
}

return_type_t ncmul::return_type_tinfo() const
{
	if (seq.empty())
		return make_return_type_t<ncmul>();

	// return type_info of first noncommutative element
	for (auto & i : seq)
		if (i.return_type() == return_types::noncommutative)
			return i.return_type_tinfo();

	// no noncommutative element found, should not happen
	return make_return_type_t<ncmul>();
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

exvector ncmul::expandchildren(unsigned options) const
{
	auto cit = this->seq.begin(), end = this->seq.end();
	while (cit != end) {
		const ex & expanded_ex = cit->expand(options);
		if (!are_ex_trivially_equal(*cit, expanded_ex)) {

			// copy first part of seq which hasn't changed
			exvector s(this->seq.begin(), cit);
			s.reserve(this->seq.size());

			// insert changed element
			s.push_back(expanded_ex);
			++cit;

			// copy rest
			while (cit != end) {
				s.push_back(cit->expand(options));
				++cit;
			}

			return s;
		}

		++cit;
	}

	return exvector(); // nothing has changed
}

const exvector & ncmul::get_factors() const
{
	return seq;
}

//////////
// friend functions
//////////

ex reeval_ncmul(const exvector & v)
{
	return dynallocate<ncmul>(v);
}

ex hold_ncmul(const exvector & v)
{
	if (v.empty())
		return _ex1;
	else if (v.size() == 1)
		return v[0];
	else
		return dynallocate<ncmul>(v).setflag(status_flags::evaluated);
}

GINAC_BIND_UNARCHIVER(ncmul);

} // namespace GiNaC
