/** @file expairseq.cpp
 *
 *  Implementation of sequences of expression pairs. */

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

#include "expairseq.h"
#include "lst.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "relational.h"
#include "wildcard.h"
#include "archive.h"
#include "operators.h"
#include "utils.h"
#include "hash_seed.h"
#include "indexed.h"
#include "compiler.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>

namespace GiNaC {

	
GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(expairseq, basic,
  print_func<print_context>(&expairseq::do_print).
  print_func<print_tree>(&expairseq::do_print_tree))


//////////
// helper classes
//////////

class epp_is_less
{
public:
	bool operator()(const epp &lh, const epp &rh) const
	{
		return (*lh).is_less(*rh);
	}
};

//////////
// default constructor
//////////

// public

expairseq::expairseq() 
{}

// protected

//////////
// other constructors
//////////

expairseq::expairseq(const ex &lh, const ex &rh)
{
	construct_from_2_ex(lh,rh);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(const exvector &v)
{
	construct_from_exvector(v);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(const epvector &v, const ex &oc, bool do_index_renaming)
  :  overall_coeff(oc)
{
	GINAC_ASSERT(is_a<numeric>(oc));
	construct_from_epvector(v, do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(epvector && vp, const ex &oc, bool do_index_renaming)
  :  overall_coeff(oc)
{
	GINAC_ASSERT(is_a<numeric>(oc));
	construct_from_epvector(std::move(vp), do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

//////////
// archiving
//////////

void expairseq::read_archive(const archive_node &n, lst &sym_lst) 
{
	inherited::read_archive(n, sym_lst);
	auto first = n.find_first("rest");
	auto last = n.find_last("coeff");
	++last;
	seq.reserve((last-first)/2);

	for (auto loc = first; loc < last;) {
		ex rest;
		ex coeff;
		n.find_ex_by_loc(loc++, rest, sym_lst);
		n.find_ex_by_loc(loc++, coeff, sym_lst);
		seq.push_back(expair(rest, coeff));
	}

	n.find_ex("overall_coeff", overall_coeff, sym_lst);

	canonicalize();
	GINAC_ASSERT(is_canonical());
}

void expairseq::archive(archive_node &n) const
{
	inherited::archive(n);
	for (auto & i : seq) {
		n.add_ex("rest", i.rest);
		n.add_ex("coeff", i.coeff);
	}
	n.add_ex("overall_coeff", overall_coeff);
}


//////////
// functions overriding virtual functions from base classes
//////////

// public

void expairseq::do_print(const print_context & c, unsigned level) const
{
	c.s << "[[";
	printseq(c, ',', precedence(), level);
	c.s << "]]";
}

void expairseq::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", nops=" << nops()
	    << std::endl;
	size_t num = seq.size();
	for (size_t i=0; i<num; ++i) {
		seq[i].rest.print(c, level + c.delta_indent);
		seq[i].coeff.print(c, level + c.delta_indent);
		if (i != num - 1)
			c.s << std::string(level + c.delta_indent, ' ') << "-----" << std::endl;
	}
	if (!overall_coeff.is_equal(default_overall_coeff())) {
		c.s << std::string(level + c.delta_indent, ' ') << "-----" << std::endl
		    << std::string(level + c.delta_indent, ' ') << "overall_coeff" << std::endl;
		overall_coeff.print(c, level + c.delta_indent);
	}
	c.s << std::string(level + c.delta_indent,' ') << "=====" << std::endl;
}

bool expairseq::info(unsigned inf) const
{
	switch(inf) {
		case info_flags::expanded:
			return (flags & status_flags::expanded);
		case info_flags::has_indices: {
			if (flags & status_flags::has_indices)
				return true;
			else if (flags & status_flags::has_no_indices)
				return false;
			for (auto & i : seq) {
				if (i.rest.info(info_flags::has_indices)) {
					this->setflag(status_flags::has_indices);
					this->clearflag(status_flags::has_no_indices);
					return true;
				}
			}
			this->clearflag(status_flags::has_indices);
			this->setflag(status_flags::has_no_indices);
			return false;
		}
	}
	return inherited::info(inf);
}

size_t expairseq::nops() const
{
	if (overall_coeff.is_equal(default_overall_coeff()))
		return seq.size();
	else
		return seq.size()+1;
}

ex expairseq::op(size_t i) const
{
	if (i < seq.size())
		return recombine_pair_to_ex(seq[i]);
	GINAC_ASSERT(!overall_coeff.is_equal(default_overall_coeff()));
	return overall_coeff;
}

ex expairseq::map(map_function &f) const
{
	epvector v;
	v.reserve(seq.size()+1);

	for (auto & it : seq)
		v.push_back(split_ex_to_pair(f(recombine_pair_to_ex(it))));

	if (overall_coeff.is_equal(default_overall_coeff()))
		return thisexpairseq(std::move(v), default_overall_coeff(), true);
	else {
		ex newcoeff = f(overall_coeff);
		if(is_a<numeric>(newcoeff))
			return thisexpairseq(std::move(v), newcoeff, true);
		else {
			v.push_back(split_ex_to_pair(newcoeff));
			return thisexpairseq(std::move(v), default_overall_coeff(), true);
		}
	}
}

/** Perform coefficient-wise automatic term rewriting rules in this class. */
ex expairseq::eval() const
{
	if (flags &status_flags::evaluated)
		return *this;

	const epvector evaled = evalchildren();
	if (!evaled.empty())
		return dynallocate<expairseq>(std::move(evaled), overall_coeff).setflag(status_flags::evaluated);
	else
		return this->hold();
}

epvector* conjugateepvector(const epvector&epv)
{
	epvector *newepv = nullptr;
	for (auto i=epv.begin(); i!=epv.end(); ++i) {
		if (newepv) {
			newepv->push_back(i->conjugate());
			continue;
		}
		expair x = i->conjugate();
		if (x.is_equal(*i)) {
			continue;
		}
		newepv = new epvector;
		newepv->reserve(epv.size());
		for (auto j=epv.begin(); j!=i; ++j) {
			newepv->push_back(*j);
		}
		newepv->push_back(x);
	}
	return newepv;
}

ex expairseq::conjugate() const
{
	std::unique_ptr<epvector> newepv(conjugateepvector(seq));
	ex x = overall_coeff.conjugate();
	if (newepv) {
		return thisexpairseq(std::move(*newepv), x);
	}
	if (are_ex_trivially_equal(x, overall_coeff)) {
		return *this;
	}
	return thisexpairseq(seq, x);
}

bool expairseq::match(const ex & pattern, exmap & repl_lst) const
{
	// This differs from basic::match() because we want "a+b+c+d" to
	// match "d+*+b" with "*" being "a+c", and we want to honor commutativity

	if (typeid(*this) == typeid(ex_to<basic>(pattern))) {

		// Check whether global wildcard (one that matches the "rest of the
		// expression", like "*" above) is present
		bool has_global_wildcard = false;
		ex global_wildcard;
		for (size_t i=0; i<pattern.nops(); i++) {
			if (is_exactly_a<wildcard>(pattern.op(i))) {
				has_global_wildcard = true;
				global_wildcard = pattern.op(i);
				break;
			}
		}

		// Even if the expression does not match the pattern, some of
		// its subexpressions could match it. For example, x^5*y^(-1)
		// does not match the pattern $0^5, but its subexpression x^5
		// does. So, save repl_lst in order to not add bogus entries.
		exmap tmp_repl = repl_lst;

		// Unfortunately, this is an O(N^2) operation because we can't
		// sort the pattern in a useful way...

		// Chop into terms
		exvector ops;
		ops.reserve(nops());
		for (size_t i=0; i<nops(); i++)
			ops.push_back(op(i));

		// Now, for every term of the pattern, look for a matching term in
		// the expression and remove the match
		for (size_t i=0; i<pattern.nops(); i++) {
			ex p = pattern.op(i);
			if (has_global_wildcard && p.is_equal(global_wildcard))
				continue;
			auto it = ops.begin(), itend = ops.end();
			while (it != itend) {
				if (it->match(p, tmp_repl)) {
					ops.erase(it);
					goto found;
				}
				++it;
			}
			return false; // no match found
found:		;
		}

		if (has_global_wildcard) {

			// Assign all the remaining terms to the global wildcard (unless
			// it has already been matched before, in which case the matches
			// must be equal)
			size_t num = ops.size();
			epvector vp;
			vp.reserve(num);
			for (size_t i=0; i<num; i++)
				vp.push_back(split_ex_to_pair(ops[i]));
			ex rest = thisexpairseq(std::move(vp), default_overall_coeff());
			for (auto & it : tmp_repl) {
				if (it.first.is_equal(global_wildcard)) {
					if (rest.is_equal(it.second)) {
						repl_lst = tmp_repl;
						return true;
					}
					return false;
				}
			}
			repl_lst = tmp_repl;
			repl_lst[global_wildcard] = rest;
			return true;

		} else {

			// No global wildcard, then the match fails if there are any
			// unmatched terms left
			if (ops.empty()) {
				repl_lst = tmp_repl;
				return true;
			}
			return false;
		}
	}
	return inherited::match(pattern, repl_lst);
}

ex expairseq::subs(const exmap & m, unsigned options) const
{
	epvector subsed = subschildren(m, options);
	if (!subsed.empty())
		return ex_to<basic>(thisexpairseq(std::move(subsed), overall_coeff, (options & subs_options::no_index_renaming) == 0));
	else if ((options & subs_options::algebraic) && is_exactly_a<mul>(*this))
		return static_cast<const mul *>(this)->algebraic_subs_mul(m, options);
	else
		return subs_one_level(m, options);
}

// protected

int expairseq::compare_same_type(const basic &other) const
{
	GINAC_ASSERT(is_a<expairseq>(other));
	const expairseq &o = static_cast<const expairseq &>(other);
	
	int cmpval;
	
	// compare number of elements
	if (seq.size() != o.seq.size())
		return (seq.size()<o.seq.size()) ? -1 : 1;
	
	// compare overall_coeff
	cmpval = overall_coeff.compare(o.overall_coeff);
	if (cmpval!=0)
		return cmpval;
	
	auto cit1 = seq.begin(), last1 = seq.end();
	auto cit2 = o.seq.begin(), last2 = o.seq.end();
	for (; (cit1!=last1) && (cit2!=last2); ++cit1, ++cit2) {
		cmpval = (*cit1).compare(*cit2);
		if (cmpval!=0) return cmpval;
	}
		
	GINAC_ASSERT(cit1==last1);
	GINAC_ASSERT(cit2==last2);
		
	return 0;
}

bool expairseq::is_equal_same_type(const basic &other) const
{
	const expairseq &o = static_cast<const expairseq &>(other);
	
	// compare number of elements
	if (seq.size()!=o.seq.size())
		return false;
	
	// compare overall_coeff
	if (!overall_coeff.is_equal(o.overall_coeff))
		return false;
	
	auto cit2 = o.seq.begin();
	for (auto & cit1 : seq) {
		if (!cit1.is_equal(*cit2))
			return false;
		++cit2;
	}

	return true;
}

unsigned expairseq::return_type() const
{
	return return_types::noncommutative_composite;
}

unsigned expairseq::calchash() const
{
	unsigned v = make_hash_seed(typeid(*this));
	for (auto & i : seq) {
		v ^= i.rest.gethash();
		v = rotate_left(v);
		v ^= i.coeff.gethash();
	}

	v ^= overall_coeff.gethash();

	// store calculated hash value only if object is already evaluated
	if (flags &status_flags::evaluated) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}
	
	return v;
}

ex expairseq::expand(unsigned options) const
{
	epvector expanded = expandchildren(options);
	if (!expanded.empty()) {
		return thisexpairseq(std::move(expanded), overall_coeff);
	}
	return (options == 0) ? setflag(status_flags::expanded) : *this;
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// protected

/** Create an object of this type.
 *  This method works similar to a constructor.  It is useful because expairseq
 *  has (at least) two possible different semantics but we want to inherit
 *  methods thus avoiding code duplication.  Sometimes a method in expairseq
 *  has to create a new one of the same semantics, which cannot be done by a
 *  ctor because the name (add, mul,...) is unknown on the expairseq level.  In
 *  order for this trick to work a derived class must of course override this
 *  definition. */
ex expairseq::thisexpairseq(const epvector &v, const ex &oc, bool do_index_renaming) const
{
	return expairseq(v, oc, do_index_renaming);
}

ex expairseq::thisexpairseq(epvector && vp, const ex &oc, bool do_index_renaming) const
{
	return expairseq(std::move(vp), oc, do_index_renaming);
}

void expairseq::printpair(const print_context & c, const expair & p, unsigned upper_precedence) const
{
	c.s << "[[";
	p.rest.print(c, precedence());
	c.s << ",";
	p.coeff.print(c, precedence());
	c.s << "]]";
}

void expairseq::printseq(const print_context & c, char delim,
                         unsigned this_precedence,
                         unsigned upper_precedence) const
{
	if (this_precedence <= upper_precedence)
		c.s << "(";
	auto it = seq.begin(), it_last = seq.end() - 1;
	for (; it!=it_last; ++it) {
		printpair(c, *it, this_precedence);
		c.s << delim;
	}
	printpair(c, *it, this_precedence);
	if (!overall_coeff.is_equal(default_overall_coeff())) {
		c.s << delim;
		overall_coeff.print(c, this_precedence);
	}
	
	if (this_precedence <= upper_precedence)
		c.s << ")";
}


/** Form an expair from an ex, using the corresponding semantics.
 *  @see expairseq::recombine_pair_to_ex() */
expair expairseq::split_ex_to_pair(const ex &e) const
{
	return expair(e,_ex1);
}


expair expairseq::combine_ex_with_coeff_to_pair(const ex &e,
                                                const ex &c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	
	return expair(e,c);
}


expair expairseq::combine_pair_with_coeff_to_pair(const expair &p,
                                                  const ex &c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	
	return expair(p.rest,ex_to<numeric>(p.coeff).mul_dyn(ex_to<numeric>(c)));
}


/** Form an ex out of an expair, using the corresponding semantics.
 *  @see expairseq::split_ex_to_pair() */
ex expairseq::recombine_pair_to_ex(const expair &p) const
{
	return lst{p.rest, p.coeff};
}

bool expairseq::expair_needs_further_processing(epp it)
{
	return false;
}

ex expairseq::default_overall_coeff() const
{
	return _ex0;
}

void expairseq::combine_overall_coeff(const ex &c)
{
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	overall_coeff = ex_to<numeric>(overall_coeff).add_dyn(ex_to<numeric>(c));
}

void expairseq::combine_overall_coeff(const ex &c1, const ex &c2)
{
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c1));
	GINAC_ASSERT(is_exactly_a<numeric>(c2));
	overall_coeff = ex_to<numeric>(overall_coeff).
	                add_dyn(ex_to<numeric>(c1).mul(ex_to<numeric>(c2)));
}

bool expairseq::can_make_flat(const expair &p) const
{
	return true;
}


//////////
// non-virtual functions in this class
//////////

void expairseq::construct_from_2_ex(const ex &lh, const ex &rh)
{
	const std::type_info& typeid_this = typeid(*this);
	if (typeid(ex_to<basic>(lh)) == typeid_this) {
		if (typeid(ex_to<basic>(rh)) == typeid_this) {
			if (is_a<mul>(lh) && lh.info(info_flags::has_indices) && 
				rh.info(info_flags::has_indices)) {
				ex newrh=rename_dummy_indices_uniquely(lh, rh);
				construct_from_2_expairseq(ex_to<expairseq>(lh),
				                           ex_to<expairseq>(newrh));
			}
			else
				construct_from_2_expairseq(ex_to<expairseq>(lh),
				                           ex_to<expairseq>(rh));
			return;
		} else {
			construct_from_expairseq_ex(ex_to<expairseq>(lh), rh);
			return;
		}
	} else if (typeid(ex_to<basic>(rh)) == typeid_this) {
		construct_from_expairseq_ex(ex_to<expairseq>(rh),lh);
		return;
	}
	
	if (is_exactly_a<numeric>(lh)) {
		if (is_exactly_a<numeric>(rh)) {
			combine_overall_coeff(lh);
			combine_overall_coeff(rh);
		} else {
			combine_overall_coeff(lh);
			seq.push_back(split_ex_to_pair(rh));
		}
	} else {
		if (is_exactly_a<numeric>(rh)) {
			combine_overall_coeff(rh);
			seq.push_back(split_ex_to_pair(lh));
		} else {
			expair p1 = split_ex_to_pair(lh);
			expair p2 = split_ex_to_pair(rh);
			
			int cmpval = p1.rest.compare(p2.rest);
			if (cmpval==0) {
				p1.coeff = ex_to<numeric>(p1.coeff).add_dyn(ex_to<numeric>(p2.coeff));
				if (!ex_to<numeric>(p1.coeff).is_zero()) {
					// no further processing is necessary, since this
					// one element will usually be recombined in eval()
					seq.push_back(p1);
				}
			} else {
				seq.reserve(2);
				if (cmpval<0) {
					seq.push_back(p1);
					seq.push_back(p2);
				} else {
					seq.push_back(p2);
					seq.push_back(p1);
				}
			}
		}
	}
}

void expairseq::construct_from_2_expairseq(const expairseq &s1,
                                           const expairseq &s2)
{
	combine_overall_coeff(s1.overall_coeff);
	combine_overall_coeff(s2.overall_coeff);

	auto first1 = s1.seq.begin(), last1 = s1.seq.end();
	auto first2 = s2.seq.begin(), last2 = s2.seq.end();

	seq.reserve(s1.seq.size()+s2.seq.size());

	bool needs_further_processing=false;
	
	while (first1!=last1 && first2!=last2) {
		int cmpval = (*first1).rest.compare((*first2).rest);

		if (cmpval==0) {
			// combine terms
			const numeric &newcoeff = ex_to<numeric>(first1->coeff).
			                           add(ex_to<numeric>(first2->coeff));
			if (!newcoeff.is_zero()) {
				seq.push_back(expair(first1->rest,newcoeff));
				if (expair_needs_further_processing(seq.end()-1)) {
					needs_further_processing = true;
				}
			}
			++first1;
			++first2;
		} else if (cmpval<0) {
			seq.push_back(*first1);
			++first1;
		} else {
			seq.push_back(*first2);
			++first2;
		}
	}
	
	while (first1!=last1) {
		seq.push_back(*first1);
		++first1;
	}
	while (first2!=last2) {
		seq.push_back(*first2);
		++first2;
	}
	
	if (needs_further_processing) {
		// Clear seq and start over.
		epvector v = std::move(seq);
		construct_from_epvector(std::move(v));
	}
}

void expairseq::construct_from_expairseq_ex(const expairseq &s,
                                            const ex &e)
{
	combine_overall_coeff(s.overall_coeff);
	if (is_exactly_a<numeric>(e)) {
		combine_overall_coeff(e);
		seq = s.seq;
		return;
	}
	
	auto first = s.seq.begin(), last = s.seq.end();
	expair p = split_ex_to_pair(e);
	
	seq.reserve(s.seq.size()+1);
	bool p_pushed = false;
	
	bool needs_further_processing=false;
	
	// merge p into s.seq
	while (first!=last) {
		int cmpval = (*first).rest.compare(p.rest);
		if (cmpval==0) {
			// combine terms
			const numeric &newcoeff = ex_to<numeric>(first->coeff).
			                           add(ex_to<numeric>(p.coeff));
			if (!newcoeff.is_zero()) {
				seq.push_back(expair(first->rest,newcoeff));
				if (expair_needs_further_processing(seq.end()-1))
					needs_further_processing = true;
			}
			++first;
			p_pushed = true;
			break;
		} else if (cmpval<0) {
			seq.push_back(*first);
			++first;
		} else {
			seq.push_back(p);
			p_pushed = true;
			break;
		}
	}
	
	if (p_pushed) {
		// while loop exited because p was pushed, now push rest of s.seq
		while (first!=last) {
			seq.push_back(*first);
			++first;
		}
	} else {
		// while loop exited because s.seq was pushed, now push p
		seq.push_back(p);
	}

	if (needs_further_processing) {
		// Clear seq and start over.
		epvector v = std::move(seq);
		construct_from_epvector(std::move(v));
	}
}

void expairseq::construct_from_exvector(const exvector &v)
{
	// simplifications: +(a,+(b,c),d) -> +(a,b,c,d) (associativity)
	//                  +(d,b,c,a) -> +(a,b,c,d) (canonicalization)
	//                  +(...,x,*(x,c1),*(x,c2)) -> +(...,*(x,1+c1+c2)) (c1, c2 numeric)
	//                  (same for (+,*) -> (*,^)

	make_flat(v);
	canonicalize();
	combine_same_terms_sorted_seq();
}

void expairseq::construct_from_epvector(const epvector &v, bool do_index_renaming)
{
	// simplifications: +(a,+(b,c),d) -> +(a,b,c,d) (associativity)
	//                  +(d,b,c,a) -> +(a,b,c,d) (canonicalization)
	//                  +(...,x,*(x,c1),*(x,c2)) -> +(...,*(x,1+c1+c2)) (c1, c2 numeric)
	//                  same for (+,*) -> (*,^)

	make_flat(v, do_index_renaming);
	canonicalize();
	combine_same_terms_sorted_seq();
}

void expairseq::construct_from_epvector(epvector &&v, bool do_index_renaming)
{
	// simplifications: +(a,+(b,c),d) -> +(a,b,c,d) (associativity)
	//                  +(d,b,c,a) -> +(a,b,c,d) (canonicalization)
	//                  +(...,x,*(x,c1),*(x,c2)) -> +(...,*(x,1+c1+c2)) (c1, c2 numeric)
	//                  same for (+,*) -> (*,^)

	make_flat(std::move(v), do_index_renaming);
	canonicalize();
	combine_same_terms_sorted_seq();
}

/** Combine this expairseq with argument exvector.
 *  It cares for associativity as well as for special handling of numerics. */
void expairseq::make_flat(const exvector &v)
{
	// count number of operands which are of same expairseq derived type
	// and their cumulative number of operands
	int nexpairseqs = 0;
	int noperands = 0;
	bool do_idx_rename = false;

	const std::type_info& typeid_this = typeid(*this);
	for (auto & cit : v) {
		if (typeid(ex_to<basic>(cit)) == typeid_this) {
			++nexpairseqs;
			noperands += ex_to<expairseq>(cit).seq.size();
		}
		if (is_a<mul>(*this) && (!do_idx_rename) &&
		    cit.info(info_flags::has_indices))
			do_idx_rename = true;
	}
	
	// reserve seq and coeffseq which will hold all operands
	seq.reserve(v.size()+noperands-nexpairseqs);
	
	// copy elements and split off numerical part
	make_flat_inserter mf(v, do_idx_rename);
	for (auto & cit : v) {
		if (typeid(ex_to<basic>(cit)) == typeid_this) {
			ex newfactor = mf.handle_factor(cit, _ex1);
			const expairseq &subseqref = ex_to<expairseq>(newfactor);
			combine_overall_coeff(subseqref.overall_coeff);
			for (auto & cit_s : subseqref.seq) {
				seq.push_back(cit_s);
			}
		} else {
			if (is_exactly_a<numeric>(cit))
				combine_overall_coeff(cit);
			else {
				ex newfactor = mf.handle_factor(cit, _ex1);
				seq.push_back(split_ex_to_pair(newfactor));
			}
		}
	}
}

/** Combine this expairseq with argument epvector.
 *  It cares for associativity as well as for special handling of numerics. */
void expairseq::make_flat(const epvector &v, bool do_index_renaming)
{
	// count number of operands which are of same expairseq derived type
	// and their cumulative number of operands
	int nexpairseqs = 0;
	int noperands = 0;
	bool really_need_rename_inds = false;

	const std::type_info& typeid_this = typeid(*this);
	for (auto & cit : v) {
		if (typeid(ex_to<basic>(cit.rest)) == typeid_this) {
			++nexpairseqs;
			noperands += ex_to<expairseq>(cit.rest).seq.size();
		}
		if ((!really_need_rename_inds) && is_a<mul>(*this) &&
		    cit.rest.info(info_flags::has_indices))
			really_need_rename_inds = true;
	}
	do_index_renaming = do_index_renaming && really_need_rename_inds;
	
	// reserve seq and coeffseq which will hold all operands
	seq.reserve(v.size()+noperands-nexpairseqs);
	make_flat_inserter mf(v, do_index_renaming);
	
	// copy elements and split off numerical part
	for (auto & cit : v) {
		if (typeid(ex_to<basic>(cit.rest)) == typeid_this &&
		    this->can_make_flat(cit)) {
			ex newrest = mf.handle_factor(cit.rest, cit.coeff);
			const expairseq &subseqref = ex_to<expairseq>(newrest);
			combine_overall_coeff(ex_to<numeric>(subseqref.overall_coeff),
			                      ex_to<numeric>(cit.coeff));
			for (auto & cit_s : subseqref.seq) {
				seq.push_back(expair(cit_s.rest,
				                     ex_to<numeric>(cit_s.coeff).mul_dyn(ex_to<numeric>(cit.coeff))));
			}
		} else {
			if (cit.is_canonical_numeric())
				combine_overall_coeff(mf.handle_factor(cit.rest, _ex1));
			else {
				ex rest = cit.rest;
				ex newrest = mf.handle_factor(rest, cit.coeff);
				if (are_ex_trivially_equal(newrest, rest))
					seq.push_back(cit);
				else
					seq.push_back(expair(newrest, cit.coeff));
			}
		}
	}
}

/** Brings this expairseq into a sorted (canonical) form. */
void expairseq::canonicalize()
{
	std::sort(seq.begin(), seq.end(), expair_rest_is_less());
}


/** Compact a presorted expairseq by combining all matching expairs to one
 *  each.  On an add object, this is responsible for 2*x+3*x+y -> 5*x+y, for
 *  instance. */
void expairseq::combine_same_terms_sorted_seq()
{
	if (seq.size()<2)
		return;

	bool needs_further_processing = false;

	auto itin1 = seq.begin();
	auto itin2 = itin1 + 1;
	auto itout = itin1;
	auto last = seq.end();
	// must_copy will be set to true the first time some combination is 
	// possible from then on the sequence has changed and must be compacted
	bool must_copy = false;
	while (itin2!=last) {
		if (itin1->rest.compare(itin2->rest)==0) {
			itin1->coeff = ex_to<numeric>(itin1->coeff).
			               add_dyn(ex_to<numeric>(itin2->coeff));
			if (expair_needs_further_processing(itin1))
				needs_further_processing = true;
			must_copy = true;
		} else {
			if (!ex_to<numeric>(itin1->coeff).is_zero()) {
				if (must_copy)
					*itout = *itin1;
				++itout;
			}
			itin1 = itin2;
		}
		++itin2;
	}
	if (!ex_to<numeric>(itin1->coeff).is_zero()) {
		if (must_copy)
			*itout = *itin1;
		++itout;
	}
	if (itout!=last)
		seq.erase(itout,last);

	if (needs_further_processing) {
		// Clear seq and start over.
		epvector v = std::move(seq);
		construct_from_epvector(std::move(v));
	}
}

/** Check if this expairseq is in sorted (canonical) form.  Useful mainly for
 *  debugging or in assertions since being sorted is an invariance. */
bool expairseq::is_canonical() const
{
	if (seq.size() <= 1)
		return 1;
	
	auto it = seq.begin(), itend = seq.end();
	auto it_last = it;
	for (++it; it!=itend; it_last=it, ++it) {
		if (!(it_last->is_less(*it) || it_last->is_equal(*it))) {
			if (!is_exactly_a<numeric>(it_last->rest) ||
				!is_exactly_a<numeric>(it->rest)) {
				// double test makes it easier to set a breakpoint...
				if (!is_exactly_a<numeric>(it_last->rest) ||
					!is_exactly_a<numeric>(it->rest)) {
					printpair(std::clog, *it_last, 0);
					std::clog << ">";
					printpair(std::clog, *it, 0);
					std::clog << "\n";
					std::clog << "pair1:" << std::endl;
					it_last->rest.print(print_tree(std::clog));
					it_last->coeff.print(print_tree(std::clog));
					std::clog << "pair2:" << std::endl;
					it->rest.print(print_tree(std::clog));
					it->coeff.print(print_tree(std::clog));
					return 0;
				}
			}
		}
	}
	return 1;
}

/** Member-wise expand the expairs in this sequence.
 *
 *  @see expairseq::expand()
 *  @return epvector containing expanded pairs, empty if no members
 *    had to be changed. */
epvector expairseq::expandchildren(unsigned options) const
{
	auto cit = seq.begin(), last = seq.end();
	while (cit!=last) {
		const ex expanded_ex = cit->rest.expand(options);
		if (!are_ex_trivially_equal(cit->rest,expanded_ex)) {
			
			// something changed, copy seq, eval and return it
			epvector s;
			s.reserve(seq.size());
			
			// copy parts of seq which are known not to have changed
			s.insert(s.begin(), seq.begin(), cit);

			// copy first changed element
			s.push_back(expair(expanded_ex, cit->coeff));
			++cit;

			// copy rest
			while (cit != last) {
				s.push_back(expair(cit->rest.expand(options), cit->coeff));
				++cit;
			}
			return s;
		}

		++cit;
	}
	
	return epvector(); // empty signalling nothing has changed
}


/** Member-wise evaluate the expairs in this sequence.
 *
 *  @see expairseq::eval()
 *  @return epvector containing evaluated pairs, empty if no members
 *    had to be changed. */
epvector expairseq::evalchildren() const
{
	auto cit = seq.begin(), last = seq.end();
	while (cit!=last) {
		const expair evaled_pair = combine_ex_with_coeff_to_pair(cit->rest, cit->coeff);
		if (unlikely(!evaled_pair.is_equal(*cit))) {

			// something changed: copy seq, eval and return it
			epvector s;
			s.reserve(seq.size());

			// copy parts of seq which are known not to have changed
			s.insert(s.begin(), seq.begin(), cit);

			// copy first changed element
			s.push_back(evaled_pair);
			++cit;

			// copy rest
			while (cit != last) {
				s.push_back(combine_ex_with_coeff_to_pair(cit->rest, cit->coeff));
				++cit;
			}
			return s;
		}

		++cit;
	}

	return epvector(); // signalling nothing has changed
}

/** Member-wise substitute in this sequence.
 *
 *  @see expairseq::subs()
 *  @return epvector containing expanded pairs, empty if no members
 *    had to be changed. */
epvector expairseq::subschildren(const exmap & m, unsigned options) const
{
	// When any of the objects to be substituted is a product or power
	// we have to recombine the pairs because the numeric coefficients may
	// be part of the search pattern.
	if (!(options & (subs_options::pattern_is_product | subs_options::pattern_is_not_product))) {

		// Search the list of substitutions and cache our findings
		for (auto & it : m) {
			if (is_exactly_a<mul>(it.first) || is_exactly_a<power>(it.first)) {
				options |= subs_options::pattern_is_product;
				break;
			}
		}
		if (!(options & subs_options::pattern_is_product))
			options |= subs_options::pattern_is_not_product;
	}

	if (options & subs_options::pattern_is_product) {

		// Substitute in the recombined pairs
		auto cit = seq.begin(), last = seq.end();
		while (cit != last) {

			const ex &orig_ex = recombine_pair_to_ex(*cit);
			const ex &subsed_ex = orig_ex.subs(m, options);
			if (!are_ex_trivially_equal(orig_ex, subsed_ex)) {

				// Something changed: copy seq, subs and return it
				epvector s;
				s.reserve(seq.size());

				// Copy parts of seq which are known not to have changed
				s.insert(s.begin(), seq.begin(), cit);

				// Copy first changed element
				s.push_back(split_ex_to_pair(subsed_ex));
				++cit;

				// Copy rest
				while (cit != last) {
					s.push_back(split_ex_to_pair(recombine_pair_to_ex(*cit).subs(m, options)));
					++cit;
				}
				return s;
			}

			++cit;
		}

	} else {

		// Substitute only in the "rest" part of the pairs
		auto cit = seq.begin(), last = seq.end();
		while (cit != last) {

			const ex subsed_rest = cit->rest.subs(m, options);
			const expair subsed_pair = combine_ex_with_coeff_to_pair(subsed_rest, cit->coeff);
			if (!subsed_pair.is_equal(*cit)) {
			
				// Something changed: copy seq, subs and return it
				epvector s;
				s.reserve(seq.size());

				// Copy parts of seq which are known not to have changed
				s.insert(s.begin(), seq.begin(), cit);
			
				// Copy first changed element
				s.push_back(subsed_pair);
				++cit;

				// Copy rest
				while (cit != last) {
					s.push_back(combine_ex_with_coeff_to_pair(cit->rest.subs(m, options), cit->coeff));
					++cit;
				}
				return s;
			}

			++cit;
		}
	}
	
	// Nothing has changed
	return epvector();
}

//////////
// static member variables
//////////

} // namespace GiNaC
extern template GiNaC::registered_class_info GiNaC::container<std::vector>::reg_info;
extern template GiNaC::registered_class_info GiNaC::container<std::list>::reg_info;
