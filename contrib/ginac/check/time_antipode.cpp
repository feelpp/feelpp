/** @file time_antipode.cpp
 *
 *  This is a beautiful example that calculates the counterterm for the
 *  overall divergence of some special sorts of Feynman diagrams in a massless
 *  Yukawa theory.  For this end it computes the antipode of the corresponding
 *  decorated rooted tree using dimensional regularization in the parameter
 *  x==-(D-4)/2, which leads to a Laurent series in x.  The renormalization
 *  scheme used is the minimal subtraction scheme (MS).  From an efficiency
 *  point of view it boils down to power series expansion.  It also has quite
 *  an intriguing check for consistency, which is why we include it here.
 *
 *  This program is based on work by
 *      Isabella Bierenbaum <bierenbaum@thep.physik.uni-mainz.de> and
 *      Dirk Kreimer <dkreimer@bu.edu>.
 *  For details, please see <http://www.arXiv.org/abs/hep-th/0111192>.
 */

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

#include "ginac.h"
#include "timer.h"
using namespace GiNaC;

#include <map>
#include <set>
#include <stdexcept>
#include <typeinfo>
#include <utility>
#include <vector>
using namespace std;

// whether to run this beast or not:
static const bool do_test = true;

// regularization parameter:
static const symbol x("x");

typedef pair<unsigned, unsigned> ijpair;
typedef pair<class node, bool> child;

const constant TrOne("Tr[One]", numeric(4));

/* Extract only the divergent part of a series and discard the rest. */
static ex div_part(const ex &exarg, const symbol &x, unsigned grad)
{
	const ex exser = exarg.series(x==0, grad);
	if (exser.degree(x)<0)
		throw runtime_error("divergent part truncation disaster");
	ex exser_trunc;
	for (int i=exser.ldegree(x); i<0; ++i)
		exser_trunc += exser.coeff(x,i)*pow(x,i);
	// NB: exser_trunc is by construction collected in x.
	return exser_trunc;
}

/* F_ab(a, i, b, j, "x") is a common pattern in all vertex evaluators. */
static ex F_ab(int a, int i, int b, int j, const symbol &x)
{
	using GiNaC::tgamma;
	if ((i==0 && a<=0) || (j==0 && b<=0))
		return 0;
	else
		return (tgamma(2-a-(i+1)*x)*
		        tgamma(2-b-(1+j)*x)*
		        tgamma(a+b-2+(1+i+j)*x)/
		        tgamma(a+i*x)/
		        tgamma(b+j*x)/tgamma(4-a-b-(2+i+j)*x));
}

/* Abstract base class (ABC) for all types of vertices. */
class vertex {
public:
	vertex(ijpair ij = ijpair(0,0)) : indices(ij) { }
	void increment_indices(const ijpair &ind) { indices.first += ind.first; indices.second += ind.second; }
	virtual ~vertex() { }
	virtual vertex* copy() const = 0;
	virtual ijpair get_increment() const { return indices; }
	virtual const ex evaluate(const symbol &x, const unsigned grad) const = 0;
	bool operator==(const vertex &v) const { return (indices==v.indices); }
	bool operator<(const vertex &v) const { return (indices<v.indices); }
protected:
	ijpair indices;
};


/*
 * Class of vertices of type Sigma.
 */
class Sigma : public vertex {
public:
	Sigma(ijpair ij = ijpair(0,0)) : vertex(ij) { }
	vertex* copy() const { return new Sigma(*this); }
	ijpair get_increment() const { return ijpair(indices.first+indices.second+1, 0); }
	const ex evaluate(const symbol &x, const unsigned grad) const;
private:
};

const ex Sigma::evaluate(const symbol &x, const unsigned grad) const
{
	// c.f. comments in node::evaluate()
	static map<Sigma,ex> catalog;
	static unsigned prev_grad = 0;
	if (grad>prev_grad) {
		catalog.clear();
		prev_grad = grad;
	}
	map<Sigma,ex>::iterator pos = catalog.find(*this);
	if (pos==catalog.end()) {
		int i = indices.first;
		int j = indices.second;
		const ex result = ((F_ab(0,i,1,j,x)+F_ab(1,i,1,j,x)-F_ab(1,i,0,j,x))/2).series(x==0, grad).expand();
		pos = catalog.insert(map<Sigma,ex>::value_type(*this,result)).first;
		if (grad<prev_grad)
			prev_grad = grad;
	}
	return pos->second;
}


/** Class of vertices of type Sigma_flipped, sitting in the upper fermionline of Vacuum; no consequences for Gamma. */
class Sigma_flipped : public Sigma {
public:
	Sigma_flipped(ijpair ij = ijpair(0,0)) : Sigma(ij) { }
	vertex* copy() const { return new Sigma_flipped(*this); }
	ijpair get_increment() const { return ijpair(0, indices.first+indices.second+1); }
	const ex evaluate(const symbol &x, const unsigned grad) const { return Sigma::evaluate(x, grad); }
private:
};


/*
 *Class of vertices of type Gamma.
 */
class Gamma : public vertex {
public:
	Gamma(ijpair ij = ijpair(0,0)) : vertex(ij) { }
	vertex* copy() const { return new Gamma(*this); }
  	ijpair get_increment() const { return ijpair(indices.first+indices.second+1, 0); }
	const ex evaluate(const symbol &x, const unsigned grad) const;
private:
};

const ex Gamma::evaluate(const symbol &x, const unsigned grad) const
{
	// c.f. comments in node::evaluate()
	static map<Gamma,ex> catalog;
	static unsigned prev_grad = 0;
	if (prev_grad<grad) {
		catalog.clear();
		prev_grad = grad;
	}
	map<Gamma,ex>::iterator pos = catalog.find(*this);
	if (pos==catalog.end()) {
		int i = indices.first;
		int j = indices.second;
		const ex result = (F_ab(1,i,1,j,x)).series(x==0,grad).expand();
		pos = catalog.insert(map<Gamma,ex>::value_type(*this,result)).first;
		if (grad<prev_grad)
			prev_grad = grad;
	}
	return pos->second;
}


/*
 * Class of vertices of type Vacuum.
 */
class Vacuum : public vertex {
public:
	Vacuum(ijpair ij = ijpair(0,0)) : vertex(ij) { }
	vertex* copy() const { return new Vacuum(*this); }
	ijpair get_increment() const { return ijpair(0, indices.first+indices.second+1); }
	const ex evaluate(const symbol &x, const unsigned grad) const;
private:
};

const ex Vacuum::evaluate(const symbol &x, const unsigned grad) const
{
	// c.f. comments in node::evaluate()
	static map<Vacuum,ex> catalog;
	static unsigned prev_grad = 0;
	if (prev_grad<grad) {
		catalog.clear();
		prev_grad = grad;
	}
	map<Vacuum,ex>::iterator pos = catalog.find(*this);
	if (pos==catalog.end()) {
		int i = indices.first;
		int j = indices.second;
		const ex result = ((-TrOne*(F_ab(0,i,1,j,x)-F_ab(1,i,1,j,x)+F_ab(1,i,0,j,x)))/2).series(x==0,grad).expand();
		pos = catalog.insert(map<Vacuum,ex>::value_type(*this,result)).first;
		if (grad<prev_grad)
			prev_grad = grad;
	}
	return pos->second;
}


/*
 * Class of nodes (or trees or subtrees), including list of children.
 */
class node {
public:
	node(const vertex &v) { vert = v.copy(); }
	node(const node &n) { vert = (n.vert)->copy(); children = n.children; }
	const node & operator=(const node &);
	~node() { delete vert; }
	void add_child(const node &, bool = false);
	const ex evaluate(const symbol &x, unsigned grad) const;
	unsigned total_edges() const;
	bool operator==(const node &) const;
	bool operator<(const node &) const;
private:
	vertex *vert;
	multiset<child> children;
};

const node & node::operator=(const node &n)
{
	if (this!=&n) {
		delete vert;
		vert = (n.vert)->copy();
		children = n.children;
	}
	return *this;
}

void node::add_child(const node &childnode, bool cut)
{
	children.insert(child(childnode, cut));
	if(!cut)
		vert->increment_indices(childnode.vert->get_increment());
}

const ex node::evaluate(const symbol &x, unsigned grad) const
{
	static map<node,ex> catalog;    // lookup table for already evaluated subnodes
	static unsigned prev_grad = 0;  // grad argument that the catalog has been build for
	if (grad>prev_grad) {
		// We haven't computed far enough last time.  Our catalog cannot cope with
		// the demands for this value of grad so let's clear it.
		catalog.clear();
		prev_grad = grad;
	}
	ex product = 1;   // accumulator for all the children
	for (multiset<child>::const_iterator i=children.begin(); i!=children.end(); ++i) {
		map<node,ex>::iterator pos = catalog.find(i->first);
		if (pos==catalog.end()) {
			pos = catalog.insert(map<node,ex>::value_type(i->first,i->first.evaluate(x,grad).series(x==0,grad).expand())).first;
			if (grad<prev_grad) {
				// We have just spoiled the catalog by inserting a series computed
				// to lower grad than the others in it.  So let's make sure next time
				// we don't use one of the newly inserted ones by bumping prev_grad
				// down to the current value of grad.
				prev_grad = grad;
			}
		}
		if (!i->second)
			product *= pos->second;
		else
			product *= -div_part(pos->second,x,grad);
	}
	return (product * vert->evaluate(x,grad));
}

unsigned node::total_edges() const
{
	unsigned accu = 0;
	for (multiset<child>::const_iterator i=children.begin(); i!=children.end(); ++i) {
		accu += i->first.total_edges();
		++accu;
	}
	return accu;
}

bool node::operator==(const node &n) const
{
	// Are the types of the top-level node vertices equal?
	if (typeid(*vert)!=typeid(*n.vert))
		return false;
	// Are the indices of the top-level nodes equal?
	if (!(*vert==*n.vert))
		return false;
	// Are the sets of children equal, one by one?
	return (children==n.children);
}

bool node::operator<(const node &n) const
{
	// Are the types of the top-level node vertices different?
	if (typeid(*vert)!=typeid(*n.vert))
		return typeid(*vert).before(typeid(*n.vert));
	// Are the indices of the top-level nodes different?
	if (!(*vert==*n.vert))
		return (vert<n.vert);
	// Are the sets of children different, one by one?
	return (children<n.children);
}

/*
 * These operators let us write down trees in an intuitive way, by adding
 * arbitrarily complex children to a given vertex.  The eye candy that can be
 * produced with it makes detection of errors much simpler than with code
 * written using calls to node's method add_child() because it allows for
 * editor-assisted indentation.
 */
const node operator+(const node &n, const child &c)
{
	node nn(n);
	nn.add_child(c.first, c.second);
	return nn;
}

void operator+=(node &n, const child &c)
{
	n.add_child(c.first, c.second);
}


/*                Gamma
 *                  |
 *                Gamma
 */
static const node tree1(unsigned cuts=0)
{
	return (Gamma()
	        + child(Gamma(),
	                bool(cuts & 1)));
}

/*                Gamma
 *               /  |  \
 *        Vacuum  Gamma  Vacuum
 *       /   |  \
 *  Sigma Sigma Sigma0
 */
static const node tree2(unsigned cuts=0)
{
	return (Gamma()
	        + child(Vacuum()
	                + child(Sigma(), bool(cuts & 1))
	                + child(Sigma(), bool(cuts & 2))
	                + child(Sigma_flipped(), bool(cuts & 4)),
	                bool(cuts & 8))
	        + child(Gamma(), bool(cuts & 16))
	        + child(Gamma(), bool(cuts & 32)));
}

/*                Gamma
 *                  |
 *                Gamma
 *                  |
 *                Gamma
 *                /   \
 *           Vacuum    Gamma
 *           /    \       \
 *        Sigma  Sigma   Sigma
 */
static const node tree3(unsigned cuts=0)
{
	return (Gamma()
	        + child(Gamma()
	                + child(Gamma()
	                        + child(Gamma()
	                                + child(Sigma(), bool(cuts & 1)),
	                                bool(cuts & 2))
	                        + child(Vacuum()
	                                + child(Sigma(), bool(cuts & 4))
	                                + child(Sigma(), bool(cuts & 8)),
	                        bool(cuts & 16)),
	                bool(cuts & 32)),
	        bool(cuts & 64)));
}

/*                Gamma
 *                /   \
 *           Sigma     Vacuum
 *          /   \       /   \
 *      Sigma Sigma  Sigma0 Sigma
 */
static const node tree4(unsigned cuts=0)
{
	return (Gamma()
	        + child(Sigma()
	                + child(Sigma(), bool(cuts & 1))
	                + child(Sigma(), bool(cuts & 2)),
	                bool(cuts & 4))
	        + child(Vacuum()
	                + child(Sigma_flipped(), bool(cuts & 8))
	                + child(Sigma(), bool(cuts & 16)),
	                bool(cuts & 32)));
}

/*                Sigma
 *               /  |  \
 *         Sigma Vacuum  Vacuum
 *               /    \       \
 *            Sigma  Sigma0   Sigma
 */
static const node tree5(unsigned cuts=0)
{
	return (Sigma()
	        + child(Sigma(), bool(cuts & 1))
	        + child(Vacuum()
	                + child(Sigma(), bool(cuts & 2))
	                + child(Sigma_flipped(), bool(cuts & 4)),
	                bool(cuts & 8))
	        + child(Vacuum()
	                + child(Sigma(), bool(cuts & 16)),
	                bool(cuts & 32)));
}

/*               Vacuum
 *               /    \
 *           Sigma    Sigma0
 *             |        |
 *           Sigma    Sigma
 *                      |
 *                    Vacuum
 */
static const node tree6(unsigned cuts=0)
{
	return (Vacuum()
	        + child(Sigma()
	                + child(Sigma(), bool(cuts & 1)),
	                bool(cuts & 2))
	        + child(Sigma_flipped()
	                + child(Sigma()
	                        + child(Vacuum(), bool(cuts & 4)),
	                        bool(cuts & 8)),
	                bool(cuts & 16)));
}

static unsigned test_tree(const node tree_generator(unsigned))
{
	const int edges = tree_generator(0).total_edges();
   	const int vertices = edges+1;
	
	// fill a vector of all possible 2^edges combinations of cuts...
	vector<node> counter;
	for (unsigned i=0; i<(1U<<edges); ++i)
		counter.push_back(tree_generator(i));
	
	// ...the sum, when evaluated and reexpanded, is the antipode...
	ex result = 0;
	for (vector<node>::iterator i=counter.begin(); i!=counter.end(); ++i)
		result = (result+i->evaluate(x,vertices-1)).series(x==0,vertices-1).expand();
	
	// ...and has the nice property that in each term all the Eulers cancel:
	if (result.has(Euler)) {
		clog << "The antipode was miscalculated\nAntipode==" << result
		     << "\nshould not have any occurrence of Euler" << endl;
		return 1;
	} else if (result.ldegree(x)!=-vertices || result.degree(x)!=0) {
		clog << "The antipode was miscalculated\nAntipode==" << result
		     << "\nshould run from " << x << "^(" << -vertices << ") to "
		     << x << "^(0)" << "but it runs from " << x << "^("
		     << result.ldegree(x) << ")" << "to " << x << "^("
		     << result.degree(x) << ")" << endl;
		return 1;
	}
	return 0;
}

unsigned time_antipode()
{
	unsigned result = 0;
	timer jaeger_le_coultre;
	
	cout << "timing computation of antipodes in Yukawa theory" << flush;
	
	if (do_test) {
		jaeger_le_coultre.start();
		result += test_tree(tree1);  cout << '.' << flush;
		result += test_tree(tree2);  cout << '.' << flush;
		result += test_tree(tree3);  cout << '.' << flush;
		result += test_tree(tree4);  cout << '.' << flush;
		result += test_tree(tree5);  cout << '.' << flush;
		result += test_tree(tree6);  cout << '.' << flush;
		
		cout << jaeger_le_coultre.read() << "s (total)" << endl;
	} else {
		cout << " disabled" << endl;
	}
	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_antipode();
}
