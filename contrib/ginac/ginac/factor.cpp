/** @file factor.cpp
 *
 *  Polynomial factorization (implementation).
 *
 *  The interface function factor() at the end of this file is defined in the
 *  GiNaC namespace. All other utility functions and classes are defined in an
 *  additional anonymous namespace.
 *
 *  Factorization starts by doing a square free factorization and making the
 *  coefficients integer. Then, depending on the number of free variables it
 *  proceeds either in dedicated univariate or multivariate factorization code.
 *
 *  Univariate factorization does a modular factorization via Berlekamp's
 *  algorithm and distinct degree factorization. Hensel lifting is used at the
 *  end.
 *  
 *  Multivariate factorization uses the univariate factorization (applying a
 *  evaluation homomorphism first) and Hensel lifting raises the answer to the
 *  multivariate domain. The Hensel lifting code is completely distinct from the
 *  code used by the univariate factorization.
 *
 *  Algorithms used can be found in
 *    [Wan] An Improved Multivariate Polynomial Factoring Algorithm,
 *          P.S.Wang,
 *          Mathematics of Computation, Vol. 32, No. 144 (1978) 1215--1231.
 *    [GCL] Algorithms for Computer Algebra,
 *          K.O.Geddes, S.R.Czapor, G.Labahn,
 *          Springer Verlag, 1992.
 *    [Mig] Some Useful Bounds,
 *          M.Mignotte, 
 *          In "Computer Algebra, Symbolic and Algebraic Computation" (B.Buchberger et al., eds.),
 *          pp. 259-263, Springer-Verlag, New York, 1982.
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

//#define DEBUGFACTOR

#include "factor.h"

#include "ex.h"
#include "numeric.h"
#include "operators.h"
#include "inifcns.h"
#include "symbol.h"
#include "relational.h"
#include "power.h"
#include "mul.h"
#include "normal.h"
#include "add.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <vector>
#ifdef DEBUGFACTOR
#include <ostream>
#endif
using namespace std;

#include <cln/cln.h>
using namespace cln;

namespace GiNaC {

#ifdef DEBUGFACTOR
#define DCOUT(str) cout << #str << endl
#define DCOUTVAR(var) cout << #var << ": " << var << endl
#define DCOUT2(str,var) cout << #str << ": " << var << endl
ostream& operator<<(ostream& o, const vector<int>& v)
{
	vector<int>::const_iterator i = v.begin(), end = v.end();
	while ( i != end ) {
		o << *i++ << " ";
	}
	return o;
}
static ostream& operator<<(ostream& o, const vector<cl_I>& v)
{
	vector<cl_I>::const_iterator i = v.begin(), end = v.end();
	while ( i != end ) {
		o << *i << "[" << i-v.begin() << "]" << " ";
		++i;
	}
	return o;
}
static ostream& operator<<(ostream& o, const vector<cl_MI>& v)
{
	vector<cl_MI>::const_iterator i = v.begin(), end = v.end();
	while ( i != end ) {
		o << *i << "[" << i-v.begin() << "]" << " ";
		++i;
	}
	return o;
}
ostream& operator<<(ostream& o, const vector<numeric>& v)
{
	for ( size_t i=0; i<v.size(); ++i ) {
		o << v[i] << " ";
	}
	return o;
}
ostream& operator<<(ostream& o, const vector< vector<cl_MI> >& v)
{
	vector< vector<cl_MI> >::const_iterator i = v.begin(), end = v.end();
	while ( i != end ) {
		o << i-v.begin() << ": " << *i << endl;
		++i;
	}
	return o;
}
#else
#define DCOUT(str)
#define DCOUTVAR(var)
#define DCOUT2(str,var)
#endif // def DEBUGFACTOR

// anonymous namespace to hide all utility functions
namespace {

////////////////////////////////////////////////////////////////////////////////
// modular univariate polynomial code

typedef std::vector<cln::cl_MI> umodpoly;
typedef std::vector<cln::cl_I> upoly;
typedef vector<umodpoly> upvec;

// COPY FROM UPOLY.HPP

// CHANGED size_t -> int !!!
template<typename T> static int degree(const T& p)
{
	return p.size() - 1;
}

template<typename T> static typename T::value_type lcoeff(const T& p)
{
	return p[p.size() - 1];
}

static bool normalize_in_field(umodpoly& a)
{
	if (a.size() == 0)
		return true;
	if ( lcoeff(a) == a[0].ring()->one() ) {
		return true;
	}

	const cln::cl_MI lc_1 = recip(lcoeff(a));
	for (std::size_t k = a.size(); k-- != 0; )
		a[k] = a[k]*lc_1;
	return false;
}

template<typename T> static void
canonicalize(T& p, const typename T::size_type hint = std::numeric_limits<typename T::size_type>::max())
{
	if (p.empty())
		return;

	std::size_t i = p.size() - 1;
	// Be fast if the polynomial is already canonicalized
	if (!zerop(p[i]))
		return;

	if (hint < p.size())
		i = hint;

	bool is_zero = false;
	do {
		if (!zerop(p[i])) {
			++i;
			break;
		}
		if (i == 0) {
			is_zero = true;
			break;
		}
		--i;
	} while (true);

	if (is_zero) {
		p.clear();
		return;
	}

	p.erase(p.begin() + i, p.end());
}

// END COPY FROM UPOLY.HPP

static void expt_pos(umodpoly& a, unsigned int q)
{
	if ( a.empty() ) return;
	cl_MI zero = a[0].ring()->zero(); 
	int deg = degree(a);
	a.resize(degree(a)*q+1, zero);
	for ( int i=deg; i>0; --i ) {
		a[i*q] = a[i];
		a[i] = zero;
	}
}

template<bool COND, typename T = void> struct enable_if
{
	typedef T type;
};

template<typename T> struct enable_if<false, T> { /* empty */ };

template<typename T> struct uvar_poly_p
{
	static const bool value = false;
};

template<> struct uvar_poly_p<upoly>
{
	static const bool value = true;
};

template<> struct uvar_poly_p<umodpoly>
{
	static const bool value = true;
};

template<typename T>
// Don't define this for anything but univariate polynomials.
static typename enable_if<uvar_poly_p<T>::value, T>::type
operator+(const T& a, const T& b)
{
	int sa = a.size();
	int sb = b.size();
	if ( sa >= sb ) {
		T r(sa);
		int i = 0;
		for ( ; i<sb; ++i ) {
			r[i] = a[i] + b[i];
		}
		for ( ; i<sa; ++i ) {
			r[i] = a[i];
		}
		canonicalize(r);
		return r;
	}
	else {
		T r(sb);
		int i = 0;
		for ( ; i<sa; ++i ) {
			r[i] = a[i] + b[i];
		}
		for ( ; i<sb; ++i ) {
			r[i] = b[i];
		}
		canonicalize(r);
		return r;
	}
}

template<typename T>
// Don't define this for anything but univariate polynomials. Otherwise
// overload resolution might fail (this actually happens when compiling
// GiNaC with g++ 3.4).
static typename enable_if<uvar_poly_p<T>::value, T>::type
operator-(const T& a, const T& b)
{
	int sa = a.size();
	int sb = b.size();
	if ( sa >= sb ) {
		T r(sa);
		int i = 0;
		for ( ; i<sb; ++i ) {
			r[i] = a[i] - b[i];
		}
		for ( ; i<sa; ++i ) {
			r[i] = a[i];
		}
		canonicalize(r);
		return r;
	}
	else {
		T r(sb);
		int i = 0;
		for ( ; i<sa; ++i ) {
			r[i] = a[i] - b[i];
		}
		for ( ; i<sb; ++i ) {
			r[i] = -b[i];
		}
		canonicalize(r);
		return r;
	}
}

static upoly operator*(const upoly& a, const upoly& b)
{
	upoly c;
	if ( a.empty() || b.empty() ) return c;

	int n = degree(a) + degree(b);
	c.resize(n+1, 0);
	for ( int i=0 ; i<=n; ++i ) {
		for ( int j=0 ; j<=i; ++j ) {
			if ( j > degree(a) || (i-j) > degree(b) ) continue;
			c[i] = c[i] + a[j] * b[i-j];
		}
	}
	canonicalize(c);
	return c;
}

static umodpoly operator*(const umodpoly& a, const umodpoly& b)
{
	umodpoly c;
	if ( a.empty() || b.empty() ) return c;

	int n = degree(a) + degree(b);
	c.resize(n+1, a[0].ring()->zero());
	for ( int i=0 ; i<=n; ++i ) {
		for ( int j=0 ; j<=i; ++j ) {
			if ( j > degree(a) || (i-j) > degree(b) ) continue;
			c[i] = c[i] + a[j] * b[i-j];
		}
	}
	canonicalize(c);
	return c;
}

static upoly operator*(const upoly& a, const cl_I& x)
{
	if ( zerop(x) ) {
		upoly r;
		return r;
	}
	upoly r(a.size());
	for ( size_t i=0; i<a.size(); ++i ) {
		r[i] = a[i] * x;
	}
	return r;
}

static upoly operator/(const upoly& a, const cl_I& x)
{
	if ( zerop(x) ) {
		upoly r;
		return r;
	}
	upoly r(a.size());
	for ( size_t i=0; i<a.size(); ++i ) {
		r[i] = exquo(a[i],x);
	}
	return r;
}

static umodpoly operator*(const umodpoly& a, const cl_MI& x)
{
	umodpoly r(a.size());
	for ( size_t i=0; i<a.size(); ++i ) {
		r[i] = a[i] * x;
	}
	canonicalize(r);
	return r;
}

static void upoly_from_ex(upoly& up, const ex& e, const ex& x)
{
	// assert: e is in Z[x]
	int deg = e.degree(x);
	up.resize(deg+1);
	int ldeg = e.ldegree(x);
	for ( ; deg>=ldeg; --deg ) {
		up[deg] = the<cl_I>(ex_to<numeric>(e.coeff(x, deg)).to_cl_N());
	}
	for ( ; deg>=0; --deg ) {
		up[deg] = 0;
	}
	canonicalize(up);
}

static void umodpoly_from_upoly(umodpoly& ump, const upoly& e, const cl_modint_ring& R)
{
	int deg = degree(e);
	ump.resize(deg+1);
	for ( ; deg>=0; --deg ) {
		ump[deg] = R->canonhom(e[deg]);
	}
	canonicalize(ump);
}

static void umodpoly_from_ex(umodpoly& ump, const ex& e, const ex& x, const cl_modint_ring& R)
{
	// assert: e is in Z[x]
	int deg = e.degree(x);
	ump.resize(deg+1);
	int ldeg = e.ldegree(x);
	for ( ; deg>=ldeg; --deg ) {
		cl_I coeff = the<cl_I>(ex_to<numeric>(e.coeff(x, deg)).to_cl_N());
		ump[deg] = R->canonhom(coeff);
	}
	for ( ; deg>=0; --deg ) {
		ump[deg] = R->zero();
	}
	canonicalize(ump);
}

#ifdef DEBUGFACTOR
static void umodpoly_from_ex(umodpoly& ump, const ex& e, const ex& x, const cl_I& modulus)
{
	umodpoly_from_ex(ump, e, x, find_modint_ring(modulus));
}
#endif

static ex upoly_to_ex(const upoly& a, const ex& x)
{
	if ( a.empty() ) return 0;
	ex e;
	for ( int i=degree(a); i>=0; --i ) {
		e += numeric(a[i]) * pow(x, i);
	}
	return e;
}

static ex umodpoly_to_ex(const umodpoly& a, const ex& x)
{
	if ( a.empty() ) return 0;
	cl_modint_ring R = a[0].ring();
	cl_I mod = R->modulus;
	cl_I halfmod = (mod-1) >> 1;
	ex e;
	for ( int i=degree(a); i>=0; --i ) {
		cl_I n = R->retract(a[i]);
		if ( n > halfmod ) {
			e += numeric(n-mod) * pow(x, i);
		} else {
			e += numeric(n) * pow(x, i);
		}
	}
	return e;
}

static upoly umodpoly_to_upoly(const umodpoly& a)
{
	upoly e(a.size());
	if ( a.empty() ) return e;
	cl_modint_ring R = a[0].ring();
	cl_I mod = R->modulus;
	cl_I halfmod = (mod-1) >> 1;
	for ( int i=degree(a); i>=0; --i ) {
		cl_I n = R->retract(a[i]);
		if ( n > halfmod ) {
			e[i] = n-mod;
		} else {
			e[i] = n;
		}
	}
	return e;
}

static umodpoly umodpoly_to_umodpoly(const umodpoly& a, const cl_modint_ring& R, unsigned int m)
{
	umodpoly e;
	if ( a.empty() ) return e;
	cl_modint_ring oldR = a[0].ring();
	size_t sa = a.size();
	e.resize(sa+m, R->zero());
	for ( size_t i=0; i<sa; ++i ) {
		e[i+m] = R->canonhom(oldR->retract(a[i]));
	}
	canonicalize(e);
	return e;
}

/** Divides all coefficients of the polynomial a by the integer x.
 *  All coefficients are supposed to be divisible by x. If they are not, the
 *  the<cl_I> cast will raise an exception.
 *
 *  @param[in,out] a  polynomial of which the coefficients will be reduced by x
 *  @param[in]     x  integer that divides the coefficients
 */
static void reduce_coeff(umodpoly& a, const cl_I& x)
{
	if ( a.empty() ) return;

	cl_modint_ring R = a[0].ring();
	umodpoly::iterator i = a.begin(), end = a.end();
	for ( ; i!=end; ++i ) {
		// cln cannot perform this division in the modular field
		cl_I c = R->retract(*i);
		*i = cl_MI(R, the<cl_I>(c / x));
	}
}

/** Calculates remainder of a/b.
 *  Assertion: a and b not empty.
 *
 *  @param[in]  a  polynomial dividend
 *  @param[in]  b  polynomial divisor
 *  @param[out] r  polynomial remainder
 */
static void rem(const umodpoly& a, const umodpoly& b, umodpoly& r)
{
	int k, n;
	n = degree(b);
	k = degree(a) - n;
	r = a;
	if ( k < 0 ) return;

	do {
		cl_MI qk = div(r[n+k], b[n]);
		if ( !zerop(qk) ) {
			for ( int i=0; i<n; ++i ) {
				unsigned int j = n + k - 1 - i;
				r[j] = r[j] - qk * b[j-k];
			}
		}
	} while ( k-- );

	fill(r.begin()+n, r.end(), a[0].ring()->zero());
	canonicalize(r);
}

/** Calculates quotient of a/b.
 *  Assertion: a and b not empty.
 *
 *  @param[in]  a  polynomial dividend
 *  @param[in]  b  polynomial divisor
 *  @param[out] q  polynomial quotient
 */
static void div(const umodpoly& a, const umodpoly& b, umodpoly& q)
{
	int k, n;
	n = degree(b);
	k = degree(a) - n;
	q.clear();
	if ( k < 0 ) return;

	umodpoly r = a;
	q.resize(k+1, a[0].ring()->zero());
	do {
		cl_MI qk = div(r[n+k], b[n]);
		if ( !zerop(qk) ) {
			q[k] = qk;
			for ( int i=0; i<n; ++i ) {
				unsigned int j = n + k - 1 - i;
				r[j] = r[j] - qk * b[j-k];
			}
		}
	} while ( k-- );

	canonicalize(q);
}

/** Calculates quotient and remainder of a/b.
 *  Assertion: a and b not empty.
 *
 *  @param[in]  a  polynomial dividend
 *  @param[in]  b  polynomial divisor
 *  @param[out] r  polynomial remainder
 *  @param[out] q  polynomial quotient
 */
static void remdiv(const umodpoly& a, const umodpoly& b, umodpoly& r, umodpoly& q)
{
	int k, n;
	n = degree(b);
	k = degree(a) - n;
	q.clear();
	r = a;
	if ( k < 0 ) return;

	q.resize(k+1, a[0].ring()->zero());
	do {
		cl_MI qk = div(r[n+k], b[n]);
		if ( !zerop(qk) ) {
			q[k] = qk;
			for ( int i=0; i<n; ++i ) {
				unsigned int j = n + k - 1 - i;
				r[j] = r[j] - qk * b[j-k];
			}
		}
	} while ( k-- );

	fill(r.begin()+n, r.end(), a[0].ring()->zero());
	canonicalize(r);
	canonicalize(q);
}

/** Calculates the GCD of polynomial a and b.
 *
 *  @param[in]  a  polynomial
 *  @param[in]  b  polynomial
 *  @param[out] c  GCD
 */
static void gcd(const umodpoly& a, const umodpoly& b, umodpoly& c)
{
	if ( degree(a) < degree(b) ) return gcd(b, a, c);

	c = a;
	normalize_in_field(c);
	umodpoly d = b;
	normalize_in_field(d);
	umodpoly r;
	while ( !d.empty() ) {
		rem(c, d, r);
		c = d;
		d = r;
	}
	normalize_in_field(c);
}

/** Calculates the derivative of the polynomial a.
 *  
 *  @param[in]  a  polynomial of which to take the derivative
 *  @param[out] d  result/derivative
 */
static void deriv(const umodpoly& a, umodpoly& d)
{
	d.clear();
	if ( a.size() <= 1 ) return;

	d.insert(d.begin(), a.begin()+1, a.end());
	int max = d.size();
	for ( int i=1; i<max; ++i ) {
		d[i] = d[i] * (i+1);
	}
	canonicalize(d);
}

static bool unequal_one(const umodpoly& a)
{
	if ( a.empty() ) return true;
	return ( a.size() != 1 || a[0] != a[0].ring()->one() );
}

static bool equal_one(const umodpoly& a)
{
	return ( a.size() == 1 && a[0] == a[0].ring()->one() );
}

/** Returns true if polynomial a is square free.
 *
 *  @param[in] a  polynomial to check
 *  @return       true if polynomial is square free, false otherwise
 */
static bool squarefree(const umodpoly& a)
{
	umodpoly b;
	deriv(a, b);
	if ( b.empty() ) {
		return false;
	}
	umodpoly c;
	gcd(a, b, c);
	return equal_one(c);
}

// END modular univariate polynomial code
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// modular matrix

typedef vector<cl_MI> mvec;

class modular_matrix
{
	friend ostream& operator<<(ostream& o, const modular_matrix& m);
public:
	modular_matrix(size_t r_, size_t c_, const cl_MI& init) : r(r_), c(c_)
	{
		m.resize(c*r, init);
	}
	size_t rowsize() const { return r; }
	size_t colsize() const { return c; }
	cl_MI& operator()(size_t row, size_t col) { return m[row*c + col]; }
	cl_MI operator()(size_t row, size_t col) const { return m[row*c + col]; }
	void mul_col(size_t col, const cl_MI x)
	{
		for ( size_t rc=0; rc<r; ++rc ) {
			std::size_t i = c*rc + col;
			m[i] = m[i] * x;
		}
	}
	void sub_col(size_t col1, size_t col2, const cl_MI fac)
	{
		for ( size_t rc=0; rc<r; ++rc ) {
			std::size_t i1 = col1 + c*rc;
			std::size_t i2 = col2 + c*rc;
			m[i1] = m[i1] - m[i2]*fac;
		}
	}
	void switch_col(size_t col1, size_t col2)
	{
		for ( size_t rc=0; rc<r; ++rc ) {
			std::size_t i1 = col1 + rc*c;
			std::size_t i2 = col2 + rc*c;
			std::swap(m[i1], m[i2]);
		}
	}
	void mul_row(size_t row, const cl_MI x)
	{
		for ( size_t cc=0; cc<c; ++cc ) {
			std::size_t i = row*c + cc; 
			m[i] = m[i] * x;
		}
	}
	void sub_row(size_t row1, size_t row2, const cl_MI fac)
	{
		for ( size_t cc=0; cc<c; ++cc ) {
			std::size_t i1 = row1*c + cc;
			std::size_t i2 = row2*c + cc;
			m[i1] = m[i1] - m[i2]*fac;
		}
	}
	void switch_row(size_t row1, size_t row2)
	{
		for ( size_t cc=0; cc<c; ++cc ) {
			std::size_t i1 = row1*c + cc;
			std::size_t i2 = row2*c + cc;
			std::swap(m[i1], m[i2]);
		}
	}
	bool is_col_zero(size_t col) const
	{
		for ( size_t rr=0; rr<r; ++rr ) {
			std::size_t i = col + rr*c;
			if ( !zerop(m[i]) ) {
				return false;
			}
		}
		return true;
	}
	bool is_row_zero(size_t row) const
	{
		for ( size_t cc=0; cc<c; ++cc ) {
			std::size_t i = row*c + cc;
			if ( !zerop(m[i]) ) {
				return false;
			}
		}
		return true;
	}
	void set_row(size_t row, const vector<cl_MI>& newrow)
	{
		for (std::size_t i2 = 0; i2 < newrow.size(); ++i2) {
			std::size_t i1 = row*c + i2;
			m[i1] = newrow[i2];
		}
	}
	mvec::const_iterator row_begin(size_t row) const { return m.begin()+row*c; }
	mvec::const_iterator row_end(size_t row) const { return m.begin()+row*c+r; }
private:
	size_t r, c;
	mvec m;
};

#ifdef DEBUGFACTOR
modular_matrix operator*(const modular_matrix& m1, const modular_matrix& m2)
{
	const unsigned int r = m1.rowsize();
	const unsigned int c = m2.colsize();
	modular_matrix o(r,c,m1(0,0));

	for ( size_t i=0; i<r; ++i ) {
		for ( size_t j=0; j<c; ++j ) {
			cl_MI buf;
			buf = m1(i,0) * m2(0,j);
			for ( size_t k=1; k<c; ++k ) {
				buf = buf + m1(i,k)*m2(k,j);
			}
			o(i,j) = buf;
		}
	}
	return o;
}

ostream& operator<<(ostream& o, const modular_matrix& m)
{
	cl_modint_ring R = m(0,0).ring();
	o << "{";
	for ( size_t i=0; i<m.rowsize(); ++i ) {
		o << "{";
		for ( size_t j=0; j<m.colsize()-1; ++j ) {
			o << R->retract(m(i,j)) << ",";
		}
		o << R->retract(m(i,m.colsize()-1)) << "}";
		if ( i != m.rowsize()-1 ) {
			o << ",";
		}
	}
	o << "}";
	return o;
}
#endif // def DEBUGFACTOR

// END modular matrix
////////////////////////////////////////////////////////////////////////////////

/** Calculates the Q matrix for a polynomial. Used by Berlekamp's algorithm.
 *
 *  @param[in]  a_  modular polynomial
 *  @param[out] Q   Q matrix
 */
static void q_matrix(const umodpoly& a_, modular_matrix& Q)
{
	umodpoly a = a_;
	normalize_in_field(a);

	int n = degree(a);
	unsigned int q = cl_I_to_uint(a[0].ring()->modulus);
	umodpoly r(n, a[0].ring()->zero());
	r[0] = a[0].ring()->one();
	Q.set_row(0, r);
	unsigned int max = (n-1) * q;
	for ( size_t m=1; m<=max; ++m ) {
		cl_MI rn_1 = r.back();
		for ( size_t i=n-1; i>0; --i ) {
			r[i] = r[i-1] - (rn_1 * a[i]);
		}
		r[0] = -rn_1 * a[0];
		if ( (m % q) == 0 ) {
			Q.set_row(m/q, r);
		}
	}
}

/** Determine the nullspace of a matrix M-1.
 *
 *  @param[in,out] M      matrix, will be modified
 *  @param[out]    basis  calculated nullspace of M-1
 */
static void nullspace(modular_matrix& M, vector<mvec>& basis)
{
	const size_t n = M.rowsize();
	const cl_MI one = M(0,0).ring()->one();
	for ( size_t i=0; i<n; ++i ) {
		M(i,i) = M(i,i) - one;
	}
	for ( size_t r=0; r<n; ++r ) {
		size_t cc = 0;
		for ( ; cc<n; ++cc ) {
			if ( !zerop(M(r,cc)) ) {
				if ( cc < r ) {
					if ( !zerop(M(cc,cc)) ) {
						continue;
					}
					M.switch_col(cc, r);
				}
				else if ( cc > r ) {
					M.switch_col(cc, r);
				}
				break;
			}
		}
		if ( cc < n ) {
			M.mul_col(r, recip(M(r,r)));
			for ( cc=0; cc<n; ++cc ) {
				if ( cc != r ) {
					M.sub_col(cc, r, M(r,cc));
				}
			}
		}
	}

	for ( size_t i=0; i<n; ++i ) {
		M(i,i) = M(i,i) - one;
	}
	for ( size_t i=0; i<n; ++i ) {
		if ( !M.is_row_zero(i) ) {
			mvec nu(M.row_begin(i), M.row_end(i));
			basis.push_back(nu);
		}
	}
}

/** Berlekamp's modular factorization.
 *  
 *  The implementation follows the algorithm in chapter 8 of [GCL].
 *
 *  @param[in]  a    modular polynomial
 *  @param[out] upv  vector containing modular factors. if upv was not empty the
 *                   new elements are added at the end
 */
static void berlekamp(const umodpoly& a, upvec& upv)
{
	cl_modint_ring R = a[0].ring();
	umodpoly one(1, R->one());

	// find nullspace of Q matrix
	modular_matrix Q(degree(a), degree(a), R->zero());
	q_matrix(a, Q);
	vector<mvec> nu;
	nullspace(Q, nu);

	const unsigned int k = nu.size();
	if ( k == 1 ) {
		// irreducible
		return;
	}

	list<umodpoly> factors;
	factors.push_back(a);
	unsigned int size = 1;
	unsigned int r = 1;
	unsigned int q = cl_I_to_uint(R->modulus);

	list<umodpoly>::iterator u = factors.begin();

	// calculate all gcd's
	while ( true ) {
		for ( unsigned int s=0; s<q; ++s ) {
			umodpoly nur = nu[r];
			nur[0] = nur[0] - cl_MI(R, s);
			canonicalize(nur);
			umodpoly g;
			gcd(nur, *u, g);
			if ( unequal_one(g) && g != *u ) {
				umodpoly uo;
				div(*u, g, uo);
				if ( equal_one(uo) ) {
					throw logic_error("berlekamp: unexpected divisor.");
				}
				else {
					*u = uo;
				}
				factors.push_back(g);
				size = 0;
				list<umodpoly>::const_iterator i = factors.begin(), end = factors.end();
				while ( i != end ) {
					if ( degree(*i) ) ++size; 
					++i;
				}
				if ( size == k ) {
					list<umodpoly>::const_iterator i = factors.begin(), end = factors.end();
					while ( i != end ) {
						upv.push_back(*i++);
					}
					return;
				}
			}
		}
		if ( ++r == k ) {
			r = 1;
			++u;
		}
	}
}

// modular square free factorization is not used at the moment so we deactivate
// the code
#if 0

/** Calculates a^(1/prime).
 *  
 *  @param[in] a      polynomial
 *  @param[in] prime  prime number -> exponent 1/prime
 *  @param[in] ap     resulting polynomial
 */
static void expt_1_over_p(const umodpoly& a, unsigned int prime, umodpoly& ap)
{
	size_t newdeg = degree(a)/prime;
	ap.resize(newdeg+1);
	ap[0] = a[0];
	for ( size_t i=1; i<=newdeg; ++i ) {
		ap[i] = a[i*prime];
	}
}

/** Modular square free factorization.
 *
 *  @param[in]  a        polynomial
 *  @param[out] factors  modular factors
 *  @param[out] mult     corresponding multiplicities (exponents)
 */
static void modsqrfree(const umodpoly& a, upvec& factors, vector<int>& mult)
{
	const unsigned int prime = cl_I_to_uint(a[0].ring()->modulus);
	int i = 1;
	umodpoly b;
	deriv(a, b);
	if ( b.size() ) {
		umodpoly c;
		gcd(a, b, c);
		umodpoly w;
		div(a, c, w);
		while ( unequal_one(w) ) {
			umodpoly y;
			gcd(w, c, y);
			umodpoly z;
			div(w, y, z);
			factors.push_back(z);
			mult.push_back(i);
			++i;
			w = y;
			umodpoly buf;
			div(c, y, buf);
			c = buf;
		}
		if ( unequal_one(c) ) {
			umodpoly cp;
			expt_1_over_p(c, prime, cp);
			size_t previ = mult.size();
			modsqrfree(cp, factors, mult);
			for ( size_t i=previ; i<mult.size(); ++i ) {
				mult[i] *= prime;
			}
		}
	}
	else {
		umodpoly ap;
		expt_1_over_p(a, prime, ap);
		size_t previ = mult.size();
		modsqrfree(ap, factors, mult);
		for ( size_t i=previ; i<mult.size(); ++i ) {
			mult[i] *= prime;
		}
	}
}

#endif // deactivation of square free factorization

/** Distinct degree factorization (DDF).
 *  
 *  The implementation follows the algorithm in chapter 8 of [GCL].
 *
 *  @param[in]  a_         modular polynomial
 *  @param[out] degrees    vector containing the degrees of the factors of the
 *                         corresponding polynomials in ddfactors.
 *  @param[out] ddfactors  vector containing polynomials which factors have the
 *                         degree given in degrees.
 */
static void distinct_degree_factor(const umodpoly& a_, vector<int>& degrees, upvec& ddfactors)
{
	umodpoly a = a_;

	cl_modint_ring R = a[0].ring();
	int q = cl_I_to_int(R->modulus);
	int nhalf = degree(a)/2;

	int i = 1;
	umodpoly w(2);
	w[0] = R->zero();
	w[1] = R->one();
	umodpoly x = w;

	while ( i <= nhalf ) {
		expt_pos(w, q);
		umodpoly buf;
		rem(w, a, buf);
		w = buf;
		umodpoly wx = w - x;
		gcd(a, wx, buf);
		if ( unequal_one(buf) ) {
			degrees.push_back(i);
			ddfactors.push_back(buf);
		}
		if ( unequal_one(buf) ) {
			umodpoly buf2;
			div(a, buf, buf2);
			a = buf2;
			nhalf = degree(a)/2;
			rem(w, a, buf);
			w = buf;
		}
		++i;
	}
	if ( unequal_one(a) ) {
		degrees.push_back(degree(a));
		ddfactors.push_back(a);
	}
}

/** Modular same degree factorization.
 *  Same degree factorization is a kind of misnomer. It performs distinct degree
 *  factorization, but instead of using the Cantor-Zassenhaus algorithm it
 *  (sub-optimally) uses Berlekamp's algorithm for the factors of the same
 *  degree.
 *
 *  @param[in]  a    modular polynomial
 *  @param[out] upv  vector containing modular factors. if upv was not empty the
 *                   new elements are added at the end
 */
static void same_degree_factor(const umodpoly& a, upvec& upv)
{
	cl_modint_ring R = a[0].ring();

	vector<int> degrees;
	upvec ddfactors;
	distinct_degree_factor(a, degrees, ddfactors);

	for ( size_t i=0; i<degrees.size(); ++i ) {
		if ( degrees[i] == degree(ddfactors[i]) ) {
			upv.push_back(ddfactors[i]);
		}
		else {
			berlekamp(ddfactors[i], upv);
		}
	}
}

// Yes, we can (choose).
#define USE_SAME_DEGREE_FACTOR

/** Modular univariate factorization.
 *
 *  In principle, we have two algorithms at our disposal: Berlekamp's algorithm
 *  and same degree factorization (SDF). SDF seems to be slightly faster in
 *  almost all cases so it is activated as default.
 *
 *  @param[in]  p    modular polynomial
 *  @param[out] upv  vector containing modular factors. if upv was not empty the
 *                   new elements are added at the end
 */
static void factor_modular(const umodpoly& p, upvec& upv)
{
#ifdef USE_SAME_DEGREE_FACTOR
	same_degree_factor(p, upv);
#else
	berlekamp(p, upv);
#endif
}

/** Calculates modular polynomials s and t such that a*s+b*t==1.
 *  Assertion: a and b are relatively prime and not zero.
 *
 *  @param[in]  a  polynomial
 *  @param[in]  b  polynomial
 *  @param[out] s  polynomial
 *  @param[out] t  polynomial
 */
static void exteuclid(const umodpoly& a, const umodpoly& b, umodpoly& s, umodpoly& t)
{
	if ( degree(a) < degree(b) ) {
		exteuclid(b, a, t, s);
		return;
	}

	umodpoly one(1, a[0].ring()->one());
	umodpoly c = a; normalize_in_field(c);
	umodpoly d = b; normalize_in_field(d);
	s = one;
	t.clear();
	umodpoly d1;
	umodpoly d2 = one;
	umodpoly q;
	while ( true ) {
		div(c, d, q);
		umodpoly r = c - q * d;
		umodpoly r1 = s - q * d1;
		umodpoly r2 = t - q * d2;
		c = d;
		s = d1;
		t = d2;
		if ( r.empty() ) break;
		d = r;
		d1 = r1;
		d2 = r2;
	}
	cl_MI fac = recip(lcoeff(a) * lcoeff(c));
	umodpoly::iterator i = s.begin(), end = s.end();
	for ( ; i!=end; ++i ) {
		*i = *i * fac;
	}
	canonicalize(s);
	fac = recip(lcoeff(b) * lcoeff(c));
	i = t.begin(), end = t.end();
	for ( ; i!=end; ++i ) {
		*i = *i * fac;
	}
	canonicalize(t);
}

/** Replaces the leading coefficient in a polynomial by a given number.
 *
 *  @param[in] poly  polynomial to change
 *  @param[in] lc    new leading coefficient
 *  @return          changed polynomial
 */
static upoly replace_lc(const upoly& poly, const cl_I& lc)
{
	if ( poly.empty() ) return poly;
	upoly r = poly;
	r.back() = lc;
	return r;
}

/** Calculates the bound for the modulus.
 *  See [Mig].
 */
static inline cl_I calc_bound(const ex& a, const ex& x, int maxdeg)
{
	cl_I maxcoeff = 0;
	cl_R coeff = 0;
	for ( int i=a.degree(x); i>=a.ldegree(x); --i ) {
		cl_I aa = abs(the<cl_I>(ex_to<numeric>(a.coeff(x, i)).to_cl_N()));
		if ( aa > maxcoeff ) maxcoeff = aa;
		coeff = coeff + square(aa);
	}
	cl_I coeffnorm = ceiling1(the<cl_R>(cln::sqrt(coeff)));
	cl_I B = coeffnorm * expt_pos(cl_I(2), cl_I(maxdeg));
	return ( B > maxcoeff ) ? B : maxcoeff;
}

/** Calculates the bound for the modulus.
 *  See [Mig].
 */
static inline cl_I calc_bound(const upoly& a, int maxdeg)
{
	cl_I maxcoeff = 0;
	cl_R coeff = 0;
	for ( int i=degree(a); i>=0; --i ) {
		cl_I aa = abs(a[i]);
		if ( aa > maxcoeff ) maxcoeff = aa;
		coeff = coeff + square(aa);
	}
	cl_I coeffnorm = ceiling1(the<cl_R>(cln::sqrt(coeff)));
	cl_I B = coeffnorm * expt_pos(cl_I(2), cl_I(maxdeg));
	return ( B > maxcoeff ) ? B : maxcoeff;
}

/** Hensel lifting as used by factor_univariate().
 *
 *  The implementation follows the algorithm in chapter 6 of [GCL].
 *
 *  @param[in]  a_   primitive univariate polynomials
 *  @param[in]  p    prime number that does not divide lcoeff(a)
 *  @param[in]  u1_  modular factor of a (mod p)
 *  @param[in]  w1_  modular factor of a (mod p), relatively prime to u1_,
 *                   fulfilling  u1_*w1_ == a mod p
 *  @param[out] u    lifted factor
 *  @param[out] w    lifted factor, u*w = a
 */
static void hensel_univar(const upoly& a_, unsigned int p, const umodpoly& u1_, const umodpoly& w1_, upoly& u, upoly& w)
{
	upoly a = a_;
	const cl_modint_ring& R = u1_[0].ring();

	// calc bound B
	int maxdeg = (degree(u1_) > degree(w1_)) ? degree(u1_) : degree(w1_);
	cl_I maxmodulus = 2*calc_bound(a, maxdeg);

	// step 1
	cl_I alpha = lcoeff(a);
	a = a * alpha;
	umodpoly nu1 = u1_;
	normalize_in_field(nu1);
	umodpoly nw1 = w1_;
	normalize_in_field(nw1);
	upoly phi;
	phi = umodpoly_to_upoly(nu1) * alpha;
	umodpoly u1;
	umodpoly_from_upoly(u1, phi, R);
	phi = umodpoly_to_upoly(nw1) * alpha;
	umodpoly w1;
	umodpoly_from_upoly(w1, phi, R);

	// step 2
	umodpoly s;
	umodpoly t;
	exteuclid(u1, w1, s, t);

	// step 3
	u = replace_lc(umodpoly_to_upoly(u1), alpha);
	w = replace_lc(umodpoly_to_upoly(w1), alpha);
	upoly e = a - u * w;
	cl_I modulus = p;

	// step 4
	while ( !e.empty() && modulus < maxmodulus ) {
		upoly c = e / modulus;
		phi = umodpoly_to_upoly(s) * c;
		umodpoly sigmatilde;
		umodpoly_from_upoly(sigmatilde, phi, R);
		phi = umodpoly_to_upoly(t) * c;
		umodpoly tautilde;
		umodpoly_from_upoly(tautilde, phi, R);
		umodpoly r, q;
		remdiv(sigmatilde, w1, r, q);
		umodpoly sigma = r;
		phi = umodpoly_to_upoly(tautilde) + umodpoly_to_upoly(q) * umodpoly_to_upoly(u1);
		umodpoly tau;
		umodpoly_from_upoly(tau, phi, R);
		u = u + umodpoly_to_upoly(tau) * modulus;
		w = w + umodpoly_to_upoly(sigma) * modulus;
		e = a - u * w;
		modulus = modulus * p;
	}

	// step 5
	if ( e.empty() ) {
		cl_I g = u[0];
		for ( size_t i=1; i<u.size(); ++i ) {
			g = gcd(g, u[i]);
			if ( g == 1 ) break;
		}
		if ( g != 1 ) {
			u = u / g;
			w = w * g;
		}
		if ( alpha != 1 ) {
			w = w / alpha;
		}
	}
	else {
		u.clear();
	}
}

/** Returns a new prime number.
 *
 *  @param[in] p  prime number
 *  @return       next prime number after p
 */
static unsigned int next_prime(unsigned int p)
{
	static vector<unsigned int> primes;
	if ( primes.size() == 0 ) {
		primes.push_back(3); primes.push_back(5); primes.push_back(7);
	}
	vector<unsigned int>::const_iterator it = primes.begin();
	if ( p >= primes.back() ) {
		unsigned int candidate = primes.back() + 2;
		while ( true ) {
			size_t n = primes.size()/2;
			for ( size_t i=0; i<n; ++i ) {
				if ( candidate % primes[i] ) continue;
				candidate += 2;
				i=-1;
			}
			primes.push_back(candidate);
			if ( candidate > p ) break;
		}
		return candidate;
	}
	vector<unsigned int>::const_iterator end = primes.end();
	for ( ; it!=end; ++it ) {
		if ( *it > p ) {
			return *it;
		}
	}
	throw logic_error("next_prime: should not reach this point!");
}

/** Manages the splitting a vector of of modular factors into two partitions.
 */
class factor_partition
{
public:
	/** Takes the vector of modular factors and initializes the first partition */
	factor_partition(const upvec& factors_) : factors(factors_)
	{
		n = factors.size();
		k.resize(n, 0);
		k[0] = 1;
		cache.resize(n-1);
		one.resize(1, factors.front()[0].ring()->one());
		len = 1;
		last = 0;
		split();
	}
	int operator[](size_t i) const { return k[i]; }
	size_t size() const { return n; }
	size_t size_left() const { return n-len; }
	size_t size_right() const { return len; }
	/** Initializes the next partition.
	    Returns true, if there is one, false otherwise. */
	bool next()
	{
		if ( last == n-1 ) {
			int rem = len - 1;
			int p = last - 1;
			while ( rem ) {
				if ( k[p] ) {
					--rem;
					--p;
					continue;
				}
				last = p - 1;
				while ( k[last] == 0 ) { --last; }
				if ( last == 0 && n == 2*len ) return false;
				k[last++] = 0;
				for ( size_t i=0; i<=len-rem; ++i ) {
					k[last] = 1;
					++last;
				}
				fill(k.begin()+last, k.end(), 0);
				--last;
				split();
				return true;
			}
			last = len;
			++len;
			if ( len > n/2 ) return false;
			fill(k.begin(), k.begin()+len, 1);
			fill(k.begin()+len+1, k.end(), 0);
		}
		else {
			k[last++] = 0;
			k[last] = 1;
		}
		split();
		return true;
	}
	/** Get first partition */
	umodpoly& left() { return lr[0]; }
	/** Get second partition */
	umodpoly& right() { return lr[1]; }
private:
	void split_cached()
	{
		size_t i = 0;
		do {
			size_t pos = i;
			int group = k[i++];
			size_t d = 0;
			while ( i < n && k[i] == group ) { ++d; ++i; }
			if ( d ) {
				if ( cache[pos].size() >= d ) {
					lr[group] = lr[group] * cache[pos][d-1];
				}
				else {
					if ( cache[pos].size() == 0 ) {
						cache[pos].push_back(factors[pos] * factors[pos+1]);
					}
					size_t j = pos + cache[pos].size() + 1;
					d -= cache[pos].size();
					while ( d ) {
						umodpoly buf = cache[pos].back() * factors[j];
						cache[pos].push_back(buf);
						--d;
						++j;
					}
					lr[group] = lr[group] * cache[pos].back();
				}
			}
			else {
				lr[group] = lr[group] * factors[pos];
			}
		} while ( i < n );
	}
	void split()
	{
		lr[0] = one;
		lr[1] = one;
		if ( n > 6 ) {
			split_cached();
		}
		else {
			for ( size_t i=0; i<n; ++i ) {
				lr[k[i]] = lr[k[i]] * factors[i];
			}
		}
	}
private:
	umodpoly lr[2];
	vector< vector<umodpoly> > cache;
	upvec factors;
	umodpoly one;
	size_t n;
	size_t len;
	size_t last;
	vector<int> k;
};

/** Contains a pair of univariate polynomial and its modular factors.
 *  Used by factor_univariate().
 */
struct ModFactors
{
	upoly poly;
	upvec factors;
};

/** Univariate polynomial factorization.
 *
 *  Modular factorization is tried for several primes to minimize the number of
 *  modular factors. Then, Hensel lifting is performed.
 *
 *  @param[in]     poly   expanded square free univariate polynomial
 *  @param[in]     x      symbol
 *  @param[in,out] prime  prime number to start trying modular factorization with,
 *                        output value is the prime number actually used
 */
static ex factor_univariate(const ex& poly, const ex& x, unsigned int& prime)
{
	ex unit, cont, prim_ex;
	poly.unitcontprim(x, unit, cont, prim_ex);
	upoly prim;
	upoly_from_ex(prim, prim_ex, x);

	// determine proper prime and minimize number of modular factors
	prime = 3;
	unsigned int lastp = prime;
	cl_modint_ring R;
	unsigned int trials = 0;
	unsigned int minfactors = 0;
	cl_I lc = lcoeff(prim) * the<cl_I>(ex_to<numeric>(cont).to_cl_N());
	upvec factors;
	while ( trials < 2 ) {
		umodpoly modpoly;
		while ( true ) {
			prime = next_prime(prime);
			if ( !zerop(rem(lc, prime)) ) {
				R = find_modint_ring(prime);
				umodpoly_from_upoly(modpoly, prim, R);
				if ( squarefree(modpoly) ) break;
			}
		}

		// do modular factorization
		upvec trialfactors;
		factor_modular(modpoly, trialfactors);
		if ( trialfactors.size() <= 1 ) {
			// irreducible for sure
			return poly;
		}

		if ( minfactors == 0 || trialfactors.size() < minfactors ) {
			factors = trialfactors;
			minfactors = trialfactors.size();
			lastp = prime;
			trials = 1;
		}
		else {
			++trials;
		}
	}
	prime = lastp;
	R = find_modint_ring(prime);

	// lift all factor combinations
	stack<ModFactors> tocheck;
	ModFactors mf;
	mf.poly = prim;
	mf.factors = factors;
	tocheck.push(mf);
	upoly f1, f2;
	ex result = 1;
	while ( tocheck.size() ) {
		const size_t n = tocheck.top().factors.size();
		factor_partition part(tocheck.top().factors);
		while ( true ) {
			// call Hensel lifting
			hensel_univar(tocheck.top().poly, prime, part.left(), part.right(), f1, f2);
			if ( !f1.empty() ) {
				// successful, update the stack and the result
				if ( part.size_left() == 1 ) {
					if ( part.size_right() == 1 ) {
						result *= upoly_to_ex(f1, x) * upoly_to_ex(f2, x);
						tocheck.pop();
						break;
					}
					result *= upoly_to_ex(f1, x);
					tocheck.top().poly = f2;
					for ( size_t i=0; i<n; ++i ) {
						if ( part[i] == 0 ) {
							tocheck.top().factors.erase(tocheck.top().factors.begin()+i);
							break;
						}
					}
					break;
				}
				else if ( part.size_right() == 1 ) {
					if ( part.size_left() == 1 ) {
						result *= upoly_to_ex(f1, x) * upoly_to_ex(f2, x);
						tocheck.pop();
						break;
					}
					result *= upoly_to_ex(f2, x);
					tocheck.top().poly = f1;
					for ( size_t i=0; i<n; ++i ) {
						if ( part[i] == 1 ) {
							tocheck.top().factors.erase(tocheck.top().factors.begin()+i);
							break;
						}
					}
					break;
				}
				else {
					upvec newfactors1(part.size_left()), newfactors2(part.size_right());
					upvec::iterator i1 = newfactors1.begin(), i2 = newfactors2.begin();
					for ( size_t i=0; i<n; ++i ) {
						if ( part[i] ) {
							*i2++ = tocheck.top().factors[i];
						}
						else {
							*i1++ = tocheck.top().factors[i];
						}
					}
					tocheck.top().factors = newfactors1;
					tocheck.top().poly = f1;
					ModFactors mf;
					mf.factors = newfactors2;
					mf.poly = f2;
					tocheck.push(mf);
					break;
				}
			}
			else {
				// not successful
				if ( !part.next() ) {
					// if no more combinations left, return polynomial as
					// irreducible
					result *= upoly_to_ex(tocheck.top().poly, x);
					tocheck.pop();
					break;
				}
			}
		}
	}

	return unit * cont * result;
}

/** Second interface to factor_univariate() to be used if the information about
 *  the prime is not needed.
 */
static inline ex factor_univariate(const ex& poly, const ex& x)
{
	unsigned int prime;
	return factor_univariate(poly, x, prime);
}

/** Represents an evaluation point (<symbol>==<integer>).
 */
struct EvalPoint
{
	ex x;
	int evalpoint;
};

#ifdef DEBUGFACTOR
ostream& operator<<(ostream& o, const vector<EvalPoint>& v)
{
	for ( size_t i=0; i<v.size(); ++i ) {
		o << "(" << v[i].x << "==" << v[i].evalpoint << ") ";
	}
	return o;
}
#endif // def DEBUGFACTOR

// forward declaration
static vector<ex> multivar_diophant(const vector<ex>& a_, const ex& x, const ex& c, const vector<EvalPoint>& I, unsigned int d, unsigned int p, unsigned int k);

/** Utility function for multivariate Hensel lifting.
 *
 *  Solves the equation
 *    s_1*b_1 + ... + s_r*b_r == 1 mod p^k
 *  with deg(s_i) < deg(a_i)
 *  and with given b_1 = a_1 * ... * a_{i-1} * a_{i+1} * ... * a_r
 *
 *  The implementation follows the algorithm in chapter 6 of [GCL].
 *
 *  @param[in]  a   vector of modular univariate polynomials
 *  @param[in]  x   symbol
 *  @param[in]  p   prime number
 *  @param[in]  k   p^k is modulus
 *  @return         vector of polynomials (s_i)
 */
static upvec multiterm_eea_lift(const upvec& a, const ex& x, unsigned int p, unsigned int k)
{
	const size_t r = a.size();
	cl_modint_ring R = find_modint_ring(expt_pos(cl_I(p),k));
	upvec q(r-1);
	q[r-2] = a[r-1];
	for ( size_t j=r-2; j>=1; --j ) {
		q[j-1] = a[j] * q[j];
	}
	umodpoly beta(1, R->one());
	upvec s;
	for ( size_t j=1; j<r; ++j ) {
		vector<ex> mdarg(2);
		mdarg[0] = umodpoly_to_ex(q[j-1], x);
		mdarg[1] = umodpoly_to_ex(a[j-1], x);
		vector<EvalPoint> empty;
		vector<ex> exsigma = multivar_diophant(mdarg, x, umodpoly_to_ex(beta, x), empty, 0, p, k);
		umodpoly sigma1;
		umodpoly_from_ex(sigma1, exsigma[0], x, R);
		umodpoly sigma2;
		umodpoly_from_ex(sigma2, exsigma[1], x, R);
		beta = sigma1;
		s.push_back(sigma2);
	}
	s.push_back(beta);
	return s;
}

/** Changes the modulus of a modular polynomial. Used by eea_lift().
 *
 *  @param[in]     R  new modular ring
 *  @param[in,out] a  polynomial to change (in situ)
 */
static void change_modulus(const cl_modint_ring& R, umodpoly& a)
{
	if ( a.empty() ) return;
	cl_modint_ring oldR = a[0].ring();
	umodpoly::iterator i = a.begin(), end = a.end();
	for ( ; i!=end; ++i ) {
		*i = R->canonhom(oldR->retract(*i));
	}
	canonicalize(a);
}

/** Utility function for multivariate Hensel lifting.
 *
 *  Solves  s*a + t*b == 1 mod p^k  given a,b.
 *
 *  The implementation follows the algorithm in chapter 6 of [GCL].
 *
 *  @param[in]  a   polynomial
 *  @param[in]  b   polynomial
 *  @param[in]  x   symbol
 *  @param[in]  p   prime number
 *  @param[in]  k   p^k is modulus
 *  @param[out] s_  output polynomial
 *  @param[out] t_  output polynomial
 */
static void eea_lift(const umodpoly& a, const umodpoly& b, const ex& x, unsigned int p, unsigned int k, umodpoly& s_, umodpoly& t_)
{
	cl_modint_ring R = find_modint_ring(p);
	umodpoly amod = a;
	change_modulus(R, amod);
	umodpoly bmod = b;
	change_modulus(R, bmod);

	umodpoly smod;
	umodpoly tmod;
	exteuclid(amod, bmod, smod, tmod);

	cl_modint_ring Rpk = find_modint_ring(expt_pos(cl_I(p),k));
	umodpoly s = smod;
	change_modulus(Rpk, s);
	umodpoly t = tmod;
	change_modulus(Rpk, t);

	cl_I modulus(p);
	umodpoly one(1, Rpk->one());
	for ( size_t j=1; j<k; ++j ) {
		umodpoly e = one - a * s - b * t;
		reduce_coeff(e, modulus);
		umodpoly c = e;
		change_modulus(R, c);
		umodpoly sigmabar = smod * c;
		umodpoly taubar = tmod * c;
		umodpoly sigma, q;
		remdiv(sigmabar, bmod, sigma, q);
		umodpoly tau = taubar + q * amod;
		umodpoly sadd = sigma;
		change_modulus(Rpk, sadd);
		cl_MI modmodulus(Rpk, modulus);
		s = s + sadd * modmodulus;
		umodpoly tadd = tau;
		change_modulus(Rpk, tadd);
		t = t + tadd * modmodulus;
		modulus = modulus * p;
	}

	s_ = s; t_ = t;
}

/** Utility function for multivariate Hensel lifting.
 *
 *  Solves the equation
 *    s_1*b_1 + ... + s_r*b_r == x^m mod p^k
 *  with given b_1 = a_1 * ... * a_{i-1} * a_{i+1} * ... * a_r
 *
 *  The implementation follows the algorithm in chapter 6 of [GCL].
 *
 *  @param a  vector with univariate polynomials mod p^k
 *  @param x  symbol
 *  @param m  exponent of x^m in the equation to solve
 *  @param p  prime number
 *  @param k  p^k is modulus
 *  @return   vector of polynomials (s_i)
 */
static upvec univar_diophant(const upvec& a, const ex& x, unsigned int m, unsigned int p, unsigned int k)
{
	cl_modint_ring R = find_modint_ring(expt_pos(cl_I(p),k));

	const size_t r = a.size();
	upvec result;
	if ( r > 2 ) {
		upvec s = multiterm_eea_lift(a, x, p, k);
		for ( size_t j=0; j<r; ++j ) {
			umodpoly bmod = umodpoly_to_umodpoly(s[j], R, m);
			umodpoly buf;
			rem(bmod, a[j], buf);
			result.push_back(buf);
		}
	}
	else {
		umodpoly s, t;
		eea_lift(a[1], a[0], x, p, k, s, t);
		umodpoly bmod = umodpoly_to_umodpoly(s, R, m);
		umodpoly buf, q;
		remdiv(bmod, a[0], buf, q);
		result.push_back(buf);
		umodpoly t1mod = umodpoly_to_umodpoly(t, R, m);
		buf = t1mod + q * a[1];
		result.push_back(buf);
	}

	return result;
}

/** Map used by function make_modular().
 *  Finds every coefficient in a polynomial and replaces it by is value in the
 *  given modular ring R (symmetric representation).
 */
struct make_modular_map : public map_function {
	cl_modint_ring R;
	make_modular_map(const cl_modint_ring& R_) : R(R_) { }
	ex operator()(const ex& e)
	{
		if ( is_a<add>(e) || is_a<mul>(e) ) {
			return e.map(*this);
		}
		else if ( is_a<numeric>(e) ) {
			numeric mod(R->modulus);
			numeric halfmod = (mod-1)/2;
			cl_MI emod = R->canonhom(the<cl_I>(ex_to<numeric>(e).to_cl_N()));
			numeric n(R->retract(emod));
			if ( n > halfmod ) {
				return n-mod;
			}
			else {
				return n;
			}
		}
		return e;
	}
};

/** Helps mimicking modular multivariate polynomial arithmetic.
 *
 *  @param e  expression of which to make the coefficients equal to their value
 *            in the modular ring R (symmetric representation)
 *  @param R  modular ring
 *  @return   resulting expression
 */
static ex make_modular(const ex& e, const cl_modint_ring& R)
{
	make_modular_map map(R);
	return map(e.expand());
}

/** Utility function for multivariate Hensel lifting.
 *
 *  Returns the polynomials s_i that fulfill
 *    s_1*b_1 + ... + s_r*b_r == c mod <I^(d+1),p^k>
 *  with given b_1 = a_1 * ... * a_{i-1} * a_{i+1} * ... * a_r
 *
 *  The implementation follows the algorithm in chapter 6 of [GCL].
 *
 *  @param a_  vector of multivariate factors mod p^k
 *  @param x   symbol (equiv. x_1 in [GCL])
 *  @param c   polynomial mod p^k
 *  @param I   vector of evaluation points
 *  @param d   maximum total degree of result
 *  @param p   prime number
 *  @param k   p^k is modulus
 *  @return    vector of polynomials (s_i)
 */
static vector<ex> multivar_diophant(const vector<ex>& a_, const ex& x, const ex& c, const vector<EvalPoint>& I,
                                    unsigned int d, unsigned int p, unsigned int k)
{
	vector<ex> a = a_;

	const cl_modint_ring R = find_modint_ring(expt_pos(cl_I(p),k));
	const size_t r = a.size();
	const size_t nu = I.size() + 1;

	vector<ex> sigma;
	if ( nu > 1 ) {
		ex xnu = I.back().x;
		int alphanu = I.back().evalpoint;

		ex A = 1;
		for ( size_t i=0; i<r; ++i ) {
			A *= a[i];
		}
		vector<ex> b(r);
		for ( size_t i=0; i<r; ++i ) {
			b[i] = normal(A / a[i]);
		}

		vector<ex> anew = a;
		for ( size_t i=0; i<r; ++i ) {
			anew[i] = anew[i].subs(xnu == alphanu);
		}
		ex cnew = c.subs(xnu == alphanu);
		vector<EvalPoint> Inew = I;
		Inew.pop_back();
		sigma = multivar_diophant(anew, x, cnew, Inew, d, p, k);

		ex buf = c;
		for ( size_t i=0; i<r; ++i ) {
			buf -= sigma[i] * b[i];
		}
		ex e = make_modular(buf, R);

		ex monomial = 1;
		for ( size_t m=1; !e.is_zero() && e.has(xnu) && m<=d; ++m ) {
			monomial *= (xnu - alphanu);
			monomial = expand(monomial);
			ex cm = e.diff(ex_to<symbol>(xnu), m).subs(xnu==alphanu) / factorial(m);
			cm = make_modular(cm, R);
			if ( !cm.is_zero() ) {
				vector<ex> delta_s = multivar_diophant(anew, x, cm, Inew, d, p, k);
				ex buf = e;
				for ( size_t j=0; j<delta_s.size(); ++j ) {
					delta_s[j] *= monomial;
					sigma[j] += delta_s[j];
					buf -= delta_s[j] * b[j];
				}
				e = make_modular(buf, R);
			}
		}
	}
	else {
		upvec amod;
		for ( size_t i=0; i<a.size(); ++i ) {
			umodpoly up;
			umodpoly_from_ex(up, a[i], x, R);
			amod.push_back(up);
		}

		sigma.insert(sigma.begin(), r, 0);
		size_t nterms;
		ex z;
		if ( is_a<add>(c) ) {
			nterms = c.nops();
			z = c.op(0);
		}
		else {
			nterms = 1;
			z = c;
		}
		for ( size_t i=0; i<nterms; ++i ) {
			int m = z.degree(x);
			cl_I cm = the<cl_I>(ex_to<numeric>(z.lcoeff(x)).to_cl_N());
			upvec delta_s = univar_diophant(amod, x, m, p, k);
			cl_MI modcm;
			cl_I poscm = cm;
			while ( poscm < 0 ) {
				poscm = poscm + expt_pos(cl_I(p),k);
			}
			modcm = cl_MI(R, poscm);
			for ( size_t j=0; j<delta_s.size(); ++j ) {
				delta_s[j] = delta_s[j] * modcm;
				sigma[j] = sigma[j] + umodpoly_to_ex(delta_s[j], x);
			}
			if ( nterms > 1 ) {
				z = c.op(i+1);
			}
		}
	}

	for ( size_t i=0; i<sigma.size(); ++i ) {
		sigma[i] = make_modular(sigma[i], R);
	}

	return sigma;
}

/** Multivariate Hensel lifting.
 *  The implementation follows the algorithm in chapter 6 of [GCL].
 *  Since we don't have a data type for modular multivariate polynomials, the
 *  respective operations are done in a GiNaC::ex and the function
 *  make_modular() is then called to make the coefficient modular p^l.
 *
 *  @param a    multivariate polynomial primitive in x
 *  @param x    symbol (equiv. x_1 in [GCL])
 *  @param I    vector of evaluation points (x_2==a_2,x_3==a_3,...)
 *  @param p    prime number (should not divide lcoeff(a mod I))
 *  @param l    p^l is the modulus of the lifted univariate field
 *  @param u    vector of modular (mod p^l) factors of a mod I
 *  @param lcU  correct leading coefficient of the univariate factors of a mod I
 *  @return     list GiNaC::lst with lifted factors (multivariate factors of a),
 *              empty if Hensel lifting did not succeed
 */
static ex hensel_multivar(const ex& a, const ex& x, const vector<EvalPoint>& I,
                          unsigned int p, const cl_I& l, const upvec& u, const vector<ex>& lcU)
{
	const size_t nu = I.size() + 1;
	const cl_modint_ring R = find_modint_ring(expt_pos(cl_I(p),l));

	vector<ex> A(nu);
	A[nu-1] = a;

	for ( size_t j=nu; j>=2; --j ) {
		ex x = I[j-2].x;
		int alpha = I[j-2].evalpoint;
		A[j-2] = A[j-1].subs(x==alpha);
		A[j-2] = make_modular(A[j-2], R);
	}

	int maxdeg = a.degree(I.front().x);
	for ( size_t i=1; i<I.size(); ++i ) {
		int maxdeg2 = a.degree(I[i].x);
		if ( maxdeg2 > maxdeg ) maxdeg = maxdeg2;
	}

	const size_t n = u.size();
	vector<ex> U(n);
	for ( size_t i=0; i<n; ++i ) {
		U[i] = umodpoly_to_ex(u[i], x);
	}

	for ( size_t j=2; j<=nu; ++j ) {
		vector<ex> U1 = U;
		ex monomial = 1;
		for ( size_t m=0; m<n; ++m) {
			if ( lcU[m] != 1 ) {
				ex coef = lcU[m];
				for ( size_t i=j-1; i<nu-1; ++i ) {
					coef = coef.subs(I[i].x == I[i].evalpoint);
				}
				coef = make_modular(coef, R);
				int deg = U[m].degree(x);
				U[m] = U[m] - U[m].lcoeff(x) * pow(x,deg) + coef * pow(x,deg);
			}
		}
		ex Uprod = 1;
		for ( size_t i=0; i<n; ++i ) {
			Uprod *= U[i];
		}
		ex e = expand(A[j-1] - Uprod);

		vector<EvalPoint> newI;
		for ( size_t i=1; i<=j-2; ++i ) {
			newI.push_back(I[i-1]);
		}

		ex xj = I[j-2].x;
		int alphaj = I[j-2].evalpoint;
		size_t deg = A[j-1].degree(xj);
		for ( size_t k=1; k<=deg; ++k ) {
			if ( !e.is_zero() ) {
				monomial *= (xj - alphaj);
				monomial = expand(monomial);
				ex dif = e.diff(ex_to<symbol>(xj), k);
				ex c = dif.subs(xj==alphaj) / factorial(k);
				if ( !c.is_zero() ) {
					vector<ex> deltaU = multivar_diophant(U1, x, c, newI, maxdeg, p, cl_I_to_uint(l));
					for ( size_t i=0; i<n; ++i ) {
						deltaU[i] *= monomial;
						U[i] += deltaU[i];
						U[i] = make_modular(U[i], R);
					}
					ex Uprod = 1;
					for ( size_t i=0; i<n; ++i ) {
						Uprod *= U[i];
					}
					e = A[j-1] - Uprod;
					e = make_modular(e, R);
				}
			}
		}
	}

	ex acand = 1;
	for ( size_t i=0; i<U.size(); ++i ) {
		acand *= U[i];
	}
	if ( expand(a-acand).is_zero() ) {
		lst res;
		for ( size_t i=0; i<U.size(); ++i ) {
			res.append(U[i]);
		}
		return res;
	}
	else {
		lst res;
		return lst();
	}
}

/** Takes a factorized expression and puts the factors in a lst. The exponents
 *  of the factors are discarded, e.g. 7*x^2*(y+1)^4 --> {7,x,y+1}. The first
 *  element of the list is always the numeric coefficient.
 */
static ex put_factors_into_lst(const ex& e)
{
	lst result;
	if ( is_a<numeric>(e) ) {
		result.append(e);
		return result;
	}
	if ( is_a<power>(e) ) {
		result.append(1);
		result.append(e.op(0));
		return result;
	}
	if ( is_a<symbol>(e) || is_a<add>(e) ) {
		result.append(1);
		result.append(e);
		return result;
	}
	if ( is_a<mul>(e) ) {
		ex nfac = 1;
		for ( size_t i=0; i<e.nops(); ++i ) {
			ex op = e.op(i);
			if ( is_a<numeric>(op) ) {
				nfac = op;
			}
			if ( is_a<power>(op) ) {
				result.append(op.op(0));
			}
			if ( is_a<symbol>(op) || is_a<add>(op) ) {
				result.append(op);
			}
		}
		result.prepend(nfac);
		return result;
	}
	throw runtime_error("put_factors_into_lst: bad term.");
}

/** Checks a set of numbers for whether each number has a unique prime factor.
 *
 *  @param[in]  f  list of numbers to check
 *  @return        true: if number set is bad, false: if set is okay (has unique
 *                 prime factors)
 */
static bool checkdivisors(const lst& f)
{
	const int k = f.nops();
	numeric q, r;
	vector<numeric> d(k);
	d[0] = ex_to<numeric>(abs(f.op(0)));
	for ( int i=1; i<k; ++i ) {
		q = ex_to<numeric>(abs(f.op(i)));
		for ( int j=i-1; j>=0; --j ) {
			r = d[j];
			do {
				r = gcd(r, q);
				q = q/r;
			} while ( r != 1 );
			if ( q == 1 ) {
				return true;
			}
		}
		d[i] = q;
	}
	return false;
}

/** Generates a set of evaluation points for a multivariate polynomial.
 *  The set fulfills the following conditions:
 *  1. lcoeff(evaluated_polynomial) does not vanish
 *  2. factors of lcoeff(evaluated_polynomial) have each a unique prime factor
 *  3. evaluated_polynomial is square free
 *  See [Wan] for more details.
 *
 *  @param[in]     u        multivariate polynomial to be factored
 *  @param[in]     vn       leading coefficient of u in x (x==first symbol in syms)
 *  @param[in]     syms     set of symbols that appear in u
 *  @param[in]     f        lst containing the factors of the leading coefficient vn
 *  @param[in,out] modulus  integer modulus for random number generation (i.e. |a_i| < modulus)
 *  @param[out]    u0       returns the evaluated (univariate) polynomial
 *  @param[out]    a        returns the valid evaluation points. must have initial size equal
 *                          number of symbols-1 before calling generate_set
 */
static void generate_set(const ex& u, const ex& vn, const exset& syms, const lst& f,
                         numeric& modulus, ex& u0, vector<numeric>& a)
{
	const ex& x = *syms.begin();
	while ( true ) {
		++modulus;
		// generate a set of integers ...
		u0 = u;
		ex vna = vn;
		ex vnatry;
		exset::const_iterator s = syms.begin();
		++s;
		for ( size_t i=0; i<a.size(); ++i ) {
			do {
				a[i] = mod(numeric(rand()), 2*modulus) - modulus;
				vnatry = vna.subs(*s == a[i]);
				// ... for which the leading coefficient doesn't vanish ...
			} while ( vnatry == 0 );
			vna = vnatry;
			u0 = u0.subs(*s == a[i]);
			++s;
		}
		// ... for which u0 is square free ...
		ex g = gcd(u0, u0.diff(ex_to<symbol>(x)));
		if ( !is_a<numeric>(g) ) {
			continue;
		}
		if ( !is_a<numeric>(vn) ) {
			// ... and for which the evaluated factors have each an unique prime factor
			lst fnum = f;
			fnum.let_op(0) = fnum.op(0) * u0.content(x);
			for ( size_t i=1; i<fnum.nops(); ++i ) {
				if ( !is_a<numeric>(fnum.op(i)) ) {
					s = syms.begin();
					++s;
					for ( size_t j=0; j<a.size(); ++j, ++s ) {
						fnum.let_op(i) = fnum.op(i).subs(*s == a[j]);
					}
				}
			}
			if ( checkdivisors(fnum) ) {
				continue;
			}
		}
		// ok, we have a valid set now
		return;
	}
}

// forward declaration
static ex factor_sqrfree(const ex& poly);

/** Multivariate factorization.
 *  
 *  The implementation is based on the algorithm described in [Wan].
 *  An evaluation homomorphism (a set of integers) is determined that fulfills
 *  certain criteria. The evaluated polynomial is univariate and is factorized
 *  by factor_univariate(). The main work then is to find the correct leading
 *  coefficients of the univariate factors. They have to correspond to the
 *  factors of the (multivariate) leading coefficient of the input polynomial
 *  (as defined for a specific variable x). After that the Hensel lifting can be
 *  performed.
 *
 *  @param[in] poly  expanded, square free polynomial
 *  @param[in] syms  contains the symbols in the polynomial
 *  @return          factorized polynomial
 */
static ex factor_multivariate(const ex& poly, const exset& syms)
{
	exset::const_iterator s;
	const ex& x = *syms.begin();

	// make polynomial primitive
	ex unit, cont, pp;
	poly.unitcontprim(x, unit, cont, pp);
	if ( !is_a<numeric>(cont) ) {
		return factor_sqrfree(cont) * factor_sqrfree(pp);
	}

	// factor leading coefficient
	ex vn = pp.collect(x).lcoeff(x);
	ex vnlst;
	if ( is_a<numeric>(vn) ) {
		vnlst = lst(vn);
	}
	else {
		ex vnfactors = factor(vn);
		vnlst = put_factors_into_lst(vnfactors);
	}

	const unsigned int maxtrials = 3;
	numeric modulus = (vnlst.nops() > 3) ? vnlst.nops() : 3;
	vector<numeric> a(syms.size()-1, 0);

	// try now to factorize until we are successful
	while ( true ) {

		unsigned int trialcount = 0;
		unsigned int prime;
		int factor_count = 0;
		int min_factor_count = -1;
		ex u, delta;
		ex ufac, ufaclst;

		// try several evaluation points to reduce the number of factors
		while ( trialcount < maxtrials ) {

			// generate a set of valid evaluation points
			generate_set(pp, vn, syms, ex_to<lst>(vnlst), modulus, u, a);

			ufac = factor_univariate(u, x, prime);
			ufaclst = put_factors_into_lst(ufac);
			factor_count = ufaclst.nops()-1;
			delta = ufaclst.op(0);

			if ( factor_count <= 1 ) {
				// irreducible
				return poly;
			}
			if ( min_factor_count < 0 ) {
				// first time here
				min_factor_count = factor_count;
			}
			else if ( min_factor_count == factor_count ) {
				// one less to try
				++trialcount;
			}
			else if ( min_factor_count > factor_count ) {
				// new minimum, reset trial counter
				min_factor_count = factor_count;
				trialcount = 0;
			}
		}

		// determine true leading coefficients for the Hensel lifting
		vector<ex> C(factor_count);
		if ( is_a<numeric>(vn) ) {
			// easy case
			for ( size_t i=1; i<ufaclst.nops(); ++i ) {
				C[i-1] = ufaclst.op(i).lcoeff(x);
			}
		}
		else {
			// difficult case.
			// we use the property of the ftilde having a unique prime factor.
			// details can be found in [Wan].
			// calculate ftilde
			vector<numeric> ftilde(vnlst.nops()-1);
			for ( size_t i=0; i<ftilde.size(); ++i ) {
				ex ft = vnlst.op(i+1);
				s = syms.begin();
				++s;
				for ( size_t j=0; j<a.size(); ++j ) {
					ft = ft.subs(*s == a[j]);
					++s;
				}
				ftilde[i] = ex_to<numeric>(ft);
			}
			// calculate D and C
			vector<bool> used_flag(ftilde.size(), false);
			vector<ex> D(factor_count, 1);
			if ( delta == 1 ) {
				for ( int i=0; i<factor_count; ++i ) {
					numeric prefac = ex_to<numeric>(ufaclst.op(i+1).lcoeff(x));
					for ( int j=ftilde.size()-1; j>=0; --j ) {
						int count = 0;
						while ( irem(prefac, ftilde[j]) == 0 ) {
							prefac = iquo(prefac, ftilde[j]);
							++count;
						}
						if ( count ) {
							used_flag[j] = true;
							D[i] = D[i] * pow(vnlst.op(j+1), count);
						}
					}
					C[i] = D[i] * prefac;
				}
			}
			else {
				for ( int i=0; i<factor_count; ++i ) {
					numeric prefac = ex_to<numeric>(ufaclst.op(i+1).lcoeff(x));
					for ( int j=ftilde.size()-1; j>=0; --j ) {
						int count = 0;
						while ( irem(prefac, ftilde[j]) == 0 ) {
							prefac = iquo(prefac, ftilde[j]);
							++count;
						}
						while ( irem(ex_to<numeric>(delta)*prefac, ftilde[j]) == 0 ) {
							numeric g = gcd(prefac, ex_to<numeric>(ftilde[j]));
							prefac = iquo(prefac, g);
							delta = delta / (ftilde[j]/g);
							ufaclst.let_op(i+1) = ufaclst.op(i+1) * (ftilde[j]/g);
							++count;
						}
						if ( count ) {
							used_flag[j] = true;
							D[i] = D[i] * pow(vnlst.op(j+1), count);
						}
					}
					C[i] = D[i] * prefac;
				}
			}
			// check if something went wrong
			bool some_factor_unused = false;
			for ( size_t i=0; i<used_flag.size(); ++i ) {
				if ( !used_flag[i] ) {
					some_factor_unused = true;
					break;
				}
			}
			if ( some_factor_unused ) {
				continue;
			}
		}
		
		// multiply the remaining content of the univariate polynomial into the
		// first factor
		if ( delta != 1 ) {
			C[0] = C[0] * delta;
			ufaclst.let_op(1) = ufaclst.op(1) * delta;
		}

		// set up evaluation points
		EvalPoint ep;
		vector<EvalPoint> epv;
		s = syms.begin();
		++s;
		for ( size_t i=0; i<a.size(); ++i ) {
			ep.x = *s++;
			ep.evalpoint = a[i].to_int();
			epv.push_back(ep);
		}

		// calc bound p^l
		int maxdeg = 0;
		for ( int i=1; i<=factor_count; ++i ) {
			if ( ufaclst.op(i).degree(x) > maxdeg ) {
				maxdeg = ufaclst[i].degree(x);
			}
		}
		cl_I B = 2*calc_bound(u, x, maxdeg);
		cl_I l = 1;
		cl_I pl = prime;
		while ( pl < B ) {
			l = l + 1;
			pl = pl * prime;
		}
		
		// set up modular factors (mod p^l)
		cl_modint_ring R = find_modint_ring(expt_pos(cl_I(prime),l));
		upvec modfactors(ufaclst.nops()-1);
		for ( size_t i=1; i<ufaclst.nops(); ++i ) {
			umodpoly_from_ex(modfactors[i-1], ufaclst.op(i), x, R);
		}

		// try Hensel lifting
		ex res = hensel_multivar(pp, x, epv, prime, l, modfactors, C);
		if ( res != lst() ) {
			ex result = cont * unit;
			for ( size_t i=0; i<res.nops(); ++i ) {
				result *= res.op(i).content(x) * res.op(i).unit(x);
				result *= res.op(i).primpart(x);
			}
			return result;
		}
	}
}

/** Finds all symbols in an expression. Used by factor_sqrfree() and factor().
 */
struct find_symbols_map : public map_function {
	exset syms;
	ex operator()(const ex& e)
	{
		if ( is_a<symbol>(e) ) {
			syms.insert(e);
			return e;
		}
		return e.map(*this);
	}
};

/** Factorizes a polynomial that is square free. It calls either the univariate
 *  or the multivariate factorization functions.
 */
static ex factor_sqrfree(const ex& poly)
{
	// determine all symbols in poly
	find_symbols_map findsymbols;
	findsymbols(poly);
	if ( findsymbols.syms.size() == 0 ) {
		return poly;
	}

	if ( findsymbols.syms.size() == 1 ) {
		// univariate case
		const ex& x = *(findsymbols.syms.begin());
		if ( poly.ldegree(x) > 0 ) {
			// pull out direct factors
			int ld = poly.ldegree(x);
			ex res = factor_univariate(expand(poly/pow(x, ld)), x);
			return res * pow(x,ld);
		}
		else {
			ex res = factor_univariate(poly, x);
			return res;
		}
	}

	// multivariate case
	ex res = factor_multivariate(poly, findsymbols.syms);
	return res;
}

/** Map used by factor() when factor_options::all is given to access all
 *  subexpressions and to call factor() on them.
 */
struct apply_factor_map : public map_function {
	unsigned options;
	apply_factor_map(unsigned options_) : options(options_) { }
	ex operator()(const ex& e)
	{
		if ( e.info(info_flags::polynomial) ) {
			return factor(e, options);
		}
		if ( is_a<add>(e) ) {
			ex s1, s2;
			for ( size_t i=0; i<e.nops(); ++i ) {
				if ( e.op(i).info(info_flags::polynomial) ) {
					s1 += e.op(i);
				}
				else {
					s2 += e.op(i);
				}
			}
			s1 = s1.eval();
			s2 = s2.eval();
			return factor(s1, options) + s2.map(*this);
		}
		return e.map(*this);
	}
};

} // anonymous namespace

/** Interface function to the outside world. It checks the arguments, tries a
 *  square free factorization, and then calls factor_sqrfree to do the hard
 *  work.
 */
ex factor(const ex& poly, unsigned options)
{
	// check arguments
	if ( !poly.info(info_flags::polynomial) ) {
		if ( options & factor_options::all ) {
			options &= ~factor_options::all;
			apply_factor_map factor_map(options);
			return factor_map(poly);
		}
		return poly;
	}

	// determine all symbols in poly
	find_symbols_map findsymbols;
	findsymbols(poly);
	if ( findsymbols.syms.size() == 0 ) {
		return poly;
	}
	lst syms;
	exset::const_iterator i=findsymbols.syms.begin(), end=findsymbols.syms.end();
	for ( ; i!=end; ++i ) {
		syms.append(*i);
	}

	// make poly square free
	ex sfpoly = sqrfree(poly.expand(), syms);

	// factorize the square free components
	if ( is_a<power>(sfpoly) ) {
		// case: (polynomial)^exponent
		const ex& base = sfpoly.op(0);
		if ( !is_a<add>(base) ) {
			// simple case: (monomial)^exponent
			return sfpoly;
		}
		ex f = factor_sqrfree(base);
		return pow(f, sfpoly.op(1));
	}
	if ( is_a<mul>(sfpoly) ) {
		// case: multiple factors
		ex res = 1;
		for ( size_t i=0; i<sfpoly.nops(); ++i ) {
			const ex& t = sfpoly.op(i);
			if ( is_a<power>(t) ) {
				const ex& base = t.op(0);
				if ( !is_a<add>(base) ) {
					res *= t;
				}
				else {
					ex f = factor_sqrfree(base);
					res *= pow(f, t.op(1));
				}
			}
			else if ( is_a<add>(t) ) {
				ex f = factor_sqrfree(t);
				res *= f;
			}
			else {
				res *= t;
			}
		}
		return res;
	}
	if ( is_a<symbol>(sfpoly) ) {
		return poly;
	}
	// case: (polynomial)
	ex f = factor_sqrfree(sfpoly);
	return f;
}

} // namespace GiNaC

#ifdef DEBUGFACTOR
#include "test.h"
#endif
