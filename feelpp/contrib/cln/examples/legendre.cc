// Compute the Legendre polynomials.

#include <cln/number.h>
#include <cln/integer.h>
#include <cln/rational.h>
#include <cln/univpoly.h>
#include <cln/modinteger.h>
#include <cln/univpoly_rational.h>
#include <cln/univpoly_modint.h>
#include <cln/io.h>
#include <cstdlib>

using namespace std;
using namespace cln;

// Computes the n-th Legendre polynomial in R[x], using the formula
// P_n(x) = 1/(2^n n!) * (d/dx)^n (x^2-1)^n. (Assume n >= 0.)

const cl_UP_RA legendre (const cl_rational_ring& R, int n)
{
	cl_univpoly_rational_ring PR = find_univpoly_ring(R);
	cl_UP_RA b = PR->create(2);
	b.set_coeff(2,1);
	b.set_coeff(1,0);
	b.set_coeff(0,-1);
	b.finalize(); // b is now x^2-1
	cl_UP_RA p = (n==0 ? PR->one() : expt_pos(b,n));
	for (int i = 0; i < n; i++)
		p = deriv(p);
	cl_RA factor = recip(factorial(n)*ash(1,n));
	for (int j = degree(p); j >= 0; j--)
		p.set_coeff(j, coeff(p,j) * factor);
	p.finalize();
	return p;
}

const cl_UP_MI legendre (const cl_modint_ring& R, int n)
{
	cl_univpoly_modint_ring PR = find_univpoly_ring(R);
	cl_UP_MI b = PR->create(2);
	b.set_coeff(2,R->canonhom(1));
	b.set_coeff(1,R->canonhom(0));
	b.set_coeff(0,R->canonhom(-1));
	b.finalize(); // b is now x^2-1
	cl_UP_MI p = (n==0 ? PR->one() : expt_pos(b,n));
	for (int i = 0; i < n; i++)
		p = deriv(p);
	cl_MI factor = recip(R->canonhom(factorial(n)*ash(1,n)));
	for (int j = degree(p); j >= 0; j--)
		p.set_coeff(j, coeff(p,j) * factor);
	p.finalize();
	return p;
}

int main (int argc, char* argv[])
{
	if (!(argc == 2 || argc == 3)) {
		cerr << "Usage: legendre n [m]" << endl;
		exit(1);
	}
	int n = atoi(argv[1]);
	if (!(n >= 0)) {
		cerr << "Usage: legendre n [m]  with n >= 0" << endl;
		exit(1);
	}
	if (argc == 2) {
		cl_UP p = legendre(cl_RA_ring,n);
		cout << p << endl;
	} else {
		cl_I m = argv[2];
		cl_UP p = legendre(find_modint_ring(m),n);
		cout << p << endl;
	}
	return 0;
}
