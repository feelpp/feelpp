// Check whether a mersenne number is prime,
// using the Lucas-Lehmer test.
// [Donald Ervin Knuth: The Art of Computer Programming, Vol. II:
//  Seminumerical Algorithms, second edition. Section 4.5.4, p. 391.]

// We work with integers.
#include <cln/integer.h>

using namespace std;
using namespace cln;

// Checks whether 2^q-1 is prime, q an odd prime.
bool mersenne_prime_p (int q)
{
	cl_I m = ((cl_I)1 << q) - 1;
	int i;
	cl_I L_i;
	for (i = 0, L_i = 4; i < q-2; i++)
		L_i = mod(L_i*L_i - 2, m);
	return (L_i==0);
}

// Same thing, but optimized.
bool mersenne_prime_p_opt (int q)
{
	cl_I m = ((cl_I)1 << q) - 1;
	int i;
	cl_I L_i;
	for (i = 0, L_i = 4; i < q-2; i++) {
		L_i = square(L_i) - 2;
		L_i = ldb(L_i,cl_byte(q,q)) + ldb(L_i,cl_byte(q,0));
		if (L_i >= m)
			L_i = L_i - m;
	}
	return (L_i==0);
}

// Now we work with modular integers.
#include <cln/modinteger.h>

// Same thing, but using modular integers.
bool mersenne_prime_p_modint (int q)
{
	cl_I m = ((cl_I)1 << q) - 1;
	cl_modint_ring R = find_modint_ring(m); // Z/mZ
	int i;
	cl_MI L_i;
	for (i = 0, L_i = R->canonhom(4); i < q-2; i++)
		L_i = R->minus(R->square(L_i),R->canonhom(2));
	return R->equal(L_i,R->zero());
}

#include <cln/io.h> // we do I/O
#include <cstdlib>  // declares exit()
#include <cln/timing.h>

int main (int argc, char* argv[])
{
	if (!(argc == 2)) {
		cerr << "Usage: lucaslehmer exponent" << endl;
		exit(1);
	}
	int q = atoi(argv[1]);
	if (!(q >= 2 && ((q % 2)==1))) {
		cerr << "Usage: lucaslehmer q  with q odd prime" << endl;
		exit(1);
	}
	bool isprime;
	{ CL_TIMING; isprime = mersenne_prime_p(q); }
	{ CL_TIMING; isprime = mersenne_prime_p_opt(q); }
	{ CL_TIMING; isprime = mersenne_prime_p_modint(q); }
	cout << "2^" << q << "-1 is ";
	if (isprime)
		cout << "prime" << endl;
	else
		cout << "composite" << endl;
}

// Computing time on a i486, 33 MHz:
//  1279: 2.02 s
//  2281: 8.74 s
// 44497: 14957 s
