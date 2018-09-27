//#define WANT_OBFUSCATING_OPERATORS
#include <cln/number.h>
#include <cln/io.h>
//#include <cln/complex.h>
//#include <cln/real.h>
//#include <cln/real_io.h>
//#include <cln/output.h>
//#include <cln/ffloat.h>
#include <cln/integer.h>
#include <cln/integer_io.h>
//#include <cln/modinteger.h>
//#include <cln/numtheory.h>
#include <cln/timing.h>
using namespace cln;
#include <iostream>
using namespace std;

#define DUMP(expr)  \
	fprint(cout, #expr" = "); fprint(cout, expr); fprint(cout, "\n");

int main (int argc, char* argv[])
{
	(void)argc; (void)argv;

#if 0
	cl_F archimedes = pi((float_format_t)10000);
	fprint(cout, archimedes);
	fprint(cout, "\n");
#endif

#if 0
	cl_FF a;
	cl_FF x1 = "-0.2173f0";
	cl_FF x2 = "5.5084f9";
	cl_FF y = "-1.19698f9";
	fprint(cout, "x1 = "); print_float_binary(cout,x1); fprint(cout, " = "); fprint(cout,x1); fprint(cout, "\n");
	fprint(cout, "x2 = "); print_float_binary(cout,x2); fprint(cout, " = "); fprint(cout,x2); fprint(cout, "\n");
	fprint(cout, "y = "); print_float_binary(cout,y); fprint(cout, " = "); fprint(cout,y); fprint(cout, "\n");
	cl_FF x = x1*x2;
	fprint(cout, "x1*x2 = "); print_float_binary(cout,x); fprint(cout, " = "); fprint(cout,x); fprint(cout, "\n");
#endif

#if 0
	cl_I x = 10;
	cl_I y = ++x;
	x *= 2;
	x++;
	fprint(cout, "x = "); fprint(cout, x); fprint(cout, "\n");
	fprint(cout, "y = "); fprint(cout, y); fprint(cout, "\n");
#endif

#if 0
	cl_I a = "77371252437321868671713407";
	cl_I w;
	bool squarep = isqrt(a,&w);
	DUMP(squarep);
	DUMP(w);
	DUMP(expt_pos(w,2) <= a);
	DUMP(a < expt_pos(w+1,2));
#endif

#if 0
	cl_I m = "79228162513111556826425457664";
	cl_I a = "19787815858762768436681494528";
	cl_modint_ring R = find_modint_ring(m);
	cl_I b = R->retract(R->canonhom(a));
	cl_I c = mod(a,abs(m));
	DUMP(b);
	DUMP(c);
	DUMP(b==c);
#endif

#if 0
	cl_N x = argv[1];
	cl_N y = sinh(x);
	{ CL_TIMING; y = sinh(x); }
	cout << y << endl;
#endif

#if 0
	cl_I x = argv[1];
	cout << x << " is " << (isprobprime(x) ? "" : "not ") << "prime" << endl;
#endif

#if 0
	float_format_t f = float_format(atoi(argv[1]));
	extern cl_LF zeta3 (uintC len);
	uintC len = (uintC)f/intDsize+1;
	{ CL_TIMING; cout << zeta(2,f) << endl; }
	{ CL_TIMING; cout << expt(pi(f),2)/6 << endl; }
	{ CL_TIMING; cout << zeta(3,f) << endl; }
	{ CL_TIMING; cout << zeta3(len) << endl; }
	{ CL_TIMING; cout << zeta(4,f) << endl; }
#endif

	cl_I a = cl_I(argv[1]);
	cl_I b = cl_I(argv[2]);
	cl_I u;
	cl_I v;
	cl_I g = xgcd(a,b,&u,&v);
	cout << "a = " << a << endl;
	cout << "b = " << b << endl;
	cout << "gcd = " << gcd(a,b) << endl;
	cout << "g = " << g << endl;
	cout << "u = " << u << endl;
	cout << "v = " << v << endl;

#if 0
	cl_F x = argv[1];
	cout << x << endl;
#endif

}
