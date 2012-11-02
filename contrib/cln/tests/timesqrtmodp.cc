#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include <cln/integer_io.h>
#include <cln/modinteger.h>
#include <cln/numtheory.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>
using namespace cln;
#include <iostream>
using namespace std;

int main (int argc, char * argv[])
{
	int repetitions = 1;
	if ((argc >= 3) && !strcmp(argv[1],"-r")) {
		repetitions = atoi(argv[2]);
		argc -= 2; argv += 2;
	}
	if (argc < 2)
		exit(1);
	cl_I len = cl_I(argv[1]);
	int e = (argc > 2 ? atoi(argv[2]) : 0);
	if (e < 1)
		e = 1;
	if (len <= e)
		exit(0);
	cl_I p;
	do {
		p = ((random_I((cl_I)1 << (len-1-e))*2+1) << e) + 1;
	} while (!isprobprime(p));
	cout << "p = " << p << endl;
	cl_modint_ring R = find_modint_ring(p);
	cl_MI x = R->random();
	cl_MI a = square(x);
	sqrt_mod_p_t sol;
#if 0
	extern int cl_sqrt_algo;
	cl_sqrt_algo = 1;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { sol = sqrt_mod_p(R,a); }
	}
	if (sol.condition)
		cerr << "p not prime!" << endl;
	else {
		if (sol.solutions == 0)
			cerr << "No sqrt found!" << endl;
		if (!(sol.solution[0] == x || sol.solution[0] == -x))
			cerr << "Wrong result!" << endl;
	}
	cl_sqrt_algo = 2;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { sol = sqrt_mod_p(R,a); }
	}
	if (sol.condition)
		cerr << "p not prime!" << endl;
	else {
		if (sol.solutions == 0)
			cerr << "No sqrt found!" << endl;
		if (!(sol.solution[0] == x || sol.solution[0] == -x))
			cerr << "Wrong result!" << endl;
	}
	cl_sqrt_algo = 3;
#endif
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { sol = sqrt_mod_p(R,a); }
	}
	if (sol.condition)
		cerr << "p not prime!" << endl;
	else {
		if (sol.solutions == 0)
			cerr << "No sqrt found!" << endl;
		if (!(sol.solution[0] == x || sol.solution[0] == -x))
			cerr << "Wrong result!" << endl;
	}
}
