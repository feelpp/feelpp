#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include <cln/modinteger.h>
#include <cln/univpoly.h>
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

	int n = atoi(argv[1]);
	cl_I m = 100001;
	int i;

	cl_modint_ring R1 = find_modint_ring(m);
	cl_univpoly_ring PR1 = find_univpoly_ring(R1);
	cl_UP p1 = PR1->create(n-1);
	for (i = 0; i < n; i++)
		p1.set_coeff(i, R1->canonhom((int)(1.618033989*i*i)));
	p1.finalize();

	cout << p1 << endl;

	cl_UP sp1 = PR1->zero();
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { sp1 = square(p1); }
	}

	cout << sp1 << endl;
}

// Time:
//    n      modint    modint2   neu
//    2      0.000123  0.000082  0.000086
//    5      0.00051   0.00031   0.00032
//   10      0.00169   0.00095   0.00100
//   25      0.0089    0.0049    0.0053
//   50      0.031     0.018     0.020
//  100      0.118     0.070     0.079
//  250      0.72      0.43      0.48
//  500      2.87      1.76      1.91
// 1000     11.4       7.0       8.0
