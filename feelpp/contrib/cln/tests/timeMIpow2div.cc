#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include <cln/modinteger.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>
using namespace cln;

int main (int argc, char * argv[])
{
	int repetitions = 1;
	if ((argc >= 3) && !strcmp(argv[1],"-r")) {
		repetitions = atoi(argv[2]);
		argc -= 2; argv += 2;
	}
	if (argc < 2)
		exit(1);
	uintL len = atoi(argv[1]);
	cl_modint_ring R = find_modint_ring((cl_I)1 << (intDsize*len));
	cl_MI a = R->random();
	cl_MI b;
	do { b = R->random(); } while (!oddp(R->retract(b)));
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { cl_MI c = R->div(a,b); }
	}
}
