#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include <cln/random.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>
#include <cl_print.h>
#include <cln/malloc.h>

int main (int argc, char * argv[])
{
	int repetitions = 1;
	if ((argc >= 3) && !strcmp(argv[1],"-r")) {
		repetitions = atoi(argv[2]);
		argc -= 2; argv += 2;
	}
	if (argc < 2)
		exit(1);
	cl_I m = cl_I(argv[1]);
	cl_I M = (cl_I)1 << (intDsize*m);
	cl_I a = random_I(M);
	extern int cl_digits_algo;
	// One run to fill the cache.
	{
		char* p = (cl_digits_algo = 0, cl_decimal_string(a));
		char* q = (cl_digits_algo = 1, cl_decimal_string(a));
		if (strcmp(p,q)) abort();
		free_hook(p);
		free_hook(q);
	}
	// Now start the timing.
	cl_digits_algo = 0;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { char* p = cl_decimal_string(a); free_hook(p); }
	}
	cl_digits_algo = 1;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { char* p = cl_decimal_string(a); free_hook(p); }
	}
}
