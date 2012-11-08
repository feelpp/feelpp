#include <cln/number.h>
#include <cln/io.h>
#include <cln/float.h>
#include <cln/lfloat.h>
#include "float/lfloat/cl_LF.h"
#include <cln/real.h>
#include <cln/random.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>

int main (int argc, char * argv[])
{
	int repetitions = 1;
	if ((argc >= 3) && !strcmp(argv[1],"-r")) {
		repetitions = atoi(argv[2]);
		argc -= 2; argv += 2;
	}
	if (argc < 3)
		exit(1);
	extern double cl_sin_factor;
	extern int cl_sin_algo;
	extern int cl_aux_algo;
	uintL len = atoi(argv[1]);
	cl_sin_factor = atof(argv[2]);
#if 0
	cl_LF one = cl_I_to_LF(1,len);
	cl_F x = scale_float(random_F(one),-1);
	cout << x << endl;
#else
	cl_F x = sqrt(cl_I_to_LF(2,len))-1;
#endif
	cl_F y;
	cl_sin_algo = 0;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { y = sin(x); }
	}
	cout << y << endl;
	cl_sin_algo = 1; cl_aux_algo = 0;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { y = sin(x); }
	}
	cout << y << endl;
	cl_sin_algo = 1; cl_aux_algo = 1;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { y = sin(x); }
	}
	cout << y << endl;
}
