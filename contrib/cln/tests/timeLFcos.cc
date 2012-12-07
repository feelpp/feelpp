#include <cln/number.h>
#include <cln/io.h>
#include <cln/float.h>
#include <cln/float_io.h>
#include <cln/lfloat.h>
#include "float/lfloat/cl_LF.h"
#include <cln/real.h>
#include <cln/random.h>
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
	extern int cl_cos_algo;
	uintL len = atoi(argv[1]);
#if 0
	cl_LF one = cl_I_to_LF(1,len);
	cl_F x = scale_float(random_F(one),-1);
	cout << x << endl;
#else
	cl_F x = sqrt(cl_I_to_LF(2,len))-1;
#endif
	cl_F y;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { y = cos(x); }
	}
	cout << y << endl;
}
