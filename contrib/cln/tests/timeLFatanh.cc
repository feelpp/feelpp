#include <cln/number.h>
#include <cln/io.h>
#include <cln/float.h>
#include <cln/lfloat.h>
#include "float/lfloat/cl_LF.h"
#include <cln/complex.h>
#include <cln/complex_io.h>
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
	uintL len = atoi(argv[1]);
#if 0
	cl_LF one = cl_I_to_LF(1,len);
	cl_F x = scale_float(random_F(one),-1);
	cout << x << endl;
#else
	cl_F x = sqrt(cl_I_to_LF(2,len))-1;
#endif
	cl_N y;
	ln(cl_I_to_LF(1000,len+10)); // fill cache
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { y = atanh(x); }
	}
	cout << y << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { y = ln((1+x)/(1-x))/2; }
	}
	cout << y << endl;
}
