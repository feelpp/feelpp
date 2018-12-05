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
#include <iostream>
using namespace cln;
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
	cl_LF one = cl_I_to_LF(1,len);
	cl_F x = random_F(one);
	cl_F y;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { y = sqrt(x); }
	}
	cout << y << endl;
}
