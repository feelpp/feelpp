#include <cln/number.h>
#include <cln/io.h>
#include <cln/float.h>
#include <cln/real.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>
#include "float/lfloat/cl_LF.h"
namespace cln
{
// FIXME: don't use internal functions.
extern cl_LF zeta (int s, uintC len);
extern cl_LF compute_zeta_exp (int s, uintC len);
extern cl_LF compute_zeta_cvz1 (int s, uintC len);
extern cl_LF compute_zeta_cvz2 (int s, uintC len);
extern cl_LF zeta3 (uintC len);
}
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
	cl_LF p;
	ln(cl_I_to_LF(1000,len+10)); // fill cache
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_zeta_exp(3,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_zeta_cvz1(3,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_zeta_cvz2(3,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = zeta3(len); }
	}
	cout << p << endl;
}
