#include <cln/number.h>
#include <cln/io.h>
#include <cln/float.h>
#include <cln/real.h>
#include <cln/random.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>
#include "float/lfloat/cl_LF.h"
#include <iostream>
namespace cln
{
// FIXME: using internal functions is a bad idea (even if it works
// on some ELF systems)
extern cl_LF compute_catalanconst_ramanujan (uintC len);
extern cl_LF compute_catalanconst_ramanujan_fast (uintC len);
extern cl_LF compute_catalanconst_expintegral1 (uintC len);
extern cl_LF compute_catalanconst_expintegral2 (uintC len);
extern cl_LF compute_catalanconst_cvz1 (uintC len);
extern cl_LF compute_catalanconst_cvz2 (uintC len);
}

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


	cl_LF p;
	ln(cl_I_to_LF(1000,len+10)); // fill cache
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_catalanconst_ramanujan(len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_catalanconst_ramanujan_fast(len); }
	}
	cout << p << endl;
#if 0
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_catalanconst_expintegral1(len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_catalanconst_expintegral2(len); }
	}
	cout << p << endl;
#endif
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_catalanconst_cvz1(len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_catalanconst_cvz2(len); }
	}
	cout << p << endl;
}
