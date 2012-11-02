#include <cln/number.h>
#include <cln/io.h>
#include <cln/float.h>
#include <cln/real.h>
#include <cln/random.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>
#include "float/lfloat/cl_LF.h"
namespace cln
{
// FIXME: don't use internal functions.
extern cl_LF compute_eulerconst (uintC len);
extern cl_LF compute_eulerconst_expintegral (uintC len);
extern cl_LF compute_eulerconst_expintegral1 (uintC len);
extern cl_LF compute_eulerconst_expintegral2 (uintC len);
extern cl_LF compute_eulerconst_besselintegral1 (uintC len);
extern cl_LF compute_eulerconst_besselintegral2 (uintC len);
extern cl_LF compute_eulerconst_besselintegral3 (uintC len);
extern cl_LF compute_eulerconst_besselintegral4 (uintC len);
}

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


	cl_LF p;
	ln(cl_I_to_LF(1000,len+10)); // fill cache
#if 0
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_eulerconst(len); }
	}
#else
#if 0
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_eulerconst_expintegral(len); }
	}
//	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_eulerconst_expintegral1(len); }
	}
//	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_eulerconst_expintegral2(len); }
	}
//	cout << p << endl;
#endif
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_eulerconst_besselintegral1(len); }
	}
//	cout << p << endl;
#if 0
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_eulerconst_besselintegral2(len); }
	}
//	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_eulerconst_besselintegral3(len); }
	}
//	cout << p << endl;
#endif
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_eulerconst_besselintegral4(len); }
	}
//	cout << p << endl;
#endif
}
