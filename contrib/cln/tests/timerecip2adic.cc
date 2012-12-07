#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include "base/digitseq/cl_DS.h"
#include "base/digitseq/cl_2DS.h"
#include <cln/random.h>
#include "base/random/cl_random_impl.h"
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
	CL_ALLOCA_STACK;
	uintD* a_MSDptr;
	uintD* a_LSDptr;
	uintD* b_MSDptr;
	uintD* b_LSDptr;
	num_stack_alloc(len,a_MSDptr=,a_LSDptr=);
	num_stack_alloc(len,b_MSDptr=,b_LSDptr=);
	random_UDS(default_random_state,a_MSDptr,len);
	lspref(a_LSDptr,0) |= 1; // force a to be odd
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { recip2adic(len,a_LSDptr,b_LSDptr); }
	}
}
