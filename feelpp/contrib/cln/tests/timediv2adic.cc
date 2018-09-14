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
	if (argc < 3)
		exit(1);
	uintL a_len = atoi(argv[1]);
	uintL b_len = atoi(argv[2]);
	if (!(a_len >= b_len && b_len > 0))
		exit(1);
	CL_ALLOCA_STACK;
	uintD* a_MSDptr;
	uintD* a_LSDptr;
	uintD* b_MSDptr;
	uintD* b_LSDptr;
	uintD* q_MSDptr;
	uintD* q_LSDptr;
	uintD* q1_MSDptr;
	uintD* q1_LSDptr;
	num_stack_alloc(a_len,a_MSDptr=,a_LSDptr=);
	num_stack_alloc(b_len,b_MSDptr=,b_LSDptr=);
	num_stack_alloc(a_len,q_MSDptr=,q_LSDptr=);
	num_stack_alloc(a_len,q1_MSDptr=,q1_LSDptr=);
	random_UDS(default_random_state,a_MSDptr,a_len);
	random_UDS(default_random_state,b_MSDptr,b_len);
	lspref(b_LSDptr,0) |= 1; // force b to be odd
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { div2adic(a_len,a_LSDptr,b_len,b_LSDptr,q_LSDptr); }
	}
}
