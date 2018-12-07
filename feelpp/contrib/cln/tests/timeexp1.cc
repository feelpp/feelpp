#include <cln/number.h>
#include <cln/io.h>
#include <cln/float.h>
#include <cln/real.h>
#include <cln/random.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>
namespace cln
{
// FIXME: don't use internal functions!
extern cl_LF compute_exp1 (uintC len);
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
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = compute_exp1(len); }
	}
//	cout << p << endl;
}
