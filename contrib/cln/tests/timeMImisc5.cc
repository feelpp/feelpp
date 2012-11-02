#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include <cln/modinteger.h>
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
	if (argc < 1)
		exit(1);
	cl_I p = "1269281897404513557783934075031171555202695168107";
	cl_modint_ring R = find_modint_ring(p);
	{
		cl_MI a = R->canonhom("1111111111111111111111111111111111111111111111111");
		cl_MI b = R->canonhom("777777777777777777777777777777777777777777777777");
		cout << "product modulo p" << endl;
		{ CL_TIMING;
		  for (int rep = repetitions; rep > 0; rep--)
		    { cl_MI c = R->mul(a,b); }
		}
		cout << "square modulo p" << endl;
		{ CL_TIMING;
		  for (int rep = repetitions; rep > 0; rep--)
		    { cl_MI c = R->square(a); }
		}
		cout << "quotient modulo p" << endl;
		{ CL_TIMING;
		  for (int rep = repetitions; rep > 0; rep--)
		    { cl_MI c = R->div(a,b); }
		}
	}
	{
		cl_MI a = R->canonhom("1234567890123456789012345678901234567890123456789");
		cl_MI b = R->canonhom("909090909090909090909090909090909090909090909090");
		cout << "product modulo p" << endl;
		{ CL_TIMING;
		  for (int rep = repetitions; rep > 0; rep--)
		    { cl_MI c = R->mul(a,b); }
		}
		cout << "square modulo p" << endl;
		{ CL_TIMING;
		  for (int rep = repetitions; rep > 0; rep--)
		    { cl_MI c = R->square(a); }
		}
		cout << "quotient modulo p" << endl;
		{ CL_TIMING;
		  for (int rep = repetitions; rep > 0; rep--)
		    { cl_MI c = R->div(a,b); }
		}
	}
}
