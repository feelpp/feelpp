#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include <cln/float.h>
#include <cln/float_io.h>
#include <cln/real.h>
#include <cln/real_io.h>
#include <cln/complex.h>
#include <cln/complex_io.h>
#include <cstdlib>
#include <cstring>
#include <cln/timing.h>

using namespace std;
using namespace cln;

int main (int argc, char * argv[])
{
	int digits = 100;
	int repetitions = 1;
	while (argc >= 3) {
		if (!strcmp(argv[1],"-r")) {
			repetitions = atoi(argv[2]);
			argc -= 2; argv += 2;
			continue;
		}
		if (!strcmp(argv[1],"-n")) {
			digits = atoi(argv[2]);
			argc -= 2; argv += 2;
			continue;
		}
		break;
	}
	if (argc < 1)
		exit(1);

	cerr << "Number of digits: " << digits << endl;
	cerr << "Number of repetitions (except for pi,euler,e): " << repetitions << endl;

	float_format_t prec = float_format(digits);
	cl_F x1 = sqrt(cl_float(2,prec));
	cl_F x2 = sqrt(cl_float(3,prec));
	cl_F x3 = The(cl_F)(log(cl_float(2,prec)));

	cerr << "multiplication" << endl;
	{ cl_F r = x1*x2;
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = x1*x2; }
	  }
	  cout << r << endl << endl;
	}

	cerr << "sqrt" << endl;
	{ cl_F r = sqrt(x3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = sqrt(x3); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "pi" << endl;
	{ cl_F r;
	  { CL_TIMING; r = pi(prec); }
	  cout << r << endl << endl;
	}

	cerr << "eulerconst" << endl;
	{ cl_F r;
	  { CL_TIMING; r = eulerconst(prec); }
	  cout << r << endl << endl;
	}

	cerr << "e" << endl;
	{ cl_F r;
	  { CL_TIMING; r = exp1(prec); }
	  cout << r << endl << endl;
	}

	cerr << "exp" << endl;
	{ cl_F r = exp(-x1);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = exp(-x1); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "log" << endl;
	{ cl_N r = log(x2);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = log(x2); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "sin" << endl;
	{ cl_R r = sin(5*x1);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = sin(5*x1); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "cos" << endl;
	{ cl_R r = cos(5*x1);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = cos(5*x1); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "asin" << endl;
	{ cl_N r = asin(x3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = asin(x3); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "acos" << endl;
	{ cl_N r = acos(x3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = acos(x3); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "atan" << endl;
	{ cl_F r = atan(x3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = atan(x3); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "sinh" << endl;
	{ cl_F r = sinh(x2);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = sinh(x2); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "cosh" << endl;
	{ cl_F r = cosh(x2);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = cosh(x2); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "asinh" << endl;
	{ cl_N r = asinh(x3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = asinh(x3); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "acosh" << endl;
	{ cl_N r = acosh(1+x3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = acosh(1+x3); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "atanh" << endl;
	{ cl_N r = atanh(x3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = atanh(x3); }
	  }
	  cout << r << endl << endl;
	}

}
