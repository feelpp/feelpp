#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include <cln/integer_io.h>
#include <cln/float.h>
#include <cln/real.h>
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
	cerr << "Number of repetitions: " << repetitions << endl;

	float_format_t prec = float_format(digits);
	float_format_t prec2 = float_format(digits*2);
	cl_I pow = expt_pos(10,digits);
	cl_I x1 = floor1((sqrt(cl_float(5,prec2))+1)/2 * expt_pos(pow,2));
	cl_I x2 = floor1(sqrt(cl_float(3,prec)) * pow);
	cl_I x3 = pow+1;

	cerr << "multiplication" << endl;
	{ cl_I r = x1*x2;
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = x1*x2; }
	  }
	  cout << r << endl << endl;
	}

	cerr << "division" << endl;
	{ cl_I_div_t qr = floor2(x1,x2);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { qr = floor2(x1,x2); }
	  }
	  cout << qr.quotient << endl << qr.remainder << endl << endl;
	}

	cerr << "isqrt" << endl;
	{ cl_I r = isqrt(x3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = isqrt(x3); }
	  }
	  cout << r << endl << endl;
	}

	cerr << "gcd" << endl;
	{ cl_I r = gcd(x1,x2);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = gcd(x1,x2); }
	  }
	  cout << r << endl << endl;
	}

}
