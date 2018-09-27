// Benchmarks from the LiDIA home page

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
#include <iostream>
using namespace cln;
using namespace std;

// Timings on Linux i486 33 MHz, 1000 decimal places = 104 32-bit words.
// Function              LiDIA       Pari       CLISP         CLN
//  pi                               0.17 / 0   0.38 / 0      0.12 / 0
//  gamma                            7.51 / 0    --           3.75 / 0
//  e                                1.20       5.06 / 0.20   0.66 / 0.10
//  multiplication                   0.018      0.010         0.010
//  sqrt(3)                          0.051      0.012         0.01
//  exp(log(2))                      3.13       0.08 / 3.94   0.04
//  log(exp(2))                      3.07       4.93 / 2.75   1.34
//  sin(pi/3)                        1.53       2.98          0.58
//  cos(pi/3)                        1.59       2.16          0.58
//  arcsin(sqrt(3)/2)                4.24       2.22          1.26
//  arccos(sqrt(3)/2)                4.26       2.22          1.26
//  sinh(log(2))                     3.16       2.02          0.03
//  cosh(log(2))                     3.17       2.09          0.04
//  arsinh(pi)                       1.93       2.62          0.65
//  arcosh(pi)                       1.95       2.26          0.69
// (Versions: Pari 1.39, clisp-1996-07-22, cln-1996-11-17.)

// Timings on Solaris Sparc 20, 75 MHz, 1000 decimal places = 104 32-bit words.
// Function              LiDIA       Pari       CLISP         CLN
//  pi                   0.06 / 0    0.04 / 0                 0.028 / 0
//  gamma                1.98 / 0    1.26 / 0                 1.99 / 0
//  e                    0           0.15                     0.15
//  multiplication
//  sqrt(3)              0.0025      0.01                     0.0025
//  exp(log(2))          0.25        0.39                     0.010
//  log(exp(2))          0.21        0.39                     0.31
//  sin(pi/3)            0.071       0.18                     0.13
//  cos(pi/3)            0.070       0.19                     0.13
//  arcsin(sqrt(3)/2)    0.30        0.55                     0.27
//  arccos(sqrt(3)/2)    0.30        0.55                     0.27
//  sinh(log(2))         0.25        0.41                     0.010
//  cosh(log(2))         0.25        0.40                     0.010
//  arsinh(pi)           0.16        0.26                     0.144
//  arcosh(pi)           0.16        0.26                     0.153
// (Versions: Pari 1.39, LiDIA 1.2.1 with libI, cln-1996-10-13.)

int main (int argc, char * argv[])
{
	int repetitions = 1;
	if ((argc >= 3) && !strcmp(argv[1],"-r")) {
		repetitions = atoi(argv[2]);
		argc -= 2; argv += 2;
	}
	if (argc < 1)
		exit(1);

	fprint(std::cerr, "Number of repetitions: ");
	fprintdecimal(std::cerr, repetitions);
	fprint(std::cerr, "\n");

	float_format_t prec = float_format(1000);

	fprint(std::cerr, "pi\n");
	{ cl_F p;
	  { CL_TIMING; p = pi(prec); }
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_F p = pi(prec); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "gamma\n");
	{ cl_F p;
	  { CL_TIMING; p = eulerconst(prec); }
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_F p = eulerconst(prec); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "e\n");
	{ cl_F p = exp1(prec);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_F p = exp1(prec); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "sqrt(3)\n");
	{ cl_R p = sqrt(cl_float(3,prec));
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_R p = sqrt(cl_float(3,prec)); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "exp(log(2))\n");
	{ cl_N p = exp(log(cl_float(2,prec)));
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_N p = exp(log(cl_float(2,prec))); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "log(exp(2))\n");
	{ cl_N p = log(exp(cl_float(2,prec)));
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_N p = log(exp(cl_float(2,prec))); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "sin(pi/3)\n");
	{ cl_R p = sin(pi(prec)/3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_R p = sin(pi(prec)/3); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "cos(pi/3)\n");
	{ cl_R p = cos(pi(prec)/3);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_R p = cos(pi(prec)/3); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "arcsin(sqrt(3)/2)\n");
	{ cl_N p = asin(sqrt(cl_float(3,prec))/2);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_N p = asin(sqrt(cl_float(3,prec))/2); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "arccos(sqrt(3)/2)\n");
	{ cl_N p = acos(sqrt(cl_float(3,prec))/2);
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_N p = acos(sqrt(cl_float(3,prec))/2); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "sinh(log(2))\n");
	{ cl_N p = sinh(log(cl_float(2,prec)));
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_N p = sinh(log(cl_float(2,prec))); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "cosh(log(2))\n");
	{ cl_N p = cosh(log(cl_float(2,prec)));
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_N p = cosh(log(cl_float(2,prec))); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "arsinh(pi)\n");
	{ cl_N p = asinh(pi(prec));
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_N p = asinh(pi(prec)); }
	  }
	  cout << p << endl << endl;
	}

	fprint(std::cerr, "arcosh(pi)\n");
	{ cl_N p = acosh(pi(prec));
	  { CL_TIMING;
	    for (int rep = repetitions; rep > 0; rep--)
	      { cl_N p = acosh(pi(prec)); }
	  }
	  cout << p << endl << endl;
	}

}
