// Computation of arctan(1/m) (m integer) to high precision.

#include "cln/integer.h"
#include "cln/rational.h"
#include "cln/real.h"
#include "cln/lfloat.h"
#include "cl_LF.h"
#include "cl_LF_tran.h"
#include "cl_alloca.h"
#include <cstdlib>
#include <cstring>
#include "cln/timing.h"

#undef floor
#include <cmath>
#define floor cln_floor

using namespace cln;

// Method 1: atan(1/m) = sum(n=0..infty, (-1)^n/(2n+1) * 1/m^(2n+1))
// Method 2: atan(1/m) = sum(n=0..infty, 4^n*n!^2/(2n+1)! * m/(m^2+1)^(n+1))
// a. Using long floats.                     [N^2]
// b. Simulating long floats using integers. [N^2]
// c. Using integers, no binary splitting.   [N^2]
// d. Using integers, with binary splitting. [FAST]
// Method 3: general built-in algorithm.     [FAST]


// Method 1: atan(1/m) = sum(n=0..infty, (-1)^n/(2n+1) * 1/m^(2n+1))

const cl_LF atan_recip_1a (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_LF eps = scale_float(cl_I_to_LF(1,actuallen),-intDsize*(sintC)actuallen);
	var cl_I m2 = m*m;
	var cl_LF fterm = cl_I_to_LF(1,actuallen)/m;
	var cl_LF fsum = fterm;
	for (var uintC n = 1; fterm >= eps; n++) {
		fterm = fterm/m2;
		fterm = cl_LF_shortenwith(fterm,eps);
		if ((n % 2) == 0)
			fsum = fsum + LF_to_LF(fterm/(2*n+1),actuallen);
		else
			fsum = fsum - LF_to_LF(fterm/(2*n+1),actuallen);
	}
	return shorten(fsum,len);
}

const cl_LF atan_recip_1b (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m;
	var cl_I fterm = floor1((cl_I)1 << (intDsize*actuallen), m);
	var cl_I fsum = fterm;
	for (var uintC n = 1; fterm > 0; n++) {
		fterm = floor1(fterm,m2);
		if ((n % 2) == 0)
			fsum = fsum + floor1(fterm,2*n+1);
		else
			fsum = fsum - floor1(fterm,2*n+1);
	}
	return scale_float(cl_I_to_LF(fsum,len),-intDsize*(sintC)actuallen);
}

const cl_LF atan_recip_1c (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m;
	var sintC N = (sintC)(0.69314718*intDsize/2*actuallen/log(double_approx(m))) + 1;
	var cl_I num = 0, den = 1; // "lazy rational number"
	for (sintC n = N-1; n>=0; n--) {
		// Multiply sum with 1/m^2:
		den = den * m2;
		// Add (-1)^n/(2n+1):
		if ((n % 2) == 0)
			num = num*(2*n+1) + den;
		else
			num = num*(2*n+1) - den;
		den = den*(2*n+1);
	}
	den = den*m;
	var cl_LF result = cl_I_to_LF(num,actuallen)/cl_I_to_LF(den,actuallen);
	return shorten(result,len);
}

const cl_LF atan_recip_1d (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m;
	var uintC N = (uintC)(0.69314718*intDsize/2*actuallen/log(double_approx(m))) + 1;
	CL_ALLOCA_STACK;
	var cl_I* bv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* qv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var uintC n;
	for (n = 0; n < N; n++) {
		new (&bv[n]) cl_I ((n % 2) == 0 ? (cl_I)(2*n+1) : -(cl_I)(2*n+1));
		new (&qv[n]) cl_I (n==0 ? m : m2);
	}
	var cl_rational_series series;
	series.av = NULL; series.bv = bv;
	series.pv = NULL; series.qv = qv; series.qsv = NULL;
	var cl_LF result = eval_rational_series(N,series,actuallen);
	for (n = 0; n < N; n++) {
		bv[n].~cl_I();
		qv[n].~cl_I();
	}
	return shorten(result,len);
}


// Method 2: atan(1/m) = sum(n=0..infty, 4^n*n!^2/(2n+1)! * m/(m^2+1)^(n+1))

const cl_LF atan_recip_2a (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_LF eps = scale_float(cl_I_to_LF(1,actuallen),-intDsize*(sintC)actuallen);
	var cl_I m2 = m*m+1;
	var cl_LF fterm = cl_I_to_LF(m,actuallen)/m2;
	var cl_LF fsum = fterm;
	for (var uintC n = 1; fterm >= eps; n++) {
		fterm = The(cl_LF)((2*n)*fterm)/((2*n+1)*m2);
		fterm = cl_LF_shortenwith(fterm,eps);
		fsum = fsum + LF_to_LF(fterm,actuallen);
	}
	return shorten(fsum,len);
}

const cl_LF atan_recip_2b (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m+1;
	var cl_I fterm = floor1((cl_I)m << (intDsize*actuallen), m2);
	var cl_I fsum = fterm;
	for (var uintC n = 1; fterm > 0; n++) {
		fterm = floor1((2*n)*fterm,(2*n+1)*m2);
		fsum = fsum + fterm;
	}
	return scale_float(cl_I_to_LF(fsum,len),-intDsize*(sintC)actuallen);
}

const cl_LF atan_recip_2c (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m+1;
	var uintC N = (uintC)(0.69314718*intDsize*actuallen/log(double_approx(m2))) + 1;
	var cl_I num = 0, den = 1; // "lazy rational number"
	for (uintC n = N; n>0; n--) {
		// Multiply sum with (2n)/(2n+1)(m^2+1):
		num = num * (2*n);
		den = den * ((2*n+1)*m2);
		// Add 1:
		num = num + den;
	}
	num = num*m;
	den = den*m2;
	var cl_LF result = cl_I_to_LF(num,actuallen)/cl_I_to_LF(den,actuallen);
	return shorten(result,len);
}

const cl_LF atan_recip_2d (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m+1;
	var uintC N = (uintC)(0.69314718*intDsize*actuallen/log(double_approx(m2))) + 1;
	CL_ALLOCA_STACK;
	var cl_I* pv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* qv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var uintC n;
	new (&pv[0]) cl_I (m);
	new (&qv[0]) cl_I (m2);
	for (n = 1; n < N; n++) {
		new (&pv[n]) cl_I (2*n);
		new (&qv[n]) cl_I ((2*n+1)*m2);
	}
	var cl_rational_series series;
	series.av = NULL; series.bv = NULL;
	series.pv = pv; series.qv = qv; series.qsv = NULL;
	var cl_LF result = eval_rational_series(N,series,actuallen);
	for (n = 0; n < N; n++) {
		pv[n].~cl_I();
		qv[n].~cl_I();
	}
	return shorten(result,len);
}


// Main program: Compute and display the timings.

int main (int argc, char * argv[])
{
	int repetitions = 1;
	if ((argc >= 3) && !strcmp(argv[1],"-r")) {
		repetitions = atoi(argv[2]);
		argc -= 2; argv += 2;
	}
	if (argc < 2)
		exit(1);
	cl_I m = (cl_I)argv[1];
	uintC len = atol(argv[2]);
	cl_LF p;
	ln(cl_I_to_LF(1000,len+10)); // fill cache
	// Method 1.
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atan_recip_1a(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atan_recip_1b(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atan_recip_1c(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atan_recip_1d(m,len); }
	}
	cout << p << endl;
	// Method 2.
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atan_recip_2a(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atan_recip_2b(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atan_recip_2c(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atan_recip_2d(m,len); }
	}
	cout << p << endl;
	// Method 3.
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = The(cl_LF)(atan(cl_RA_to_LF(1/(cl_RA)m,len))); }
	}
	cout << p << endl;
}


// Timings of the above algorithms, on an i486 33 MHz, running Linux.
// m = 390112. (For Jörg Arndt's formula (4.15).)
//    N      1a     1b     1c     1d      2a     2b     2c     2d      3
//    10     0.0027 0.0018 0.0019 0.0019  0.0032 0.0022 0.0019 0.0019  0.0042
//    25     0.0085 0.0061 0.0058 0.0061  0.0095 0.0069 0.0056 0.0061  0.028
//    50     0.024  0.018  0.017  0.017   0.026  0.020  0.016  0.017   0.149
//   100     0.075  0.061  0.057  0.054   0.079  0.065  0.052  0.052   0.71
//   250     0.41   0.33   0.32   0.26    0.42   0.36   0.28   0.24    3.66
//   500     1.57   1.31   1.22   0.88    1.57   1.36   1.10   0.83   13.7
//  1000     6.08   5.14   4.56   2.76    6.12   5.36   4.06   2.58   45.5
//  2500    36.5   32.2   25.8   10.2    38.4   33.6   22.2    9.1   191
//  5000
// 10000
// asymp.    N^2    N^2    N^2    FAST    N^2    N^2    N^2    FAST    FAST
//
// m = 319. (For Jörg Arndt's formula (4.7).)
//    N      1a     1b     1c     1d      2a     2b     2c     2d      3
//  1000     6.06   4.40   9.17   3.82    5.29   3.90   7.50   3.53   50.3
//
// m = 18. (For Jörg Arndt's formula (4.4).)
//    N      1a     1b     1c     1d      2a     2b     2c     2d      3
//  1000    11.8    9.0   22.3    6.0    10.2    7.7   17.1    5.7    54.3
