// Computation of artanh(1/m) (m integer) to high precision.

#include "cln/integer.h"
#include "cln/rational.h"
#include "cln/real.h"
#include "cln/complex.h"
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


// Method 1: atanh(1/m) = sum(n=0..infty, 1/(2n+1) * 1/m^(2n+1))
// Method 2: atanh(1/m) = sum(n=0..infty, (-4)^n*n!^2/(2n+1)! * m/(m^2-1)^(n+1))
// a. Using long floats.                     [N^2]
// b. Simulating long floats using integers. [N^2]
// c. Using integers, no binary splitting.   [N^2]
// d. Using integers, with binary splitting. [FAST]
// Method 3: general built-in algorithm.     [FAST]
// Method 4: atanh(x) = 1/2 ln((1+x)/(1-x)),
// using the general built-in algorithm      [FAST]


// Method 1: atanh(1/m) = sum(n=0..infty, 1/(2n+1) * 1/m^(2n+1))

const cl_LF atanh_recip_1a (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_LF eps = scale_float(cl_I_to_LF(1,actuallen),-intDsize*(sintC)actuallen);
	var cl_I m2 = m*m;
	var cl_LF fterm = cl_I_to_LF(1,actuallen)/m;
	var cl_LF fsum = fterm;
	for (var uintC n = 1; fterm >= eps; n++) {
		fterm = fterm/m2;
		fterm = cl_LF_shortenwith(fterm,eps);
		fsum = fsum + LF_to_LF(fterm/(2*n+1),actuallen);
	}
	return shorten(fsum,len);
}

const cl_LF atanh_recip_1b (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m;
	var cl_I fterm = floor1((cl_I)1 << (intDsize*actuallen), m);
	var cl_I fsum = fterm;
	for (var uintC n = 1; fterm > 0; n++) {
		fterm = floor1(fterm,m2);
		fsum = fsum + floor1(fterm,2*n+1);
	}
	return scale_float(cl_I_to_LF(fsum,len),-intDsize*(sintC)actuallen);
}

const cl_LF atanh_recip_1c (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m;
	var sintC N = (sintC)(0.69314718*intDsize/2*actuallen/log(double_approx(m))) + 1;
	var cl_I num = 0, den = 1; // "lazy rational number"
	for (sintC n = N-1; n>=0; n--) {
		// Multiply sum with 1/m^2:
		den = den * m2;
		// Add 1/(2n+1):
		num = num*(2*n+1) + den;
		den = den*(2*n+1);
	}
	den = den*m;
	var cl_LF result = cl_I_to_LF(num,actuallen)/cl_I_to_LF(den,actuallen);
	return shorten(result,len);
}

const cl_LF atanh_recip_1d (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m;
	var uintC N = (uintC)(0.69314718*intDsize/2*actuallen/log(double_approx(m))) + 1;
	CL_ALLOCA_STACK;
	var cl_I* bv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* qv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var uintC n;
	for (n = 0; n < N; n++) {
		new (&bv[n]) cl_I ((cl_I)(2*n+1));
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


// Method 2: atanh(1/m) = sum(n=0..infty, (-4)^n*n!^2/(2n+1)! * m/(m^2-1)^(n+1))

const cl_LF atanh_recip_2a (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_LF eps = scale_float(cl_I_to_LF(1,actuallen),-intDsize*(sintC)actuallen);
	var cl_I m2 = m*m-1;
	var cl_LF fterm = cl_I_to_LF(m,actuallen)/m2;
	var cl_LF fsum = fterm;
	for (var uintC n = 1; fterm >= eps; n++) {
		fterm = The(cl_LF)((2*n)*fterm)/((2*n+1)*m2);
		fterm = cl_LF_shortenwith(fterm,eps);
		if ((n % 2) == 0)
			fsum = fsum + LF_to_LF(fterm,actuallen);
		else
			fsum = fsum - LF_to_LF(fterm,actuallen);
	}
	return shorten(fsum,len);
}

const cl_LF atanh_recip_2b (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m-1;
	var cl_I fterm = floor1((cl_I)m << (intDsize*actuallen), m2);
	var cl_I fsum = fterm;
	for (var uintC n = 1; fterm > 0; n++) {
		fterm = floor1((2*n)*fterm,(2*n+1)*m2);
		if ((n % 2) == 0)
			fsum = fsum + fterm;
		else
			fsum = fsum - fterm;
	}
	return scale_float(cl_I_to_LF(fsum,len),-intDsize*(sintC)actuallen);
}

const cl_LF atanh_recip_2c (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m-1;
	var uintC N = (uintC)(0.69314718*intDsize*actuallen/log(double_approx(m2))) + 1;
	var cl_I num = 0, den = 1; // "lazy rational number"
	for (uintC n = N; n>0; n--) {
		// Multiply sum with -(2n)/(2n+1)(m^2+1):
		num = num * (2*n);
		den = - den * ((2*n+1)*m2);
		// Add 1:
		num = num + den;
	}
	num = num*m;
	den = den*m2;
	var cl_LF result = cl_I_to_LF(num,actuallen)/cl_I_to_LF(den,actuallen);
	return shorten(result,len);
}

const cl_LF atanh_recip_2d (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m-1;
	var uintC N = (uintC)(0.69314718*intDsize*actuallen/log(double_approx(m2))) + 1;
	CL_ALLOCA_STACK;
	var cl_I* pv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* qv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var uintC n;
	new (&pv[0]) cl_I (m);
	new (&qv[0]) cl_I (m2);
	for (n = 1; n < N; n++) {
		new (&pv[n]) cl_I (-(cl_I)(2*n));
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
	    { p = atanh_recip_1a(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atanh_recip_1b(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atanh_recip_1c(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atanh_recip_1d(m,len); }
	}
	cout << p << endl;
	// Method 2.
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atanh_recip_2a(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atanh_recip_2b(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atanh_recip_2c(m,len); }
	}
	cout << p << endl;
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = atanh_recip_2d(m,len); }
	}
	cout << p << endl;
	// Method 3.
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = The(cl_LF)(atanh(cl_RA_to_LF(1/(cl_RA)m,len))); }
	}
	cout << p << endl;
	// Method 4.
	{ CL_TIMING;
	  for (int rep = repetitions; rep > 0; rep--)
	    { p = The(cl_LF)(scale_float(ln(cl_RA_to_LF((cl_RA)(m+1)/(cl_RA)(m-1),len)),-1)); }
	}
	cout << p << endl;
}


// Timings of the above algorithms, on an i486 33 MHz, running Linux.
// m = 3 -> 1/2 ln(2)
//   N    1a     1b     1c     1d     2a     2b     2c     2d      3
//   10   0.021  0.014  0.019  0.012  0.029  0.015  0.023  0.015   0.015
//   25   0.060  0.041  0.073  0.041  0.082  0.046  0.086  0.051   0.066
//   50   0.164  0.110  0.258  0.120  0.203  0.124  0.295  0.142
//  100   0.49   0.35   1.05   0.37   0.60   0.35   1.19   0.42
//  250   2.5    1.9    7.2    1.7    2.9    1.9    8.0    1.8
//  500  10.1    7.2   33.4    5.5   10.7    7.3   36.5    5.9
// 1000  38     30    145     16.1   39     29    158     16.8
// 2500 231    188    976     53    237    186   1081     58
// asymp. N^2    N^2    N^2    FAST   N^2    N^2    N^2    FAST
//
// m = 9 -> 1/2 ln(5/4)
//   N    1a     1b     1c     1d     2a     2b     2c     2d     3      4
//   10   0.0106 0.0072 0.0084 0.0061 0.0139 0.0073 0.0098 0.0073 0.0140 0.0211
//   25   0.031  0.021  0.029  0.019  0.039  0.022  0.031  0.022  0.063  0.081
//   50   0.083  0.057  0.091  0.056  0.098  0.058  0.098  0.060  0.232  0.212
//  100   0.25   0.17   0.32   0.16   0.28   0.17   0.34   0.17   0.60   0.59
//  250   1.28   0.94   2.11   0.77   1.40   0.91   2.18   0.76   2.76   2.76
//  500   5.1    3.6    9.4    2.5    5.2    3.4    9.3    2.4   10.4    9.7
// 1000  19.1   14.7   42      7.8   18.5   13.6   42      7.4   31     30
// 2500 116     93    279     29.6  113     86    278     30.0  129    125
// asymp. N^2    N^2    N^2    FAST   N^2    N^2    N^2    FAST   FAST   FAST
