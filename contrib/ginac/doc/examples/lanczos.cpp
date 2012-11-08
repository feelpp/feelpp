/*
 * This program can be used to find the coefficients needed to approximate
 * the gamma function using the Lanczos approximation.
 *
 * Usage: lanczos -n order -D digits
 *
 * The order defaults to 10. digits defaults to GiNaCs default.
 * It is recommended to run the program several times with an increasing
 * value for digits until numerically stablilty for the values has been
 * reached. The program will also print the number of digits for which
 * the approximation is still reliable. This will be lower then the
 * number of digits given at command line. It is determined by comparing
 * Gamma(1/2) to sqrt(Pi). Note that the program may crash if the number of
 * digits is unreasonably small for the given order. Another thing that
 * can happen if the number of digits is too small is that the program will
 * print "Forget it, this is waaaaaaay too inaccurate." at the top of the
 * output.
 *
 * The gamma function can be (for real_part(z) > -1/2) calculated using
 *
 *    Gamma(z+1) = sqrt(2*Pi)*power(z+g+ex(1)/2, z+ex(1)/2)*exp(-(z+g+ex(1)/2))
 *                  *A_g(z),
 * where,
 * 
 *    A_g(z) = coeff[0] + coeff[1]/(z+1) + coeff[2]/(z+2) + ...
 *                + coeff[N-1]/(z+N-1).
 *
 * The value of g is taken to be equal to the order N.
 *
 * More details can be found at Wikipedia:
 * http://en.wikipedia.org/wiki/Lanczos_approximation.
 *
 * (C) 2006 Chris Dams
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301 USA
 */

#include <vector>
#include <cstddef> // for size_t
#include <iostream>
#include <cln/io.h>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/rational.h>
#include <cln/real.h>
#include <cln/real_io.h>
#include <cln/complex.h>

using namespace std;
using namespace cln;

/*
 * Chebyshev polynomial coefficient matrix as far as is required for
 * the Lanczos approximation.
 */
void calc_chebyshev(vector<vector<cl_I> >& C, const size_t size)
{	
	C.reserve(size);
	for (size_t i=0; i<size; ++i)
		C.push_back(vector<cl_I>(size, 0));
	C[1][1] = cl_I(1);
	C[2][2] = cl_I(1);
	for (size_t i=3; i<size; i+=2)
		C[i][1] = -C[i-2][1];
	for (size_t i=3; i<size; ++i)
		C[i][i] = 2*C[i-1][i-1];
	for (size_t j=2; j<size; ++j)
		for (size_t i=j+2; i<size; i+=2)
			C[i][j] = 2*C[i-1][j-1] - C[i-2][j];
}

/*
 * The coefficients p_n(g) that occur in the Lanczos approximation.
 */
const cl_F p(const size_t k, const cl_F& g, const vector<vector<cln::cl_I> >& C)
{	
	const float_format_t prec = float_format(g);
	const cl_F sqrtPi = sqrt(pi(prec));
	cl_F result = cl_float(0, prec);

	result = result + 
		(2*C[2*k+1][1])*recip(sqrtPi)*
		recip(sqrt(2*g+1))*exp(g+cl_RA(1)/2);

	for (size_t a=1; a <= k; ++a)
		result = result + 
			(2*C[2*k+1][2*a+1]*doublefactorial(2*a-1))*recip(sqrtPi)*
			recip(expt(2*cl_I(a)+2*g+1, a))*recip(sqrt(2*cl_I(a)+2*g+1))*
			exp(cl_I(a)+g+cl_RA(1)/2);
	return result;
}

/*
 * Calculate coefficients that occurs in the expression for the
 * Lanczos approximation of order n.
 */
void calc_lanczos_coeffs(vector<cl_F>& lanc, const cln::cl_F& g)
{	
	const size_t n = lanc.size();
	vector<vector<cln::cl_I> > C;
	calc_chebyshev(C, 2*n+2);
	
	// \Pi_{i=1}^n (z-i+1)/(z+i) = \Pi_{i=1}^n (1 - (2i-1)/(z+i))
	// Such a product can be rewritten as multivariate polynomial \in
	// Q[1/(z+1),...1/(z+n)] of degree 1. To store coefficients of this
	// polynomial we use vector<cl_I>, so the set of such polynomials is
	// stored as
	vector<vector<cln::cl_I> > fractions(n);
	
	// xi = 1/(z+i)
	fractions[0] = vector<cl_I>(1);
	fractions[0][0] = cl_I(1);
	fractions[1] = fractions[0];
	fractions[1].push_back(cl_I(-1)); // 1 - 1/(z+1)

	for (size_t i=2; i<n; ++i)
	{	
		// next = previous*(1 - (2*i-1)*xi) = 
		// previous - (2*i-1)*xi - (2i-1)*(previous-1)*xi
		// Such representation is convenient since previous contains only
		// only x1...x[i-1]

		fractions[i] = fractions[i-1];
		fractions[i].push_back(-cl_I(2*i-1)); // - (2*i-1)*xi
		// Now add -(2i+1)*(previous-1)*xi using formula
		// 1/(z+j)*1/(z+i) = 1/(i-j)*(1/(z+j) - 1/(z+i)), that is
		// xj*xi = 1/(i-j)(xj - xi)
		for (size_t j=1; j < i; j++) {
			cl_I curr = - exquo(cl_I(2*i-1)*fractions[i-1][j], cl_I(i-j));
			fractions[i][j] = fractions[i][j] + curr;
			fractions[i][i] = fractions[i][i] - curr; 
			// N.B. i-th polynomial has (i+1) non-zero coefficients
		}
	}

	vector<cl_F> p_cache(n);
	for (size_t i=0; i < n; i++)
		p_cache[i] = p(i, g, C);

	lanc[0] = p_cache[0]/2;

	float_format_t prec = float_format(g);
	// A = p(i, g, C)*fraction[i] = p(i, g, C)*F[i][j]*xj,
	for (size_t j=1; j < n; ++j) {
		lanc[j] = cl_float(0, prec);
		lanc[0] = lanc[0] + p_cache[j];
		for (size_t i=j; i < n; i++)
			lanc[j] = lanc[j] + p_cache[i]*fractions[i][j];
	}
}

/*
 * Calculate Gamma(z) using the Lanczos approximation with parameter g and
 * coefficients stored in the exvector coeffs.
 */
const cl_N calc_gamma(const cl_N& z, const cl_F& g, const vector<cl_F>& lanc)
{	
	const cl_F thePi = pi(float_format(g)); 
	// XXX: need to check if z is floating-point and precision of z
	// matches one of g
	if (realpart(z) < 0.5)
		return (thePi/sin(thePi*z))/calc_gamma(1-z, g, lanc);
	cl_N A = lanc[0];
	for (size_t i=1; i < lanc.size(); ++i)
		A = A + lanc[i]/(z-1+i);
	cl_N result = sqrt(2*thePi)*expt(z+g-cl_RA(1)/2, z-cl_RA(1)/2)*
		            exp(-(z+g-cl_RA(1)/2))*A;
	return result;
}

void usage(char *progname)
{	cout << "Usage: " << progname << " -n order -D digits" << endl;
	exit(0);
}

void read_options(int argc, char**argv, int &order, int& digits)
{  int c;
   while((c=getopt(argc,argv,"n:D:"))!=-1)
   {  if(c=='n')
         order = atoi(optarg);
      else if (c=='D')
         digits = atoi(optarg);
      else
         usage(argv[0]);
   }
   if(optind!=argc)
      usage(argv[0]);
}

int main(int argc, char *argv[])
{	
/*
 * Handle command line options.
 */
	int order = 10;
	int digits_ini = 17;
	read_options(argc, argv, order, digits_ini);
	float_format_t prec = float_format(digits_ini);
	const cl_F thePi = pi(prec);
/*
 * Calculate coefficients.
 */
	const cl_F g_val = cl_float(order, prec);
	vector<cl_F> coeffs(order);
	calc_lanczos_coeffs(coeffs, g_val);
/*
 * Determine the accuracy by comparing Gamma(1/2) to sqrt(Pi).
 */
	cl_N gamma_half = calc_gamma(cl_float(cl_RA(1)/2, prec), g_val, coeffs);
	cl_F digits = (ln(thePi)/2 - ln(abs(gamma_half - sqrt(thePi))))/ln(cl_float(10, prec));
	int i_digits = cl_I_to_int(floor1(digits));
	if (digits < cl_I(1))
		cout << "Forget it, this is waaaaaaay too inaccurate." << endl;
	else
		cout << "Reliable digits: " << i_digits << endl;
/*
 * Print the coefficients.
 */
	for (size_t i=0; i<coeffs.size(); ++i) {
		cout << "coeffs_" << order << "[" << i << "] = \"";
		cout << coeffs[i];
		cout << "\";" << endl;
	}
	return 0;
}
