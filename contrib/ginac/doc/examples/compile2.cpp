#include <iostream>
using namespace std;
#include <ginac/ginac.h>
using namespace GiNaC;
// Yes, we are using CUBA (should be installed on the system!)
#include <cuba.h>

/*
 * Demonstrates the use of compile_ex with the CUBA library.
 *
 * The user can enter an expression on the command line. This expression is
 * compiled via compile_ex and integrated over the region 0 <= x,y <= 1 with
 * the help of the CUBA library (http://www.feynarts.de/cuba).
 *
 */

int main()
{
	// Let the user enter a expression
	symbol x("x"), y("y");
	string s;
	cout << "Enter an expression containing 'x' and/or 'y': ";
	cin >> s;
	// Expression now in expr
	ex expr(s, lst(x,y));
 
	cout << "start integration of " << expr << " ..." << endl;

	// Some definitions for VEGAS 
	#define NDIM 2
	#define NCOMP 1
	#define EPSREL 1e-3
	#define EPSABS 1e-12
	#define VERBOSE 0
	#define MINEVAL 0
	#define MAXEVAL 50000
	#define NSTART 1000
	#define NINCREASE 500

	// Some variables for VEGAS
	int comp, nregions, neval, fail;
	double integral[NCOMP], error[NCOMP], prob[NCOMP];

	// Our function pointer that points to the compiled ex
	FUNCP_CUBA fp;
	compile_ex(lst(expr), lst(x,y), fp);

	// Starting VEGAS
	// By invocation of compile() the expression in expr is converted into the
	// appropriate function pointer
	Vegas(NDIM, NCOMP, fp,
	      EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NSTART, NINCREASE,
	      &neval, &fail, integral, error, prob);

	// Show the result
	cout << "result: " << integral[0] << endl;

	return 0;
}
