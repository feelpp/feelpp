/**
 * @file factor_upoly_q_bug.cpp Check for a bug in factor_univariate().
 *
 * factor_univariate() didn't check if the argument is an integer polynomial,
 * the result was a segfault.
 */
#include "ginac.h"
#include <iostream>
using namespace GiNaC;
using namespace std;

int main(int argc, char** argv)
{
	cout << "checking if factor() handles rational polynomials. ";
	parser p;
	ex e = p("174247781*x^2-1989199947807987/200000000000000*x");
	ex ef = factor(e);
	ex diff = (e - ef).expand();
	cout << "yes" << endl;
	return 0;
}
