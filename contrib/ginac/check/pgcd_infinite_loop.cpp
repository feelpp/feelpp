#include <iostream>
#include <string>
#include "ginac.h"
using namespace GiNaC;
using namespace std;

static const string srep("\
792*z^8*w^4*x^3*y^4*u^7 + 24*z^4*w^4*x^2*y^3*u^4	\
+ 264*z^8*w^3*x^2*y^7*u^5 + 198*z^4*w^5*x^5*y*u^6	\
+ 110*z^2*w^3*x^5*y^4*u^6 - 120*z^8*w*x^4*u^6		\
- 480*z^5*w*x^4*y^6*u^8 - 720*z^7*x^3*y^3*u^7		\
+ 165*z^4*w^2*x^4*y*u^5 + 450*z^8*w^6*x^2*y*u^8		\
+ 40*z^2*w^3*x^3*y^3*u^6 - 288*z^7*w^2*x^3*y^6*u^6	\
+ 250*z^6*w^4*x^2*y^4*u^8 + 576*z^7*w^7*x^2*y^4*u^8	\
- 80*z^6*w^2*x^5*y^3*u^7 - 144*z^8*w^4*x^5*u^7		\
+ 120*z^4*w*x^2*y^6*u^6 + 320*z^5*w^5*x^2*y^7*u^8	\
+ 192*z^7*w^6*x*y^7*u^6 - 12*z^4*w^3*x^3*y^5*u^6	\
- 36*z^4*w^4*x^4*y^2*u^8 + 72*z^4*w^5*x^3*u^6		\
- 20*z^2*w^2*x^4*y^5*u^8 + 660*z^8*w*x^2*y^4*u^6	\
+ 66*z^4*w^4*x^4*y^4*u^4 + 440*z^6*w^2*x^3*y^7*u^7	\
- 30*z^4*w*x^3*y^2*u^7 - 48*z^8*w^3*x^4*y^3*u^5		\
+ 72*z^6*w^2*x*y^6*u^4 - 864*z^7*w^3*x^4*y^3*u^8	\
+ 480*z^7*w^4*x*y^4*u^7 + 60*z^4*w^2*x^2*u^5		\
+ 375*z^8*w^3*x*y*u^7 + 150*z^8*w^5*x*y^4*u^6		\
+ 180*z^6*x*y^3*u^5 + 216*z^6*w^3*x^2*y^3*u^6");

int main(int argc, char** argv)
{
	cout << "Checking for more pgcd() bugs (infinite loop, miscalculation) ... " << flush;
	parser the_parser;
	ex e = the_parser(srep);
	const symbol x = ex_to<symbol>(the_parser.get_syms()["x"]);
	ex g = gcd(e, e.diff(x));
	ex should_be = the_parser(string("u^4*z^2"));
	if (!(g-should_be).expand().is_zero()) {
		cout << "GCD was miscomputed. " << flush;
		return 1;
	}
	cout << "not found. " << flush;
	return 0;
}
