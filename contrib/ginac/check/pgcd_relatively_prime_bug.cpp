/** @file pgcd_relatively_prime_bug.cpp
 *
 * A program exposing historical bug in the pgcd() function. 
 */
#include <string>
#include <iostream>
#include <utility>
#include "ginac.h"
using namespace std;
using namespace GiNaC;

int main(int argc, char** argv)
{
	cout << "Checking for pgcd() bug regarding relatively prime polynomials: " << flush;
	const symbol q("q");
	parser reader;
	reader.get_syms().insert(make_pair(string("q"), q));

	ex t = reader("-E20^16*E4^8*E5^8*E1^16*q^4"
		      "-(E10^24-E20^8*E5^16)*E4^16*E1^8"
		      "+E2^24*E20^16*E5^8*q^4");
	ex g = gcd(t.expand(), t.diff(q).expand()) - 1;
	if (!g.is_zero()) {
		cout << " oops!" << endl <<
			"** Error: should be 0, got " << g << endl << flush;
		throw std::logic_error("gcd was miscalculated");
	}
	cout << "not found" << endl << flush;
	return 0;
}

