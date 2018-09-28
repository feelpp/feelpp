// This program prints the smallest probable prime >= x, x being given on the
// command line.

// We work with real numbers and integers.
#include <cln/real.h>
#include <cln/integer.h>

// We do I/O.
#include <cln/io.h>
#include <cln/integer_io.h>

// The function nextprobprime() is part of the number theory package.
#include <cln/numtheory.h>

int main (int argc, char* argv[])
{
	if (argc != 2) {
		std::cerr << "Usage: nextprime x" << std::endl;
		return(1);
	}
	cln::cl_R x = (cln::cl_R)argv[1];
	cln::cl_I p = cln::nextprobprime(x);
	std::cout << p << std::endl;
	return(0);
}
