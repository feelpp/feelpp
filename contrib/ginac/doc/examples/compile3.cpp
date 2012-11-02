#include <ctime>
#include <iostream>
using namespace std;
#include <ginac/ginac.h>
using namespace GiNaC;

/*
 * Demonstrates the use of link_ex.
 *
 * When run for the first time link_ex will fail. This is a rude way of
 * checking whether the needed .so file is available.  The .so is then created
 * by compile_ex using the filename parameter. When run again link_ex will use
 * the existing .so file.
 *
 */
int main()
{
	FUNCP_2P fp;
	try {
		link_ex("compile3_testprg.so", fp);
		cout << "Using existing 'compile3_testprg.so'." << endl;
	}
	catch (const std::exception& e) {
		// hope the exception is just raised because of missing 'compile2_testprg.so' file,
		// so being lazy no error management here ...
		cout << "Error: " << e.what() << endl;
		cout << "Building new 'compile3_testprg.so'." << endl;
		symbol a, b;
		ex expr = a*b;
		compile_ex(expr, a, b, fp, "compile3_testprg");
	}

	cout << "result of 2.3*1.5 is " << fp(2.3, 1.5) << endl;

	return 0;
}
