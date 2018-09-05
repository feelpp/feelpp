// Input expression containing variable 'x' and compute its derivative
// with respect to 'x'.
// Example from the tutorial (chapter Input/Output, section `Expression
// input').
#include <iostream>
#include <string>
#include <stdexcept>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

int main()
{
	cout << "Enter an expression containing 'x': " << flush;
	parser reader;

	try {
		ex e = reader(cin);
		symtab table = reader.get_syms();

		symbol x = table.find("x") != table.end() ? 
			   ex_to<symbol>(table["x"]) : symbol("x");

		cout << "The derivative of " << e << " with respect to x is ";
		cout << e.diff(x) << "." << endl;
	} catch (exception &p) {
		cerr << p.what() << endl;
	}
}

