#include <ctime>
#include <iostream>
using namespace std;
#include <ginac/ginac.h>
using namespace GiNaC;

/*
 * Demonstrates the use of compile_ex.
 *
 * Compiles a small expression as C code via compile_ex and evaluates the
 * expression numerically. The evalation speed is timed and compared to the
 * evaluation of the original GiNaC expression.
 *
 */

int main()
{
	// Define some expression
	symbol x("x");
	ex expr = sin(x);
 
	// Some variables for timing
	time_t start, end;
	double cpu_time_used;

	// Our function pointer that points to the compiled ex
	FUNCP_1P fp;
	compile_ex(expr, x, fp);

	// Do some (not necessarily meaningful ;-)) numerical stuff ...

	cout << "Doing numerics with compile_ex ..." << endl;

	// First using compile_ex
	{
		double result;
		double point = 0.2;
		start = clock();
		for (int i=0; i<100000; ++i) {
			point += 0.001;
			result += fp(point);
		}
		end = clock();

		// Show the result
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		cout << "result = " << result << " in " << cpu_time_used << " seconds" << endl;
	}

	cout << "Doing numerics without compile_ex ..." << endl;

	// Then without compile_ex
	{
		ex result;
		ex point = 0.2;
		start = clock();
		for (int i=0; i<100000; ++i) {
			point += 0.001;
			result += sin(point);
		}
		end = clock();

		// Show the other result
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		cout << "result = " << result << " in " << cpu_time_used << " seconds" << endl;
	}

	return 0;
}
