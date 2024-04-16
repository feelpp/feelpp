// This program prints the largest now known even perfect number.

#include <cln/integer.h>
#include <cln/integer_io.h>

using namespace std;
using namespace cln;

int main ()
{
	// previous ones were 1257787, 1398269, 2976221, 3021377, 6972593,
	// 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667,
	// 42643801, 43112609, 57885161, 74207281, 77232917, 82589933
	int p = 82589933;
	cl_I x = (((cl_I)1 << p) - 1) << (p-1);
	cout << x << endl;
}
