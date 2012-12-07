// This program prints the largest now known even perfect number.

#include <cln/integer.h>
#include <cln/integer_io.h>

using namespace std;
using namespace cln;

int main ()
{
	// previous ones were 1257787, 1398269, 2976221, 3021377, 6972593,
	// 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667
	int p = 43112609;
	cl_I x = (((cl_I)1 << p) - 1) << (p-1);
	cout << x << endl;
}
