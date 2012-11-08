#include <cstdlib>
#include <cln/io.h>
#include <cln/real.h>

using namespace std;
using namespace cln;

extern int test_integer();
extern int test_rational();
extern int test_sfloat();
extern int test_ffloat();
extern int test_dfloat();
extern int test_lfloat();

int test_elementary (void)
{
	int error = 0;
	error |= test_integer();
	error |= test_rational();
	error |= test_sfloat();
	error |= test_ffloat();
	error |= test_dfloat();
	error |= test_lfloat();
	return error;
}

extern int test_gcd (void);
extern int test_xgcd (void);
extern int test_sqrtp (void);

int test_all (void)
{
	int error = 0;
	error |= test_elementary();
	error |= test_gcd();
	error |= test_xgcd();
	error |= test_sqrtp();
	return error;
}

int main ()
{
	if (!test_all()) {
		cout << "Tests passed." << endl;
		exit(0);
	} else {
		cout << "Tests failed" << endl;
		exit(1);
	}
}
