#include <iostream>

// Elementary operations.
extern int test_I_abs (int iterations);
extern int test_I_compare (int iterations);
extern int test_I_plus (int iterations);
extern int test_I_minus (int iterations);
extern int test_I_plus1 (int iterations);
extern int test_I_minus1 (int iterations);
extern int test_I_mul (int iterations);
extern int test_I_div (int iterations);
// Euclidean ring operations.
extern int test_I_gcd (int iterations);
extern int test_I_xgcd (int iterations);
// Bit vector operations.
extern int test_I_ash (int iterations);
extern int test_I_evenp (int iterations);
extern int test_I_oddp (int iterations);
extern int test_I_lognot (int iterations);
extern int test_I_logand (int iterations);
extern int test_I_logandc1 (int iterations);
extern int test_I_logandc2 (int iterations);
extern int test_I_logior (int iterations);
extern int test_I_logorc1 (int iterations);
extern int test_I_logorc2 (int iterations);
extern int test_I_logxor (int iterations);
extern int test_I_lognand (int iterations);
extern int test_I_lognor (int iterations);
extern int test_I_logeqv (int iterations);
extern int test_I_boole (int iterations);
extern int test_I_logbitp (int iterations);
extern int test_I_logtest (int iterations);
extern int test_I_ldb (int iterations);
extern int test_I_ldb_test (int iterations);
extern int test_I_mask_field (int iterations);
extern int test_I_dpb (int iterations);
extern int test_I_deposit_field (int iterations);
extern int test_I_logcount (int iterations);
extern int test_I_integer_length (int iterations);
extern int test_I_ord2 (int iterations);
extern int test_I_power2p (int iterations);
// More complex operations.
extern int test_I_isqrt (int iterations);
extern int test_I_sqrtp (int iterations);
// Miscellaneous.
extern int test_I_io (int iterations);
extern int test_I_GV (int iterations);

#define RUN(tester,iterations)  \
	std::cout << "Testing "#tester"..." << std::endl; \
	error |= tester (iterations);


int test_I (int iterations)
{
	int error = 0;
	// Elementary operations.
	RUN(test_I_abs,iterations);
	RUN(test_I_compare,iterations);
	RUN(test_I_plus,iterations);
	RUN(test_I_minus,iterations);
	RUN(test_I_plus1,iterations);
	RUN(test_I_minus1,iterations);
	RUN(test_I_mul,iterations);
	RUN(test_I_div,iterations);
	// Euclidean ring operations.
	RUN(test_I_gcd,iterations);
	RUN(test_I_xgcd,iterations);
	// Bit vector operations.
	RUN(test_I_ash,iterations);
	RUN(test_I_evenp,iterations);
	RUN(test_I_oddp,iterations);
	RUN(test_I_lognot,iterations);
	RUN(test_I_logand,iterations);
	RUN(test_I_logandc1,iterations);
	RUN(test_I_logandc2,iterations);
	RUN(test_I_logior,iterations);
	RUN(test_I_logorc1,iterations);
	RUN(test_I_logorc2,iterations);
	RUN(test_I_logxor,iterations);
	RUN(test_I_lognand,iterations);
	RUN(test_I_lognor,iterations);
	RUN(test_I_logeqv,iterations);
	RUN(test_I_boole,iterations);
	RUN(test_I_logbitp,iterations);
	RUN(test_I_logtest,iterations);
	RUN(test_I_ldb,iterations);
	RUN(test_I_ldb_test,iterations);
	RUN(test_I_mask_field,iterations);
	RUN(test_I_dpb,iterations);
	RUN(test_I_deposit_field,iterations);
	RUN(test_I_logcount,iterations);
	RUN(test_I_integer_length,iterations);
	RUN(test_I_ord2,iterations);
	RUN(test_I_power2p,iterations);
	// More complex operations.
	RUN(test_I_isqrt,iterations);
	RUN(test_I_sqrtp,iterations);
	// Miscellaneous.
	RUN(test_I_io,iterations);
	RUN(test_I_GV,iterations);
	return error;
}
