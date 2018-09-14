#ifndef _EXAM_H
#define _EXAM_H

#include <cln/number.h>
#include <cln/io.h>
using namespace std;
using namespace cln;

// Michael Stoll  23. 3. 1993
// C++ version: Bruno Haible 1.11.1995

struct plus_test {
	const char * arg1;
	const char * arg2;
	const char * result;
};

struct minus_test {
	const char * arg1;
	const char * arg2;
	const char * result;
};

struct mul_test {
	const char * arg1;
	const char * arg2;
	const char * result;
};

struct floor_test {
	const char * arg1;
	const char * arg2;
	const char * result1;
	const char * result2;
};

struct div_test {
	const char * arg1;
	const char * arg2;
	const char * result;
};

#define num_elements(array)  (sizeof(array)/sizeof(array[0]))

#define DO_BINOP_TEST(typename,type,rtype,opname,op)  \
static int test_##typename##_##opname (void)				\
{									\
	int error = 0;							\
	for (unsigned int i = 0; i < num_elements(typename##_##opname##_tests); i++) { \
		opname##_test& test = typename##_##opname##_tests[i];	\
		type arg1 = type(test.arg1);				\
		type arg2 = type(test.arg2);				\
		rtype computed_result = arg1 op arg2;			\
		rtype result = rtype(test.result);			\
		if (computed_result != result) {			\
			std::cerr << "Error in " #typename "_" #opname "_tests[" << i << "] !" << endl; \
			std::cerr << "Result should be: " << result << endl;	\
			std::cerr << "Result computed : " << computed_result << endl << endl;	\
			error = 1;					\
		}							\
	}								\
	return error;							\
}

#define DO_FLOOR_TEST(typename,type)  \
static int test_##typename##_floor (void)				\
{									\
	int error = 0;							\
	for (unsigned int i = 0; i < num_elements(typename##_floor_tests); i++) { \
		floor_test& test = typename##_floor_tests[i];		\
		type arg1 = type(test.arg1);				\
		type arg2 = type(test.arg2);				\
		type##_div_t computed_result = floor2(arg1,arg2);	\
		cl_I result1 = cl_I(test.result1);			\
		type result2 = type(test.result2);			\
		if ((computed_result.quotient != result1) || (computed_result.remainder != result2)) { \
			std::cerr << "Error in " #typename "_floor_tests[" << i << endl; \
			std::cerr << "Results should be: " << result1 << ", " << result2 << endl;	\
			std::cerr << "Results computed : " << computed_result.quotient << ", " << computed_result.remainder << endl << endl;	\
			error = 1;					\
		}							\
	}								\
	return error;							\
}

#define DO_TESTS(typename,type,qtype)  \
  DO_BINOP_TEST(typename,type,type,plus,+)				\
  DO_BINOP_TEST(typename,type,type,minus,-)				\
  DO_BINOP_TEST(typename,type,type,mul,*)				\
  DO_FLOOR_TEST(typename,type)						\
  DO_BINOP_TEST(typename,type,qtype,div,/)				\
int test_##typename (void)						\
{									\
	int error = 0;							\
	error |= test_##typename##_plus();				\
	error |= test_##typename##_minus();				\
	error |= test_##typename##_mul();				\
	error |= test_##typename##_floor();				\
	error |= test_##typename##_div();				\
	return error;							\
}

#endif /* _EXAM_H */
