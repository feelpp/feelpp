#include <cln/integer.h>
#include <cln/integer_io.h>
#include <cln/io.h>

using namespace std;
using namespace cln;

#define ASSERT(expr)  \
  if (!(expr)) {                                        \
	std::cerr << "Assertion failed! File " << __FILE__ << ", line " << __LINE__ << endl;    \
	error = 1;                                      \
  }

struct gcd_test {
	const char * arg1;
	const char * arg2;
	const char * result;
};

#define num_elements(array)  (sizeof(array)/sizeof(array[0]))

// Note: This macro differs slightly from the one in "exam.h".
#define DO_BINOP_TEST(typename,type,rtype,opname)  \
static int test_##typename##_##opname (void)				\
{									\
	int error = 0;							\
	for (unsigned int i = 0; i < num_elements(typename##_##opname##_tests); i++) { \
		opname##_test& test = typename##_##opname##_tests[i];	\
		type arg1 = type(test.arg1);				\
		type arg2 = type(test.arg2);				\
		rtype computed_result = opname(arg1,arg2);		\
		rtype result = rtype(test.result);			\
		if (computed_result != result) {			\
			std::cerr << "Error in " #typename "_" #opname "_tests[" << i << "] !" << endl;	\
			std::cerr << "Result should be: " << result << endl;	\
			std::cerr << "Result computed : " << computed_result << endl << endl;	\
			error = 1;					\
		}							\
	}								\
	return error;							\
}

static gcd_test integer_gcd_tests[] = {
	{ "123456789", "345", "3" },
	{ "345", "123456789", "3" },
	{ "10", "0", "10" },
	{ "0", "10", "10" },
	{ "2523533737", "855322739", "1" },
	{ "855322739", "2523533737", "1" },
	{ "101611479673163974026724715741235467160607959655653420075620", "533177863832047932237026621580126811198495699416238676294977", "1" },
	{ "30729415811", "323233683197", "31071199" },
	{ "77874422", "32223899", "1" },
	{ "974507656412513757857315037382926980395082974811562770185617915360", "-1539496810360685510909469177732386446833404488164283", "1" },
	{ "2823618260546496405819033080103700734250203999069672146446", "18374686479688400895", "1" }
};

DO_BINOP_TEST(integer,cl_I,cl_I,gcd)

int test_gcd (void)
{
	int error = 0;
	error |= test_integer_gcd();
	return error;
}

int test_xgcd (void)
{
	int error = 0;
	{
		cl_I a = 77874422;
		cl_I b = 32223899;
		cl_I u;
		cl_I v;
		cl_I g = xgcd(a,b, &u,&v);
		ASSERT(g == 1);
		ASSERT(g == a*u+b*v);
		ASSERT(u == -9206830);
		ASSERT(v == 22249839);
	}
	{
		cl_I a = "560014183";
		cl_I b = 312839871;
		cl_I u;
		cl_I v;
		cl_I g = xgcd(a,b, &u,&v);
		ASSERT(g == 1);
		ASSERT(g == a*u+b*v);
		ASSERT(u == 77165803);
		ASSERT(v == -138134388);
	}
	{
		cl_I a = "#x80000000";
		cl_I b = "#x-C0000000";
		cl_I u;
		cl_I v;
		cl_I g = xgcd(a,b, &u,&v);
		ASSERT(g == (cl_I)"#x40000000");
		ASSERT(g == a*u+b*v);
		ASSERT(u == -1);
		ASSERT(v == -1);
	}
	{
		cl_I a = "974507656412513757857315037382926980395082974811562770185617915360";
		cl_I b = "-1539496810360685510909469177732386446833404488164283";
		cl_I u;
		cl_I v;
		cl_I g = xgcd(a,b, &u,&v);
		ASSERT(g == 1);
		ASSERT(g == a*u+b*v);
	}
	return error;
}
