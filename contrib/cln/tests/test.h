#include <cln/io.h>
using namespace std;
using namespace cln;

#define ASSERT(expr)  \
  if (!(expr)) {					\
	std::cerr << "Assertion failed! File " << __FILE__ << ", line " << __LINE__ << endl;	\
	error = 1;					\
  }

#define ASSERT1(expr,a)  \
  if (!(expr)) {					\
	std::cerr << "Assertion failed! File " << __FILE__ << ", line " << __LINE__ << endl;	\
	std::cerr << #a" = " << a << endl;			\
	error = 1;					\
  }

#define ASSERT2(expr,a,b)  \
  if (!(expr)) {					\
	std::cerr << "Assertion failed! File " << __FILE__ << ", line " << __LINE__ << endl;	\
	std::cerr << #a" = " << a << endl;			\
	std::cerr << #b" = " << b << endl;			\
	error = 1;					\
  }

#define ASSERT3(expr,a,b,c)  \
  if (!(expr)) {					\
	std::cerr << "Assertion failed! File " << __FILE__ << ", line " << __LINE__ << endl;	\
	std::cerr << #a" = " << a << endl;			\
	std::cerr << #b" = " << b << endl;			\
	std::cerr << #c" = " << c << endl;			\
	error = 1;					\
  }

#define ASSERT4(expr,a,b,c,d)  \
  if (!(expr)) {					\
	std::cerr << "Assertion failed! File " << __FILE__ << ", line " << __LINE__ << endl;	\
	std::cerr << #a" = " << a << endl;			\
	std::cerr << #b" = " << b << endl;			\
	std::cerr << #c" = " << c << endl;			\
	std::cerr << #d" = " << d << endl;			\
	error = 1;					\
  }
