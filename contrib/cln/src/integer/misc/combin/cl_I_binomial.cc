// binomial().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/misc/combin/cl_I_combin.h"

namespace cln {

const cl_I binomial (uintL n, uintL k)
{
	// Method:
	// Ensure 0 <= k <= n: If n<k, return 0.
	// Ensure 0 <= k <= n/2: If k>n/2, replace k by n-k.
	// Compute the product (n-k+1)...n and divide by k!.
	// When computing the product m < i <= n, split into the odd part
	// and the even part. The even part is 2^(ord2(n!)-ord2(m!)),
	// and recall that ord2(n!) = n - logcount(n). The odd part is
	// computed recursively: It is the product of the odd part of
	// m/2 < i <= n/2, times the product of m < 2*i+1 <= n. The
	// recursion goes until floor(m/2^j) = floor(n/2^j) or floor(n/2^j) < 2.
	if (n < k)
		return 0;
	// Now 0 <= k <= n.
	if (2*k > n)
		k = n-k;
	// Now 0 <= k <= n/2.
	var cl_I prod = 1;
	var uintL m = n-k;
	var uintL j = 0;
	{
		var uintL a = m;
		var uintL b = n;
		while (a < b && b > 1) {
			a = floor(a,2); b = floor(b,2);
			j++;
		}
	}
	while (j > 0) {
		j--;
		var uintL a = m>>j;
		var uintL b = n>>j;
		// Compute product(a < i <= b and i odd, i)
		a = floor(a-1,2);
		b = floor(b-1,2);
		// = product(a < i <= b, 2i+1)
		if (a < b)
			prod = prod * cl_I_prod_ungerade(a,b);
	}
	prod = prod << (k + logcount(m) - logcount(n));
	return exquopos(prod,factorial(k));
}

}  // namespace cln
