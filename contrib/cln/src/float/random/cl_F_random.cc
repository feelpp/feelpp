// random_F().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "base/random/cl_random_impl.h"
#include "base/digitseq/cl_DS.h"
#include "integer/cl_I.h"

namespace cln {

const cl_F random_F (random_state& randomstate, const cl_F& n)
{
	var uintC d = float_digits(n); // d = (float-digits n) > 0
	// Bilde neue UDS mit d Zufallsbits:
	CL_ALLOCA_STACK;
	var uintC len = ceiling(d,intDsize);
	var uintD* MSDptr;
	num_stack_alloc_1(len,MSDptr=,);
	random_UDS(randomstate,MSDptr,len); // len (>0) Zufallsdigits
	// von intDsize*ceiling(d/intDsize) auf d Bits herunterschneiden:
	{ var uintL dr = d % intDsize; if (dr>0) { mspref(MSDptr,0) &= (bit(dr)-1); } }
	// in Integer umwandeln:
	var cl_I mant = UDS_to_I(MSDptr,len);
	// Bilde  Zufalls-Float zwischen 0 und 1
	//        = (scale-float (float Zufalls-Integer,d_Bits n) (- d)) :
	var cl_F result = scale_float(cl_float(mant,n),-(sintC)d) * n;
	// result ist ein Zufalls-Float >=0, <=n.
	if (result == n)
		// falls (durch Rundung) result=n, durch 0 ersetzen:
		{ result = cl_float(0,result); }
	return result;
}

}  // namespace cln
