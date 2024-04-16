// scale_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "float/cl_F.h"
#include "integer/cl_I.h"

namespace cln {

const cl_LF scale_float (const cl_LF& x, const cl_I& delta)
{
  // Methode:
  // delta=0 -> x als Ergebnis
  // x=0.0 -> x als Ergebnis
  // delta muß ein Integer betragsmäßig <= LF_exp_high-LF_exp_low sein.
  // Neues LF mit um delta vergrößertem Exponenten bilden.
      if (eq(delta,0)) { return x; } // delta=0 -> x als Ergebnis
      var uintE uexp = TheLfloat(x)->expo;
      if (uexp==0) { return x; }
      var uintE udelta;
      // |delta| muß <= LF_exp_high-LF_exp_low < 2^intEsize sein.
	if (fixnump(delta)) {
		// Fixnum
		var sintV sdelta = FN_to_V(delta);
		if (sdelta >= 0)
			{ udelta = sdelta; goto pos; }
		else
			{ udelta = sdelta; goto neg; }
	} else {
		// Bignum
		var cl_heap_bignum* bn = TheBignum(delta);
		if ((sintD)mspref(arrayMSDptr(bn->data,bn->length),0) >= 0) {
			// delta >= 0
			try {
				udelta = cl_I_to_UE(delta);
				goto pos;
			} catch (const runtime_exception&) {
				goto overflow;
			}
		} else {
			// delta < 0
			try {
				udelta = cl_I_to_E(delta);
				goto neg;
			} catch (const runtime_exception&) {
				goto underflow;
			}
		}
	}

      pos: // udelta = delta >=0
	if (   ((uexp = uexp+udelta) < udelta) // Exponent-Überlauf?
	    || (uexp > LF_exp_high) // oder Exponent zu groß?
	   )
	  overflow:
	  { throw floating_point_overflow_exception(); }
	goto ok;

      neg: // delta <0, udelta = 2^intEsize+delta
	if (   ((uexp = uexp+udelta) >= udelta) // oder Exponent-Unterlauf?
	    || (uexp < LF_exp_low) // oder Exponent zu klein?
	   )
	  underflow:
	  { if (underflow_allowed())
	    { throw floating_point_underflow_exception(); }
	    else
	    { return encode_LF0(TheLfloat(x)->len); } // Ergebnis 0.0
	  }
	goto ok;

      ok:
        var uintC len = TheLfloat(x)->len;
        return encode_LFu(TheLfloat(x)->sign,uexp,arrayMSDptr(TheLfloat(x)->data,len),len);
}

}  // namespace cln
