// scale_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "float/cl_F.h"

namespace cln {

const cl_LF scale_float (const cl_LF& x, sintC delta)
{
  // Methode:
  // delta=0 -> x als Ergebnis
  // x=0.0 -> x als Ergebnis
  // delta muß ein Integer betragsmäßig <= LF_exp_high-LF_exp_low sein.
  // Neues LF mit um delta vergrößertem Exponenten bilden.
      if (delta == 0) { return x; } // delta=0 -> x als Ergebnis
      var uintE uexp = TheLfloat(x)->expo;
      if (uexp==0) { return x; }
      var uintE udelta = delta;
      if (delta >= 0) {
        // udelta = delta >=0
	if (   ((uexp = uexp+udelta) < udelta) // Exponent-Überlauf?
	    || (uexp > LF_exp_high) // oder Exponent zu groß?
	   )
	  { throw floating_point_overflow_exception(); }
      } else {
        // delta <0, udelta = 2^intEsize+delta
	if (   ((uintE)(-(uexp = uexp+udelta)) <= (uintE)(-udelta)) // oder Exponent-Unterlauf?
	    || (uexp < LF_exp_low) // oder Exponent zu klein?
	   )
	  { if (underflow_allowed())
	    { throw floating_point_underflow_exception(); }
	    else
	    { return encode_LF0(TheLfloat(x)->len); } // Ergebnis 0.0
	  }
      }
      var uintC len = TheLfloat(x)->len;
      return encode_LFu(TheLfloat(x)->sign,uexp,arrayMSDptr(TheLfloat(x)->data,len),len);
}

}  // namespace cln
