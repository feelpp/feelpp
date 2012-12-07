// cl_LF_I_div().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/lfloat/cl_LF_impl.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "float/cl_F.h"
#include "base/cl_N.h"

namespace cln {

const cl_LF cl_LF_I_div (const cl_LF& x, const cl_I& y)
{
// Method:
// If x = 0.0, return x.
// If y is longer than x, convert y to a float and divide.
// Else divide the mantissa of x by the absolute value of y, then round.
	if (TheLfloat(x)->expo == 0) {
		if (zerop(y))
			throw division_by_0_exception();
		else
			return x;
	}
	var cl_signean sign = -(cl_signean)minusp(y); // Vorzeichen von y
	var cl_I abs_y = (sign==0 ? y : -y);
	var uintC y_exp = integer_length(abs_y);
	var uintC len = TheLfloat(x)->len;
#ifndef CL_LF_PEDANTIC
	if (ceiling(y_exp,intDsize) > len)
		return x / cl_I_to_LF(y,len);
#endif
	// x länger als y, direkt dividieren.
	CL_ALLOCA_STACK;
	var const uintD* y_MSDptr;
	var uintC y_len;
	var const uintD* y_LSDptr;
	I_to_NDS_nocopy(abs_y, y_MSDptr=,y_len=,y_LSDptr=,false,); // NDS zu y bilden, y_len>0
	// y nicht zu einer NUDS normalisieren! (Damit ein Bit Spielraum ist.)
	// Zähler bilden: x * 2^(intDsize*y_len)
	var uintD* z_MSDptr;
	var uintC z_len;
	var uintD* z_LSDptr;
	z_len = len + y_len;
	num_stack_alloc(z_len, z_MSDptr=,z_LSDptr=);
	{ var uintD* ptr = copy_loop_msp(arrayMSDptr(TheLfloat(x)->data,len),z_MSDptr,len); // len Digits
	  clear_loop_msp(ptr,y_len); // und y_len Null-Digits
	}
	// Quotienten bilden:
	var DS q;
	var DS r;
	UDS_divide(z_MSDptr,z_len,z_LSDptr,
	           y_MSDptr,y_len,y_LSDptr,
	           &q, &r
	          );
	// q ist der Quotient,
	// q = floor(x/y*2^(intDsize*y_len)) >= 2^(intDsize*len),
	// q <= x/y*2^(intDsize*y_len) < 2^(1+intDsize+intDsize*len),
	// also mit len+1 oder len+2 Digits. r der Rest.
	var uintD* MSDptr = q.MSDptr;
	var uintL shiftcount;
	integerlengthD(mspref(MSDptr,0),shiftcount=);
	// Bei q.len = len+2 : shiftcount=1,
	// bei q.len = len+1 : 1<=shiftcount<=intDsize.
	var uintD carry_rechts;
	if (shiftcount==intDsize) {
		carry_rechts = mspref(MSDptr,len);
	} else {
		carry_rechts = shiftright_loop_msp(MSDptr,len+1,shiftcount%intDsize);
		if (q.len > len+1) {
			shiftcount += intDsize;
			carry_rechts |= (mspref(MSDptr,len+1)==0 ? 0 : 1);
		}
		msshrink(MSDptr);
	}
	// Quotient MSDptr/len/.. ist nun normalisiert: höchstes Bit =1.
	// exponent := exponent(x) - intDsize*y_len + shiftcount
	var uintE uexp = TheLfloat(x)->expo;
	var uintE dexp = intDsize*y_len - shiftcount; // >= 0 !
	if ((uexp < dexp) || ((uexp = uexp - dexp) < LF_exp_low)) {
		if (underflow_allowed())
			{ throw floating_point_underflow_exception(); }
		else
			{ return encode_LF0(len); }
	}
	// Runden:
	if ( ((sintD)carry_rechts >= 0) // herausgeschobenes Bit =0 -> abrunden
	     || ( (carry_rechts == (uintD)bit(intDsize-1)) // =1 und weitere Bits >0 oder Rest >0 -> aufrunden
	          && (r.len==0)
	          // round-to-even
		  && ((mspref(MSDptr,len-1) & bit(0)) ==0)
	   )    )
	  // abrunden
	  {}
	  else
	  // aufrunden
	  { if ( inc_loop_lsp(MSDptr mspop len,len) )
	      // Übertrag durchs Aufrunden
	      { mspref(MSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
	        // Exponenten incrementieren:
	        if (++uexp ==  LF_exp_high+1) { throw floating_point_overflow_exception(); }
	  }   }
	return encode_LFu(TheLfloat(x)->sign ^ sign, uexp, MSDptr, len);
}
// Bit complexity (N := max(length(x),length(y))): O(M(N)).

}  // namespace cln
