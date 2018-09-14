// cl_LF_I_mul().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "float/lfloat/cl_LF_impl.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "float/cl_F.h"

namespace cln {

const cl_R cl_LF_I_mul (const cl_LF& x, const cl_I& y)
{
// Method:
// If y=0, return 0.
// If x=0.0, return x.
// If y is longer than x, convert y to a float and multiply.
// Else multiply the mantissa of x with the absolute value of y, then round.
	if (eq(y,0)) { return 0; }
	if (TheLfloat(x)->expo == 0) { return x; }
	var cl_signean sign = -(cl_signean)minusp(y); // Vorzeichen von y
	var cl_I abs_y = (sign==0 ? y : -y);
	var uintC y_exp = integer_length(abs_y);
	var uintC len = TheLfloat(x)->len;
#ifndef CL_LF_PEDANTIC
	if (ceiling(y_exp,intDsize) > len)
		return x * cl_I_to_LF(y,len);
#endif
	// x länger als y, direkt multiplizieren.
	CL_ALLOCA_STACK;
	var const uintD* y_MSDptr;
	var uintC y_len;
	var const uintD* y_LSDptr;
	I_to_NDS_nocopy(abs_y, y_MSDptr=,y_len=,y_LSDptr=,false,); // NDS zu y bilden, y_len>0
	if (mspref(y_MSDptr,0)==0) y_len--; // NUDS zu y bilden, y_len>0
	// Multiplizieren.
	var uintD* prodMSDptr;
	var uintC prodlen;
	UDS_UDS_mul_UDS(len,arrayLSDptr(TheLfloat(x)->data,len),
	                y_len,y_LSDptr,
	                prodMSDptr=,prodlen=,);
	// x fing mit 0 Nullbits an, y mit maximal intDsize-1 Nullbits,
	// daher fängt das Produkt mit maximal intDsize Nullbits an.
	var uintL shiftcount;
	if (mspref(prodMSDptr,0)==0) {
		shiftcount = intDsize;
		msshrink(prodMSDptr); prodlen--;
	} else {
		integerlengthD(mspref(prodMSDptr,0), shiftcount = intDsize -);
		if (shiftcount > 0)
			shiftleft_loop_lsp(prodMSDptr mspop (len+1),len+1,shiftcount,0);
	}
	// Produkt ist nun normalisiert: höchstes Bit =1.
	// exponent := exponent(x) + intDsize*y_len - shiftcount
	var uintE uexp = TheLfloat(x)->expo;
	var uintE iexp = intDsize*y_len - shiftcount; // >= 0 !
	uexp = uexp + iexp;
	if ((uexp < iexp) || (uexp > LF_exp_high))
		throw floating_point_overflow_exception();
	// Runden:
	var uintD* midptr = prodMSDptr mspop len;
	var uintC restlen = prodlen - len;
	if ( (restlen==0)
	     || ((sintD)mspref(midptr,0) >= 0) // nächstes Bit =0 -> abrunden
	     || ( ((mspref(midptr,0) & ((uintD)bit(intDsize-1)-1)) ==0) // Bit =1, weitere Bits >0 -> aufrunden
	          && !test_loop_msp(midptr mspop 1,restlen-1)
	          // round-to-even
	          && ((lspref(midptr,0) & bit(0)) ==0)
	   )    )
	  // abrunden
	  {}
	  else
	  // aufrunden
	  { if ( inc_loop_lsp(midptr,len) )
	      // Übertrag durchs Aufrunden
	      { mspref(prodMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
	        if (++uexp == LF_exp_high+1) { throw floating_point_overflow_exception(); }
	  }   }
	return encode_LFu(TheLfloat(x)->sign ^ sign, uexp, prodMSDptr, len);
}
// Bit complexity (N = max(length(x),length(y))): O(M(N)).

}  // namespace cln
