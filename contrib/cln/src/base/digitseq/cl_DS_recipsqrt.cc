// cl_UDS_recipsqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/digitseq/cl_DS.h"


// Implementation.

#include "base/cl_low.h"
#include "cln/exception.h"

namespace cln {

// Compute the reciprocal square root of a digit sequence.
// Input: UDS a_MSDptr/a_len/.. of length a_len,
//        with 1/4 <= a < 1.
//        [i.e. 1/4*beta^a_len <= a < beta^a_len]
// Output: UDS b_MSDptr/b_len+2/.. of length b_len+1 (b_len>1), plus 1 more bit
//         in the last limb) such that
//         1 <= b <= 2  [i.e. beta^b_len <= b <= 2*beta^b_len]
//         and  | 1/sqrt(a) - b | < 1/2*beta^(-b_len).
// If a_len > b_len, only the most significant b_len+1 limbs of a are used.
  extern void cl_UDS_recipsqrt (const uintD* a_MSDptr, uintC a_len,
                                uintD* b_MSDptr, uintC b_len);
// Method:
// Using Newton iteration for computation of x^-1/2.
// The Newton iteration for f(y) = x-1/y^2 reads:
//   y --> y - (x-1/y^2)/(2/y^3) = y + y*(1-x*y^2)/2 =: g(y).
// We have  T^3-3*T+2 = (T-1)^2*(T+2), hence
//   1/sqrt(x) - g(y) = 1/(2*sqrt(x)) * (sqrt(x)*y-1)^2 * (sqrt(x)*y+2).
// Hence g(y) <= 1/sqrt(x).
// If we choose 0 < y_0 <= 1/sqrt(x), then set y_(n+1) := g(y_n), we will
// always have 0 < y_n <= 1/sqrt(x).
// Since
//   1/sqrt(x) - g(y) = sqrt(x)*(sqrt(x)*y+2)/2 * (1/sqrt(x) - y)^2,
// which is >= 0 and < 3/2 * (1/sqrt(x) - y)^2, we have a quadratically
// convergent iteration.
// For n = 1,2,...,b_len we compute approximations y with 1 <= yn <= 2
// and  | 1/sqrt(x) - yn | < 1/2*beta^(-n).
// Step n=1:
//   Compute the isqrt of the leading two digits of x, yields one digit.
//   Compute its reciprocal, then do one iteration as below (n=0 -> m=1).
// Step n -> m with n < m <= 2*n:
//   Write x = xm + xr with 0 <= xr < beta^-(m+1).
//   Set ym' = yn + (yn*(1-xm*yn*yn))/2, round down to a multiple ym
//   of beta^-(m+1).
//   (Actually, compute yn*yn, round up to a multiple of beta^-(m+1),   [1]
//    multiply with xm,        round up to a multiple of beta^-(m+1),   [2]
//    subtract from 1,         no rounding needed,                      [2]
//    multiply with yn,        round down to a multiple of beta^-(m+1), [5]
//    divide by 2,             round down to a multiple of beta^-(m+1), [3]
//    add to yn,               no rounding needed.  [Max rounding error: ^])
//   The exact value ym' (no rounding) would satisfy
//     0 <= 1/sqrt(xm) - ym' < 3/2 * (1/sqrt(xm) - yn)^2
//                           < 3/8 * beta^(-2*n)          by hypothesis,
//                           <= 3/8 * beta^-m.
//   The rounding errors all go into the same direction, so
//     0 <= ym' - ym < 3 * beta^-(m+1) < 1/4 * beta^-m.
//   Combine both inequalities:
//     0 <= 1/sqrt(xm) - ym < 1/2 * beta^-m.
//   Neglecting xr can introduce a small error in the opposite direction:
//     0 <= 1/sqrt(xm) - 1/sqrt(x) = (sqrt(x) - sqrt(xm))/(sqrt(x)*sqrt(xm))
//        = xr / (sqrt(x)*sqrt(xm)*(sqrt(x)+sqrt(xm)))
//        <= 4*xr < 4*beta^-(m+1) < 1/2*beta^-m.
//   Combine both inequalities:
//     | 1/sqrt(x) - ym | < 1/2 * beta^-m.
//   (Actually, choosing the opposite rounding direction wouldn't hurt either.)
// Choice of n:
//   So that the computation is minimal, e.g. in the case b_len=10:
//   1 -> 2 -> 3 -> 5 -> 10 and not 1 -> 2 -> 4 -> 8 -> 10.
  void cl_UDS_recipsqrt (const uintD* a_MSDptr, uintC a_len,
                         uintD* b_MSDptr, uintC b_len)
    {
	var uintC y_len = b_len+2;
	var uintC x_len = (a_len <= b_len ? a_len : b_len+1);
	var const uintD* const x_MSDptr = a_MSDptr;
	var uintD* y_MSDptr;
	var uintD* y2_MSDptr;
	var uintD* y3_MSDptr;
	var uintD* y4_MSDptr;
	CL_ALLOCA_STACK;
	num_stack_alloc(y_len,y_MSDptr=,);
	num_stack_alloc(2*y_len,y2_MSDptr=,);
	num_stack_alloc(2*y_len,y3_MSDptr=,);
	num_stack_alloc(2*y_len,y4_MSDptr=,);
	// Step n = 1.
	{ var uintD x1 = mspref(x_MSDptr,0);
	  var uintD x2 = (a_len > 1 ? mspref(x_MSDptr,1) : 0);
	  var uintD y0;
	  var uintD y1;
	  var bool sqrtp;
	  isqrtD(x1,x2, y1=,sqrtp=);
	  // 2^31 <= y1 < 2^32.
	  y0 = 1;
	  if (!sqrtp) // want to compute 1/sqrt(x) rounded down
		if (++y1 == 0)
			goto step1_done; // 1/1.0000 = 1.0000
	  // Set y0|y1 := 2^(2*intDsize)/y1
	  //            = 2^intDsize + (2^(2*intDsize)-2^intDsize*y1)/y1.
	  if ((uintD)(-y1) >= y1) {
		y0 = 2; y1 = 0;
	  } else {
		#if HAVE_DD
		divuD(highlowDD_0((uintD)(-y1)),y1, y1=,);
		#else
		divuD((uintD)(-y1),0,y1, y1=,);
		#endif
	  }
	step1_done:
	  mspref(y_MSDptr,0) = y0;
	  mspref(y_MSDptr,1) = y1;
	}
	// Other steps.
	var int k;
	integerlengthC(b_len-1,k=);
	// 2^(k-1) < b_len <= 2^k, so we need k steps, plus one
	// one more step at the beginning (because step 1 was not complete).
	var uintC n = 0;
	for (; k>=0; k--)
	  { var uintC m = ((b_len-1)>>k)+1; // = ceiling(b_len/2^k)
	    // Compute ym := yn + (yn*(1-xm*yn*yn))/2, rounded.
	    // Storage: at y_MSDptr: (1 + n+1) limbs, yn.
	    //          at y2_MSDptr: (2 + 2*n+2) limbs, yn^2.
	    //          at y3_MSDptr: (1 + m+1) limbs, xm*yn*yn, 1-xm*yn*yn.
	    //          at y4_MSDptr: (2-n + m+n+2) limbs, yn*(1-xm*yn*yn).
	    clear_loop_msp(y_MSDptr mspop (n+2),m-n);
	    cl_UDS_mul_square(y_MSDptr mspop (n+2),n+2,
	                      y2_MSDptr mspop 2*(n+2));
	    var uintC xm_len = (m < x_len ? m+1 : x_len);
	    var uintC y2_len = m+2; // = (m+1 <= 2*n+2 ? m+2 : 2*n+3);
	    cl_UDS_mul(x_MSDptr mspop xm_len,xm_len,
	               y2_MSDptr mspop (y2_len+1),y2_len,
	               y3_MSDptr mspop (xm_len+y2_len));
	    if (mspref(y3_MSDptr,0)==0)
	      // xm*yn*yn < 1
	      { neg_loop_lsp(y3_MSDptr mspop (m+2),m+2);
	        mspref(y3_MSDptr,0) += 1;
	        if (test_loop_msp(y3_MSDptr,n)) throw runtime_exception(); // check 0 <= y3 < beta^-(n-1)
	        cl_UDS_mul(y_MSDptr mspop (n+2),n+2,
	                   y3_MSDptr mspop (m+2),m+2-n,
	                   y4_MSDptr mspop (m+4));
	        shift1right_loop_msp(y4_MSDptr,m+3-n,0);
	        if (addto_loop_lsp(y4_MSDptr mspop (m+3-n),y_MSDptr mspop (m+2),m+3-n))
	          if ((n<1) || inc_loop_lsp(y_MSDptr mspop (n-1),n-1)) throw runtime_exception();
	      }
	      else
	      // xm*yn*yn >= 1 (this can happen since xm >= xn)
	      { mspref(y3_MSDptr,0) -= 1;
	        if (test_loop_msp(y3_MSDptr,n)) throw runtime_exception(); // check 0 >= y3 > -beta^-(n-1)
	        cl_UDS_mul(y_MSDptr mspop (n+2),n+2,
	                   y3_MSDptr mspop (m+2),m+2-n,
	                   y4_MSDptr mspop (m+4));
	        shift1right_loop_msp(y4_MSDptr,m+3-n,0);
	        if (subfrom_loop_lsp(y4_MSDptr mspop (m+3-n),y_MSDptr mspop (m+2),m+3-n))
	          if ((n<1) || dec_loop_lsp(y_MSDptr mspop (n-1),n-1)) throw runtime_exception();
	      }
	    n = m;
	    // n = ceiling(b_len/2^k) limbs of y have now been computed.
	  }
	copy_loop_msp(y_MSDptr,b_MSDptr,b_len+2);
}
// Bit complexity (N := b_len): O(M(N)).

}  // namespace cln
