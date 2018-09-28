// cl_UDS_recip().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/digitseq/cl_DS.h"


// Implementation.

namespace cln {

// Compute the reciprocal value of a digit sequence.
// Input: UDS a_MSDptr/a_len/.. of length a_len,
//        with  1/2*beta^a_len <= a < beta^a_len.
// Output: UDS b_MSDptr/b_len+2/.. of length b_len+1 (b_len>1), plus 1 more bit
//        in the last limb, such that
//        beta^b_len <= b <= 2*beta^b_len  and
//        | beta^(a_len+b_len)/a - b | < 1.
// If a_len > b_len, only the most significant b_len limbs + 3 bits of a
// are used.
  extern void cl_UDS_recip (const uintD* a_MSDptr, uintC a_len,
                            uintD* b_MSDptr, uintC b_len);
// Method:
// Using Newton/Heron iteration.
// Write x = a/beta^a_len and y = b/beta^b_len.
// So we start out with 1/2 <= x < 1 and search an y with 1 <= y <= 2
// and  | 1/x - y | < beta^(-b_len).
// For n = 1,2,...,b_len we compute approximations y with 1 <= yn <= 2
// and  | 1/x - yn | < beta^(-n). The first n limbs of x, plus the
// next 3 bits (of the (n+1)st limb) enter the computation of yn. Apart
// from that, yn remains valid for any x which shares the same n+1
// most significant limbs.
// Step n = 1:
//   Write  x = x1/beta + x2/beta^2 + xr with 0 <= xr < 1/(8*beta).
//   Divide (beta^2-beta*x1-x2) by x1, gives  beta^2-x1*beta-x2 = q*x1+r.
//   If this division overflows, i.e. q >= beta, then x1 = beta/2, x2 = 0,
//   and we just return y1 = 2.
//   Else set qd := ceiling(max(q*x2/beta - r, 0) / (x1+1)) and return
//   y1 = (beta+q-qd)/beta.
//   Rationale: Obviously 0 <= qd <= q and 0 <= qd <= 2. We have
//     beta^2 - beta*(beta+q-qd)*x
//     <= beta^2 - (beta+q-qd)*(x1 + x2/beta)
//     = beta^2 - beta*x1 - x2 - q*(x1 + x2/beta) + qd*(x1 + x2/beta)
//     = q*x1 + r - q*(x1 + x2/beta) + qd*(x1 + x2/beta)
//     = r - q*x2/beta + qd*(x1 + x2/beta)
//     if qd=0: <= r <= x1-1 < x1
//     if qd>0: < r - q*x2/beta + qd*(x1+1) <= x1
//     hence always < x1 <= beta*x, hence
//     1 - x*y1 <= x/beta, hence 1/x - y1 <= 1/beta.
//     And on the other hand
//     beta^2 - beta*(beta+q-qd)*x
//     = beta^2 - (beta+q-qd)*(x1 + x2/beta) - beta*(beta+q-qd)*xr
//     where the third term is
//       <= 2*beta^2*xr < beta/4 <= x1/2 <= beta*x/2.
//     Hence
//     beta^2 - beta*(beta+q-qd)*x >
//     > beta^2 - (beta+q-qd)*(x1 + x2/beta) - beta*x/2
//     = r - q*x2/beta + qd*(x1 + x2/beta) - beta*x/2
//     >= - qd - beta*x/2 > - beta*x, hence
//     1 - x*y1 >= -x/beta, hence 1/x - y1 >= -1/beta.
// Step n -> m with n < m <= 2*n:
//   Write x = xm + xr with 0 <= xr < 1/(8*beta^m).
//   Set ym' = 2*yn - xm*yn*yn,
//   ym = ym' rounded up to be a multiple of 1/(2*beta^m).
//   Rationale:
//     1/x - ym <= 1/x - ym' = 1/x - 2*yn + (x-xr)*yn*yn
//     <= 1/x - 2*yn + x*yn*yn = x * (1/x - yn)^2 < x*beta^(-2n)
//     < beta^(-2n) <= beta^(-m), and
//     1/x - ym' = 1/x - 2*yn + (x-xr)*yn*yn
//     > 1/x - 2*yn + x*yn*yn - 1/(2*beta^m)
//     = x * (1/x - yn)^2 - 1/(2*beta^m) >= - 1/(2*beta^m), hence
//     1/x - ym > 1/x - ym' - 1/(2*beta^m) >= -1/beta^m.
//   Since it is needed to compute ym as a multiple of 1/(2*beta^m),
//   not only as a multiple of 1/beta^m, we compute with zn = 2*yn.
//   The iteration now reads  zm = round_up(2*zn - xm*zn*zn/2). 
// Choice of n:
//   So that the computation is minimal, e.g. in the case b_len=10:
//   1 -> 2 -> 3 -> 5 -> 10 and not 1 -> 2 -> 4 -> 8 -> 10.
  void cl_UDS_recip (const uintD* a_MSDptr, uintC a_len,
                     uintD* b_MSDptr, uintC b_len)
    {
      var uintC y_len = b_len+1;
      var uintC x_len = (a_len <= b_len ? a_len+1 : y_len);
      var uintD* x_MSDptr;
      var uintD* y_MSDptr;
      var uintD* y2_MSDptr;
      var uintD* y3_MSDptr;
      CL_ALLOCA_STACK;
      num_stack_alloc(x_len,x_MSDptr=,);
      num_stack_alloc(y_len,y_MSDptr=,);
      num_stack_alloc(2*y_len,y2_MSDptr=,);
      num_stack_alloc(x_len+2*y_len,y3_MSDptr=,);
      // Prepare x/2 at x_MSDptr by shifting a right by 1 bit.
      if (a_len <= b_len)
        { mspref(x_MSDptr,a_len) =
            shiftrightcopy_loop_msp(a_MSDptr,x_MSDptr,a_len,1,0);
        }
        else
        { mspref(x_MSDptr,b_len) =
            shiftrightcopy_loop_msp(a_MSDptr,x_MSDptr,b_len,1,0)
            | ((mspref(a_MSDptr,b_len) & -bit(intDsize-3)) >> 1);
        }
      // Step n = 1.
      { var uintD x1 = mspref(a_MSDptr,0);
        var uintD x2 = (a_len > 1 ? (mspref(a_MSDptr,1) & -bit(intDsize-3)) : 0);
        if ((x1 == (uintD)bit(intDsize-1)) && (x2 == 0))
          { mspref(y_MSDptr,0) = 4; mspref(y_MSDptr,1) = 0; }
          else
          { var uintD q;
            var uintD r;
            var uintD chi;
            var uintD clo;
            #if HAVE_DD
              divuD((uintDD)(-highlowDD(x1,x2)),x1, q=,r=);
              var uintDD c = muluD(q,x2);
              chi = highD(c); clo = lowD(c);
            #else
              divuD((uintD)(-x1 - (x2>0 ? 1 : 0)),(uintD)(-x2),x1, q=,r=);
              muluD(q,x2,chi=,clo=);
            #endif
            if (clo > 0)
              chi++;
            // qd := ceiling(max(chi-r,0)/(x1+1))
            if (chi > r)
              { chi -= r;
                if (chi > x1)
                  { q--; }
                q--;
              }
            mspref(y_MSDptr,0) = 2 + (q>>(intDsize-1));
            mspref(y_MSDptr,1) = q<<1;
          }
      }
      // Other steps.
      var int k;
      integerlengthC(b_len-1,k=);
      // 2^(k-1) < b_len <= 2^k, so we need k steps.
      var uintC n = 1;
      for (; k>0; k--)
        { // n = ceiling(b_len/2^k) limbs of y have already been computed.
          var uintC m = ((b_len-1)>>(k-1))+1; // = ceiling(b_len/2^(k-1))
          // Compute zm := 2*zn - round_down(xm/2*zn*zn).
          cl_UDS_mul_square(y_MSDptr mspop (n+1),n+1,y2_MSDptr mspop 2*(n+1));
          var uintC xm_len = (m < x_len ? m+1 : x_len);
          cl_UDS_mul(x_MSDptr mspop xm_len,xm_len,
                     y2_MSDptr mspop 2*(n+1),2*n+1,
                     y3_MSDptr mspop (xm_len+2*n+1));
          // Round down by just taking the first m+1 limbs at y3_MSDptr.
          shift1left_loop_lsp(y_MSDptr mspop (n+1),n+1);
          clear_loop_msp(y_MSDptr mspop (n+1),m-n);
          subfrom_loop_lsp(y3_MSDptr mspop (m+1),y_MSDptr mspop (m+1),m+1);
          n = m;
        }
      // All n = b_len limbs of y have been computed. Divide by 2.
      mspref(b_MSDptr,b_len+1) = 
        shiftrightcopy_loop_msp(y_MSDptr,b_MSDptr,b_len+1,1,0);
    }
// Bit complexity (N := b_len): O(M(N)).

}  // namespace cln
