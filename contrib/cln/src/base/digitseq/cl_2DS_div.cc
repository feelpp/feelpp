// div2adic().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/digitseq/cl_2DS.h"


// Implementation.

#include "base/digit/cl_2D.h"
#include "base/digitseq/cl_DS.h"
#include "cln/exception.h"

namespace cln {

// Time for dividing a n word number by a n word number, this is the common
// case and therefore the important one:
// OS: Linux 2.2, intDsize==32,        OS: TRU64/4.0, intDsize==64,
// Machine: P-III/450MHz               Machine: EV5/300MHz:
//      n   standard  Newton             standard  Newton
//      30   0.00002   0.00006            0.00004   0.00020
//     100   0.00009   0.00045            0.00033   0.0015
//     300   0.00069   0.0028             0.0028    0.0085
//    1000   0.018     0.019              0.031     0.065
//    2000   0.028     0.057              0.12      0.20
//    3000   0.078     0.11  <-(~4500)    0.28      0.23  <-(~2700)
//   10000   1.09      0.48               3.14      1.13
//   30000  10.1       1.21              29.7       2.70
// Time for dividing a 2*n word number by a n word number:
// OS: Linux 2.2, intDsize==32,        OS: TRU64/4.0, intDsize==64,
// Machine: P-III/450MHz               Machine: EV5/300MHz:
//      n   standard  Newton             standard  Newton
//      30   0.00004   0.00019            0.00013   0.00067
//     100   0.00032   0.0014             0.0013    0.0046
//     300   0.0027    0.0084             0.011     0.025
//    1000   0.029     0.057              0.12      0.20
//    2000   0.16      0.18  <-(~2400)    0.50      0.46  <-(~1800)
//    3000   0.38      0.22               1.1       0.50
//   10000   4.5       1.05              13.0       2.48
//   30000  51.7       2.67             120.0       6.31
//          Newton faster for:         Newton faster for:
// 1.0*N / N  3300<N<3800, 4400<N        2700<N<3600, N<3800
// 1.1*N / N  3100<N<3700, 4100<N        2450<N<3300, N<3450
// 1.2*N / N  2850<N<3200, 3700<N        2250<N
// 1.3*N / N  2650<N<3000, 3450<N        2050<N
// 1.4*N / N  3400<N                     1850<N
// 1.5*N / N  3100<N                     1750<N
// 1.6*N / N  2850<N                     1650<N
// 1.7*N / N  2650<N                     1600<N
// 1.8*N / N  2550<N                     1500<N
// 1.9*N / N  2450<N                     1400<N
// 2.0*N / N  2400<N                     1350<N
//
// Break-even-point. When in doubt, prefer to choose the standard algorithm.
#if CL_USE_GMP
  static inline bool cl_recip_suitable (uintC m, uintC n) // n <= m
    { if (n < 2000)
        return false;
      else // when n >= 4400/(m/n)^2, i.e. (m/66)^2 > n
        { var uintC mq = floor(m,66);
          if ((mq >= bit(intCsize/2)) || (mq*mq > n))
            return true;
          else
            return false;
        }
    }
#else
// Use the old default values from CLN version <= 1.0.3 as a crude estimate.
// They came from timings on a i486 33 MHz running Linux:
// Divide N digits by N digits          Divide 2*N digits by N digits
//     N   standard  Newton                 N   standard  Newton
//     10    0.00015 0.00054                10    0.00023 0.00054
//     25    0.00065 0.00256                25    0.00116 0.00256
//     50    0.0024  0.0083                 50    0.0044  0.0082
//    100    0.0089  0.027                 100    0.0172  0.027
//    250    0.054   0.130                 250    0.107   0.130
//    500    0.22    0.42                  500    0.425   0.42  <-(~500)
//   1000    0.86    1.30                 1000    1.72    1.30
//   2500    5.6     4.1  <-(~2070)       2500   11.0     4.1
//   5000   22.3     9.4                  5000   44.7     9.3
//  10000   91.2    20.6                 10000  182      20.5
//
// 1.0*N / N : Newton for N >= 2070 or 1790 >= N >= 1460
// 1.1*N / N : Newton for N >= 1880 or 1790 >= N >= 1320
// 1.2*N / N : Newton for N >= 1250
// 1.3*N / N : Newton for N >= 1010
// 1.4*N / N : Newton for N >=  940
// 1.5*N / N : Newton for N >=  750
// 1.6*N / N : Newton for N >=  625
// 1.7*N / N : Newton for N >=  550
// 1.8*N / N : Newton for N >=  500
// 1.9*N / N : Newton for N >=  500
// 2.0*N / N : Newton for N >=  500
  static inline bool cl_recip_suitable (uintC m, uintC n) // n <= m
    { if (n < 500)
        return false;
      else // when n >= 2100/(m/n)^2, i.e. (m/46)^2 > n
        { var uintC mq = floor(m,46);
          if ((mq >= bit(intCsize/2)) || (mq*mq > n))
            return true;
          else
            return false;
        }
    }
#endif

void div2adic (uintC a_len, const uintD* a_LSDptr, uintC b_len, const uintD* b_LSDptr, uintD* dest_LSDptr)
{
  var uintC lendiff = a_len - b_len;
  if (cl_recip_suitable(a_len,b_len))
    { // Division using reciprocal (Newton-Hensel algorithm).
      CL_ALLOCA_STACK;
      // Bestimme Kehrwert c von b mod 2^(intDsize*b_len).
      var uintD* c_LSDptr;
      num_stack_alloc(b_len,,c_LSDptr=);
      recip2adic(b_len,b_LSDptr,c_LSDptr);
      // Bestimme q := a * c mod 2^(intDsize*b_len).
      var uintD* q_LSDptr;
      num_stack_alloc(2*b_len,,q_LSDptr=);
      cl_UDS_mul(a_LSDptr,b_len,c_LSDptr,b_len,q_LSDptr);
      // Zur Bestimmung des Restes wieder mit b multiplizieren:
      var uintD* p_LSDptr;
      num_stack_alloc(2*b_len,,p_LSDptr=);
      cl_UDS_mul(q_LSDptr,b_len,b_LSDptr,b_len,p_LSDptr);
      // Überprüfen, daß p == a mod 2^(intDsize*b_len):
      if (compare_loop_msp(a_LSDptr lspop b_len,p_LSDptr lspop b_len,b_len))
        throw runtime_exception();
      // Quotient q und "Rest" (a-b*q)/2^(intDsize*b_len) ablegen:
      copy_loop_lsp(q_LSDptr,dest_LSDptr,b_len);
      if (lendiff <= b_len)
        { sub_loop_lsp(a_LSDptr lspop b_len,p_LSDptr lspop b_len,dest_LSDptr lspop b_len,lendiff); }
        else
        { var uintD carry = sub_loop_lsp(a_LSDptr lspop b_len,p_LSDptr lspop b_len,dest_LSDptr lspop b_len,b_len);
          copy_loop_lsp(a_LSDptr lspop 2*b_len,dest_LSDptr lspop 2*b_len,lendiff-b_len);
          if (carry) { dec_loop_lsp(dest_LSDptr lspop 2*b_len,lendiff-b_len); }
        }
    }
    else
    { // Standard division.
      var uintD b0inv = div2adic(1,lspref(b_LSDptr,0)); // b'
      copy_loop_lsp(a_LSDptr,dest_LSDptr,a_len); // d := a
      do { var uintD digit = lspref(dest_LSDptr,0); // nächstes d[j]
           digit = mul2adic(b0inv,digit);
           // digit = nächstes c[j]
           if (a_len <= b_len)
             { mulusub_loop_lsp(digit,b_LSDptr,dest_LSDptr,a_len); } // d := d - b * c[j] * beta^j
             else
             // a_len > b_len, b wird als durch Nullen fortgesetzt gedacht.
             { var uintD carry = mulusub_loop_lsp(digit,b_LSDptr,dest_LSDptr,b_len);
               if (lspref(dest_LSDptr,b_len) >= carry)
                 { lspref(dest_LSDptr,b_len) -= carry; }
               else
                 { lspref(dest_LSDptr,b_len) -= carry;
                   dec_loop_lsp(dest_LSDptr lspop (b_len+1),a_len-(b_len+1));
             }   }
           // Nun ist lspref(dest_LSDptr,0) = 0.
           lspref(dest_LSDptr,0) = digit; // c[j] ablegen
           lsshrink(dest_LSDptr); a_len--; // nächstes j
         }
         until (a_len==lendiff);
    }
}
// Bit complexity (N = max(a_len,b_len)): O(M(N)).

}  // namespace cln
