// binary operator /

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"
#include "float/cl_F.h"
#include "base/cl_N.h"

namespace cln {

const cl_LF operator/ (const cl_LF& x1, const cl_LF& x2)
{
// Methode:
// x2 = 0.0 -> Error
// x1 = 0.0 -> Ergebnis 0.0
// Sonst:
// Ergebnis-Vorzeichen = xor der beiden Vorzeichen von x1 und x2
// Ergebnis-Exponent = Differenz der beiden Exponenten von x1 und x2
// Ergebnis-Mantisse = Mantisse mant1 / Mantisse mant2, gerundet.
//   mant1/mant2 > 1/2, mant1/mant2 < 2;
//   nach Rundung mant1/mant2 >=1/2, <=2*mant1<2.
//   Bei mant1/mant2 >=1 brauche 16n-1 Nachkommabits,
//   bei mant1/mant2 <1 brauche 16n Nachkommabits.
//   Fürs Runden: brauche ein Rundungsbit (Rest gibt an, ob exakt).
//   Brauche daher insgesamt 16n+1 Nachkommabits von mant1/mant2.
//   Dividiere daher (als Unsigned Integers)
//     2^16(n+1)*(2^16n*m0) durch (2^16n*m1).
//   Falls der Quotient >=2^16(n+1) ist, schiebe ihn um 1 Bit nach rechts,
//     erhöhe den Exponenten um 1 und runde das letzte Digit weg.
//   Falls der Quotient <2^16(n+1) ist, runde das letzte Digit weg. Bei rounding
//     overflow schiebe um 1 Bit nach rechts und erhöhe den Exponenten um 1.
      var uintC len1 = TheLfloat(x1)->len;
      var uintC len2 = TheLfloat(x2)->len;
      var uintC len = (len1 < len2 ? len1 : len2); // min. Länge n von x1 und x2
      var uintE uexp2 = TheLfloat(x2)->expo;
      if (uexp2==0) { throw division_by_0_exception(); } // x2=0.0 -> Error
      var uintE uexp1 = TheLfloat(x1)->expo;
      if (uexp1==0) // x1=0.0 -> Ergebnis 0.0
        { if (len < len1) return shorten(x1,len); else return x1; }
      // Exponenten subtrahieren:
      // (uexp1-LF_exp_mid) - (uexp2-LF_exp_mid) = (uexp1-uexp2+LF_exp_mid)-LF_exp_mid
      if (uexp1 >= uexp2)
        { uexp1 = uexp1 - uexp2; // kein Carry
          if (uexp1 > LF_exp_high-LF_exp_mid) { throw floating_point_overflow_exception(); }
          uexp1 = uexp1 + LF_exp_mid;
        }
        else
        { uexp1 = uexp1 - uexp2; // Carry
          if (uexp1 < (uintE)(LF_exp_low-1-LF_exp_mid))
            { if (underflow_allowed())
                { throw floating_point_underflow_exception(); }
                else
                { return encode_LF0(len); } // Ergebnis 0.0
            }
          uexp1 = uexp1 + LF_exp_mid;
        }
      // Nun ist LF_exp_low-1 <= uexp1 <= LF_exp_high.
      // neues Long-Float allozieren:
      var Lfloat y = allocate_lfloat(len,uexp1,
                                     TheLfloat(x1)->sign ^ TheLfloat(x2)->sign // Vorzeichen kombinieren
                                    );
      // Nenner bilden:
      var uintC n_len;
      n_len = len2;
#ifndef CL_LF_PEDANTIC
      if (n_len > len) { n_len = len+1; }
#endif
      // Zähler bilden:
      CL_ALLOCA_STACK;
      var uintD* z_MSDptr;
      var uintC z_len;
      var uintD* z_LSDptr;
      z_len = n_len + len + 1;
      num_stack_alloc(z_len, z_MSDptr=,z_LSDptr=);
      if (z_len > len1)
        { var uintD* ptr =
            copy_loop_msp(arrayMSDptr(TheLfloat(x1)->data,len1),z_MSDptr,len1); // n Digits kopieren
          clear_loop_msp(ptr,z_len-len1); // und n+1 Null-Digits
        }
        else
        { copy_loop_msp(arrayMSDptr(TheLfloat(x1)->data,len1),z_MSDptr,z_len); }
      // Quotienten bilden: 2n+1-Digit-Zahl durch n-Digit-Zahl dividieren
      {var DS q;
       var DS r;
       {var uintD* x2_mantMSDptr = arrayMSDptr(TheLfloat(x2)->data,len2);
        UDS_divide(z_MSDptr,z_len,z_LSDptr,
                   x2_mantMSDptr,n_len,x2_mantMSDptr mspop n_len,
                   &q, &r
                  );
       }
       // q ist der Quotient mit n+1 oder n+2 Digits, r der Rest.
       if (q.len > len+1)
         // Quotient hat n+2 Digits -> um 1 Bit nach rechts schieben:
         { var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
           var uintD carry_rechts =
             shiftrightcopy_loop_msp(q.MSDptr mspop 1,y_mantMSDptr,len,1,
                                     /* carry links = mspref(q.MSDptr,0) = 1 */ 1 );
           // Exponenten incrementieren:
           if (++(TheLfloat(y)->expo) == LF_exp_high+1) { throw floating_point_overflow_exception(); }
           // Runden:
           if ( (carry_rechts == 0) // herausgeschobenes Bit =0 -> abrunden
                || ( (lspref(q.LSDptr,0)==0) // =1 und weitere Bits >0 oder Rest >0 -> aufrunden
                     && (r.len==0)
                     // round-to-even
                     && ((lspref(q.LSDptr,1) & bit(1)) ==0)
              )    )
             // abrunden
             {}
             else
             // aufrunden
             { inc_loop_lsp(y_mantMSDptr mspop len,len); }
         }
         else
         // Quotient hat n+1 Digits -> nur kopieren:
         { var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
           copy_loop_msp(q.MSDptr,y_mantMSDptr,len);
           // Runden:
           if ( ((sintD)lspref(q.LSDptr,0) >= 0) // nächstes Bit =0 -> abrunden
                || ( ((lspref(q.LSDptr,0) & ((uintD)bit(intDsize-1)-1)) ==0) // =1 und weitere Bits >0 oder Rest >0 -> aufrunden
                     && (r.len==0)
                     // round-to-even
                     && ((lspref(q.LSDptr,1) & bit(0)) ==0)
              )    )
             // abrunden
             {}
             else
             // aufrunden
             { if ( inc_loop_lsp(y_mantMSDptr mspop len,len) )
                 // Übertrag durchs Aufrunden
                 { mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
                   // Exponenten incrementieren:
                   if (++(TheLfloat(y)->expo) == LF_exp_high+1) { throw floating_point_overflow_exception(); }
             }   }
         }
      }
      // LF_exp_low <= exp <= LF_exp_high sicherstellen:
      if (TheLfloat(y)->expo == LF_exp_low-1)
        { if (underflow_allowed())
            { throw floating_point_underflow_exception(); }
            else
            { return encode_LF0(len); } // Ergebnis 0.0
        }
      return y;
}
// Bit complexity (N := max(length(x1),length(x2))): O(M(N)).

}  // namespace cln
