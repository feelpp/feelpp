// binary operator *

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"
#include "float/cl_F.h"

namespace cln {

const cl_LF operator* (const cl_LF& x1, const cl_LF& x2)
{
// Methode:
// Falls x1=0.0 oder x2=0.0 -> Ergebnis 0.0
// Sonst: Ergebnis-Vorzeichen = VZ von x1 xor VZ von x2.
//        Ergebnis-Exponent = Summe der Exponenten von x1 und x2.
//        Produkt der Mantissen bilden (2n Digits).
//        Falls das fhrende Bit =0 ist: Mantissenprodukt um 1 Bit nach links
//          schieben (die vorderen n+1 Digits gengen)
//          und Exponent decrementieren.
//        Runden auf n Digits liefert die Ergebnis-Mantisse.
      var uintC len1 = TheLfloat(x1)->len;
      var uintC len2 = TheLfloat(x2)->len;
      var uintC len = (len1 < len2 ? len1 : len2); // min. L�ge n von x1 und x2
      var uintE uexp1 = TheLfloat(x1)->expo;
      if (uexp1==0) // x1=0.0 -> Ergebnis 0.0
        { if (len < len1) return shorten(x1,len); else return x1; }
      var uintE uexp2 = TheLfloat(x2)->expo;
      if (uexp2==0) // x2=0.0 -> Ergebnis 0.0
        { if (len < len2) return shorten(x2,len); else return x2; }
      // Exponenten addieren:
      // (uexp1-LF_exp_mid) + (uexp2-LF_exp_mid) = (uexp1+uexp2-LF_exp_mid)-LF_exp_mid
      uexp1 = uexp1 + uexp2;
      if (uexp1 >= uexp2)
        // kein Carry
        { if (uexp1 < LF_exp_mid+LF_exp_low)
            { if (underflow_allowed())
                { throw floating_point_underflow_exception(); }
                else
                { return encode_LF0(len); } // Ergebnis 0.0
        }   }
        else
        // Carry
        { if (uexp1 > (uintE)(LF_exp_mid+LF_exp_high+1)) { throw floating_point_overflow_exception(); } }
      uexp1 = uexp1 - LF_exp_mid;
      // Nun ist LF_exp_low <= uexp1 <= LF_exp_high+1.
      // neues Long-Float allozieren:
      var Lfloat y = allocate_lfloat(len,uexp1,
                                     TheLfloat(x1)->sign ^ TheLfloat(x2)->sign // Vorzeichen kombinieren
                                    );
      // Produkt bilden:
      var const uintD* x1_LSDptr = arrayLSDptr(TheLfloat(x1)->data,len1);
      var const uintD* x2_LSDptr = arrayLSDptr(TheLfloat(x2)->data,len2);
#ifndef CL_LF_PEDANTIC
      if (len1 > len2)
        { x1_LSDptr = x1_LSDptr lspop (len1-(len2+1)); len1 = len2+1; }
      else if (len1 < len2)
        { x2_LSDptr = x2_LSDptr lspop (len2-(len1+1)); len2 = len1+1; }
#endif
      var uintD* MSDptr;
      CL_ALLOCA_STACK;
      UDS_UDS_mul_UDS(len1,x1_LSDptr,
                      len2,x2_LSDptr,
                      MSDptr=,,);
      {var uintD* midptr = MSDptr mspop len; // Pointer in die Mitte der len1+len2 Digits
       if ((sintD)mspref(MSDptr,0) >= 0) // fhrendes Bit abtesten
         { // erste n+1 Digits um 1 Bit nach links schieben:
           shift1left_loop_lsp(midptr mspop 1,len+1);
           // Exponenten decrementieren:
           if (--(TheLfloat(y)->expo) == LF_exp_low-1)
             { if (underflow_allowed())
                 { throw floating_point_underflow_exception(); }
                 else
                 { return encode_LF0(len); } // Ergebnis 0.0
             }
         }
       // erste H�fte des Mantissenprodukts bertragen:
       {var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
        var uintD* y_mantLSDptr = copy_loop_msp(MSDptr,y_mantMSDptr,len);
        // Runden:
        if ( ((sintD)mspref(midptr,0) >= 0) // n�hstes Bit =0 -> abrunden
             || ( ((mspref(midptr,0) & ((uintD)bit(intDsize-1)-1)) ==0) // Bit =1, weitere Bits >0 -> aufrunden
                  && !test_loop_msp(midptr mspop 1,len1+len2-len-1)
                  // round-to-even
                  && ((lspref(midptr,0) & bit(0)) ==0)
           )    )
          // abrunden
          {}
          else
          // aufrunden
          { if ( inc_loop_lsp(y_mantLSDptr,len) )
              { // �ertrag durchs Aufrunden (kann nur auftreten,
                // wenn vorhin um 1 Bit nach links geschoben wurde)
                mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
                (TheLfloat(y)->expo)++; // Exponent wieder zurck-erhhen
          }   }
        // LF_exp_low <= exp <= LF_exp_high sicherstellen:
        if (TheLfloat(y)->expo == LF_exp_high+1) { throw floating_point_overflow_exception(); }
      }}
      return y;
}
// Bit complexity (N = max(length(x1),length(x2))): O(M(N)).

}  // namespace cln
