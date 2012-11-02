// square().

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

const cl_LF square (const cl_LF& x)
{
// Methode: wie operator*(x,x).
      var uintC len = TheLfloat(x)->len;
      var uintE uexp = TheLfloat(x)->expo;
      if (uexp==0) // x=0.0 -> Ergebnis 0.0
        { return x; }
      // Exponenten addieren:
      // (uexp-LF_exp_mid) + (uexp-LF_exp_mid) = (2*uexp-LF_exp_mid)-LF_exp_mid
      if ((sintE)uexp >= 0)
        // kein Carry
        { uexp = 2*uexp;
          if (uexp < LF_exp_mid+LF_exp_low)
            { if (underflow_allowed())
                { throw floating_point_underflow_exception(); }
                else
                { return encode_LF0(len); } // Ergebnis 0.0
        }   }
        else
        // Carry
        { uexp = 2*uexp;
          if (uexp > (uintE)(LF_exp_mid+LF_exp_high+1)) { throw floating_point_overflow_exception(); }
        }
      uexp = uexp - LF_exp_mid;
      // Nun ist LF_exp_low <= uexp <= LF_exp_high+1.
      // neues Long-Float allozieren:
      var Lfloat y = allocate_lfloat(len,uexp,0);
      // Produkt bilden:
      var const uintD* x_LSDptr = arrayLSDptr(TheLfloat(x)->data,len);
      var uintD* MSDptr;
      var uintD* LSDptr;
      CL_ALLOCA_STACK;
      num_stack_alloc(2*len,MSDptr=,LSDptr=);
      cl_UDS_mul_square(x_LSDptr,len,LSDptr);
      {var uintD* midptr = MSDptr mspop len; // Pointer in die Mitte der 2*len Digits
       if ((sintD)mspref(MSDptr,0) >= 0) // führendes Bit abtesten
         { // erste n+1 Digits um 1 Bit nach links schieben:
           shift1left_loop_lsp(midptr mspop 1,len+1);
           // Exponenten decrementieren:
           if ((TheLfloat(y)->expo)-- == LF_exp_low-1)
             { if (underflow_allowed())
                 { throw floating_point_underflow_exception(); }
                 else
                 { return encode_LF0(len); } // Ergebnis 0.0
             }
         }
       // erste Hälfte des Mantissenprodukts übertragen:
       {var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
        var uintD* y_mantLSDptr = copy_loop_msp(MSDptr,y_mantMSDptr,len);
        // Runden:
        if ( ((sintD)mspref(midptr,0) >= 0) // nächstes Bit =0 -> abrunden
             || ( ((mspref(midptr,0) & ((uintD)bit(intDsize-1)-1)) ==0) // Bit =1, weitere Bits >0 -> aufrunden
                  && !test_loop_msp(midptr mspop 1,len-1)
                  // round-to-even
                  && ((lspref(midptr,0) & bit(0)) ==0)
           )    )
          // abrunden
          {}
          else
          // aufrunden
          { if ( inc_loop_lsp(y_mantLSDptr,len) )
              { // Übertrag durchs Aufrunden (kann nur auftreten,
                // wenn vorhin um 1 Bit nach links geschoben wurde)
                mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
                (TheLfloat(y)->expo)++; // Exponent wieder zurück-erhöhen
          }   }
        // LF_exp_low <= exp <= LF_exp_high sicherstellen:
        if (TheLfloat(y)->expo == LF_exp_high+1) { throw floating_point_overflow_exception(); }
      }}
      return y;
}
// Bit complexity (N = length(x)): O(M(N)).

}  // namespace cln
