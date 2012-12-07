// cl_I_to_LF().

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

const cl_LF cl_I_to_LF (const cl_I& x, uintC len)
{
// Methode:
// x=0 -> Ergebnis 0.0
// Merke Vorzeichen von x.
// x:=(abs x)
// Exponent:=(integer-length x)
// Mantisse enthalte die höchstwertigen 16n Bits des Integers x (wobei die
//   führenden 16-(e mod 16) Nullbits zu streichen sind).
// Runde die weiteren Bits weg:
//   Kommen keine mehr -> abrunden,
//   nächstes Bit = 0 -> abrunden,
//   nächstes Bit = 1 und Rest =0 -> round-to-even,
//   nächstes Bit = 1 und Rest >0 -> aufrunden.
// Bei Aufrundung: rounding overflow -> Mantisse um 1 Bit nach rechts schieben
//   und Exponent incrementieren.
      if (eq(x,0)) { return encode_LF0(len); } // x=0 -> Ergebnis 0.0
      var cl_signean sign = -(cl_signean)minusp(x); // Vorzeichen von x
      var cl_I abs_x = (sign==0 ? x : -x);
      var uintC exp = integer_length(abs_x); // (integer-length x) < intDsize*2^intCsize
      // Teste, ob exp <= LF_exp_high-LF_exp_mid :
      if (   (log2_intDsize+intCsize < 32)
          && ((uintE)(intDsize*bitc(intCsize)-1) <= (uintE)(LF_exp_high-LF_exp_mid))
         )
        {} // garantiert exp <= intDsize*2^intCsize-1 <= LF_exp_high-LF_exp_mid
        else
        { if (!(exp <= (uintE)(LF_exp_high-LF_exp_mid))) { throw floating_point_overflow_exception(); } }
      // Long-Float bauen:
      var Lfloat y = allocate_lfloat(len,exp+LF_exp_mid,sign);
      var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
      var const uintD* x_MSDptr;
      var uintC x_len;
      I_to_NDS_nocopy(abs_x, x_MSDptr=,x_len=,,false,); // NDS zu x bilden, x_len>0
      // x_MSDptr/x_len/.. um (exp mod 16) Bits nach rechts shiften und in
      // y einfüllen (genauer: nur maximal len Digits davon):
      {var uintL shiftcount = exp % intDsize;
       // Die NDS fängt mit intDsize-shiftcount Nullbits an, dann kommt eine 1.
       if (x_len > len)
         { x_len -= 1+len;
           if (shiftcount>0)
             { var uintD carry_rechts =
                 shiftrightcopy_loop_msp(x_MSDptr mspop 1,y_mantMSDptr,len,shiftcount,mspref(x_MSDptr,0));
               // Mantisse ist gefüllt. Runden:
               if ( ((sintD)carry_rechts >= 0) // nächstes Bit =0 -> abrunden
                    || ( ((carry_rechts & ((uintD)bit(intDsize-1)-1)) ==0) // =1, Rest >0 -> aufrunden
                         && !test_loop_msp(x_MSDptr mspop 1 mspop len,x_len)
                         // round-to-even
                         && ((mspref(y_mantMSDptr,len-1) & bit(0)) ==0)
                  )    )
                 goto ab; // aufrunden
                 else
                 goto auf; // aufrunden
             }
             else
             { copy_loop_msp(x_MSDptr mspop 1,y_mantMSDptr,len);
               // Mantisse ist gefüllt. Runden:
               var const uintD* ptr = x_MSDptr mspop 1 mspop len;
               if ( (x_len==0) // keine Bits mehr -> abrunden
                    || ((sintD)mspref(ptr,0) >= 0) // nächstes Bit =0 -> abrunden
                    || ( ((mspref(ptr,0) & ((uintD)bit(intDsize-1)-1)) ==0) // =1, Rest >0 -> aufrunden
                         && !test_loop_msp(ptr mspop 1,x_len-1)
                         // round-to-even
                         && ((lspref(ptr,0) & bit(0)) ==0)
                  )    )
                 goto ab; // aufrunden
                 else
                 goto auf; // aufrunden
             }
           auf: // aufrunden
             if ( inc_loop_lsp(y_mantMSDptr mspop len,len) )
               // Übertrag durchs Aufrunden
               { mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
                 // Exponenten incrementieren:
                 if (   (log2_intDsize+intCsize < 32)
                     && ((uintE)(intDsize*bitc(intCsize)-1) < (uintE)(LF_exp_high-LF_exp_mid))
                    )
                   // garantiert exp < intDsize*2^intCsize-1 <= LF_exp_high-LF_exp_mid
                   { (TheLfloat(y)->expo)++; } // jetzt exp <= LF_exp_high-LF_exp_mid
                   else
                   { if (++(TheLfloat(y)->expo) == LF_exp_high+1) { throw floating_point_overflow_exception(); } }
               }
           ab: // abrunden
             ;
         }
         else // x_len <= len
         { var uintD carry_rechts;
           len -= x_len;
           x_len -= 1;
           if (shiftcount>0)
             { carry_rechts = shiftrightcopy_loop_msp(x_MSDptr mspop 1,y_mantMSDptr,x_len,shiftcount,mspref(x_MSDptr,0)); }
             else
             { copy_loop_msp(x_MSDptr mspop 1,y_mantMSDptr,x_len); carry_rechts = 0; }
          {var uintD* y_ptr = y_mantMSDptr mspop x_len;
           msprefnext(y_ptr) = carry_rechts; // Carry als nächstes Digit
           clear_loop_msp(y_ptr,len); // dann len-x_len Nulldigits
         }}
      }
      return y;
}

}  // namespace cln
