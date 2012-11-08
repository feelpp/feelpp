// sqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "float/cl_F.h"
#include "base/digitseq/cl_DS.h"
#include "cln/exception.h"

namespace cln {

const cl_LF sqrt (const cl_LF& x)
{
// Methode:
// x = 0.0 -> Ergebnis 0.0
// Ergebnis-Vorzeichen := positiv,
// Ergebnis-Exponent := ceiling(e/2),
// Ergebnis-Mantisse:
//   Erweitere die Mantisse (n Digits) um n+2 Nulldigits nach hinten.
//   Bei ungeradem e schiebe dies (oder nur die ersten n+1 Digits davon)
//     um 1 Bit nach rechts.
//   Bilde daraus die Ganzzahl-Wurzel, eine n+1-Digit-Zahl mit einer
//     fhrenden 1.
//   Runde das letzte Digit weg:
//     Bit 15 = 0 -> abrunden,
//     Bit 15 = 1, Rest =0 und Wurzel exakt -> round-to-even,
//     sonst aufrunden.
//   Bei rounding overflow Mantisse um 1 Bit nach rechts schieben
//     und Exponent incrementieren.
      var uintE uexp = TheLfloat(x)->expo;
      if (uexp==0) { return x; } // x=0.0 -> 0.0 als Ergebnis
      var uintC len = TheLfloat(x)->len;
      // Radikanden bilden:
      CL_ALLOCA_STACK;
      var uintD* r_MSDptr;
      var uintD* r_LSDptr;
      var uintC r_len = 2*len+2; // Länge des Radikanden
      num_stack_alloc(r_len, r_MSDptr=,r_LSDptr=);
      if ((uexp & bit(0)) == (LF_exp_mid & bit(0)))
        // Exponent gerade
        {var uintD* ptr =
           copy_loop_msp(arrayMSDptr(TheLfloat(x)->data,len),r_MSDptr,len); // n Digits kopieren
         clear_loop_msp(ptr,len+2); // n+2 Nulldigits anhï¿½gen
        }
        else
        // Exponent ungerade
        {var uintD carry_rechts = // n Digits kopieren und um 1 Bit rechts shiften
           shiftrightcopy_loop_msp(arrayMSDptr(TheLfloat(x)->data,len),r_MSDptr,len,1,0);
         var uintD* ptr = r_MSDptr mspop len;
         msprefnext(ptr) = carry_rechts; // ï¿½ertrag und
         clear_loop_msp(ptr,len+1); // n+1 Nulldigits anhï¿½gen
        }
      // Compute ((uexp - LF_exp_mid + 1) >> 1) + LF_exp_mid without risking
      // uintE overflow.
      uexp = ((uexp - ((LF_exp_mid - 1) & 1)) >> 1) - ((LF_exp_mid - 1) >> 1)
             + LF_exp_mid;
      // Ergebnis allozieren:
      var Lfloat y = allocate_lfloat(len,uexp,0);
      var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
      // Wurzel ziehen:
#ifndef CL_LF_PEDANTIC
      if (len > 2900) // This is about 15% faster
        { // Kehrwert der Wurzel errechnen:
          var uintD* s_MSDptr;
          var uintD* s_LSDptr;
          num_stack_alloc(len+2, s_MSDptr=,s_LSDptr=);
          cl_UDS_recipsqrt(r_MSDptr,r_len, s_MSDptr,len);
          // Mit dem Radikanden multiplizieren:
          var uintD* p_MSDptr;
          var uintD* p_LSDptr;
          num_stack_alloc(2*len+3, p_MSDptr=,p_LSDptr=);
          cl_UDS_mul(r_MSDptr mspop (len+1),len+1,s_LSDptr,len+2,p_LSDptr);
          // Ablegen und runden:
          copy_loop_msp(p_MSDptr mspop 1,y_mantMSDptr,len); // NUDS nach y kopieren
          if (mspref(p_MSDptr,0) == 0)
            { if ( ((sintD)mspref(p_MSDptr,len+1) >= 0) // nï¿½hstes Bit =0 -> abrunden
                   || ( ((mspref(p_MSDptr,len+1) & ((uintD)bit(intDsize-1)-1)) ==0) // =1 und weitere Bits >0 -> aufrunden
                        && !test_loop_msp(p_MSDptr mspop (len+2),len+1)
                        // round-to-even (etwas witzlos, da eh alles ungenau ist)
                        && ((mspref(p_MSDptr,len) & bit(0)) ==0)
                 )    )
                // abrunden
                {}
                else
                // aufrunden
                { if ( inc_loop_lsp(y_mantMSDptr mspop len,len) )
                    // ï¿½ertrag durchs Aufrunden
                    { mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
                      (TheLfloat(y)->expo)++; // Exponenten incrementieren
                }   }
            }
            else
            // ï¿½ertrag durch Rundungsfehler
            { if (test_loop_msp(y_mantMSDptr,len)) throw runtime_exception();
              mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
              (TheLfloat(y)->expo)++; // Exponenten incrementieren
            }
          return y;
        }
#endif
      var DS w;
      var bool exactp;
      UDS_sqrt(r_MSDptr,r_len,r_LSDptr, &w, exactp=);
      // w ist die Ganzzahl-Wurzel, eine n+1-Digit-Zahl.
      copy_loop_msp(w.MSDptr,y_mantMSDptr,len); // NUDS nach y kopieren
      // Runden:
      if ( ((sintD)lspref(w.LSDptr,0) >= 0) // nï¿½hstes Bit =0 -> abrunden
           || ( ((lspref(w.LSDptr,0) & ((uintD)bit(intDsize-1)-1)) ==0) // =1 und weitere Bits >0 oder Rest >0 -> aufrunden
                && exactp
                // round-to-even
                && ((lspref(w.LSDptr,1) & bit(0)) ==0)
         )    )
        // abrunden
        {}
        else
        // aufrunden
        { if ( inc_loop_lsp(y_mantMSDptr mspop len,len) )
            // ï¿½ertrag durchs Aufrunden
            { mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
              (TheLfloat(y)->expo)++; // Exponenten incrementieren
        }   }
      return y;
}
// Bit complexity (N := length(x)): O(M(N)).

}  // namespace cln
