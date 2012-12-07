// futruncate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_LF futruncate (const cl_LF& x)
{
// Methode:
// x = 0.0 -> Ergebnis 0.0
// e<=0 -> Ergebnis 1.0 oder -1.0, je nach Vorzeichen von x.
// 1<=e<16n -> Greife die letzten (16n-e) Bits von x heraus.
//             Sind sie alle =0 -> Ergebnis x.
//             Sonst setze sie alle auf 0 und erhöhe dann die vorderen e Bits
//             um 1.
//             Kein Überlauf -> fertig.
//             Sonst (Ergebnis eine Zweierpotenz): Mantisse := .1000...000,
//               e:=e+1. (Test auf Überlauf wegen e<=16n überflüssig)
// e>=16n -> Ergebnis x.
#if 0
      var cl_signean sign;
      var sintE exp;
      var const uintD* mantMSDptr;
      var uintC mantlen;
      LF_decode(x, { return x; }, sign=,exp=,mantMSDptr=,mantlen=,);
      if (exp<=0) { return encode_LF1s(sign,mantlen); } // e<=0 -> Ergebnis +-1.0
      if ((uintE)exp >= intDsize*mantlen) // e>=16n -> x als Ergebnis
        { return x; }
        else
        // 0 < e < 16n
        { // Testen, ob alle hinteren 16n-e Bits =0 sind:
          var uintC count = floor((uintE)exp,intDsize); // zu kopierende Digits, < mantlen
          var uintC bitcount = ((uintE)exp) % intDsize; // zu kopierende Bits danach, >=0, <intDsize
          var uintD mask = minus_bitm(intDsize-bitcount); // Maske mit bitcount Bits
          var uintD* mantptr = mantMSDptr mspop count;
          if (   ((mspref(mantptr,0) & ~mask) ==0)
              && !test_loop_msp(mantptr mspop 1,mantlen-count-1)
             )
            { return x; }
          // neue NUDS erzeugen mit e Bits aus mant mit Increment
          // und 16n-e Nullbits:
          CL_ALLOCA_STACK;
          var uintD* MSDptr;
          num_stack_alloc(mantlen, MSDptr=,);
          { var uintD* ptr =
              copy_loop_msp(mantMSDptr,MSDptr,count); // count ganze Digits kopieren
            if ((mspref(ptr,0) = ((mspref(mantptr,0) & mask) - mask)) == 0) // dann bitcount Bits kopieren und incrementieren
              { if (!( inc_loop_lsp(ptr,count) ==0)) // evtl. weiterincrementieren
                  { mspref(MSDptr,0) = bit(intDsize-1); exp = exp+1; } // evtl. Exponenten erhöhen
              }
            clear_loop_msp(ptr mspop 1,mantlen-count-1); // Rest mit Nullen füllen
          }
          return encode_LF(sign,exp,MSDptr,mantlen);
        }
#else
      var uintC len = TheLfloat(x)->len;
      var uintE uexp = TheLfloat(x)->expo;
      if (uexp <= LF_exp_mid)
        { if (uexp == 0) { return x; } // x=0.0 -> Ergebnis 0.0
          return encode_LF1s(TheLfloat(x)->sign,len); // e<=0 -> Ergebnis +-1.0
        }
      var uintE exp = uexp - LF_exp_mid;
      if (exp >= intDsize*len) // e>=16n -> x als Ergebnis
        { return x; }
      // 0 < e < 16n
      // Testen, ob alle hinteren 16n-e Bits =0 sind:
      var uintC count = floor(exp,intDsize); // zu kopierende Digits, < mantlen
      var uintC bitcount = exp % intDsize; // zu kopierende Bits danach, >=0, <intDsize
      var uintD mask = minus_bitm(intDsize-bitcount); // Maske mit bitcount Bits
      {var const uintD* mantptr = LF_MSDptr(x) mspop count;
       if (   ((mspref(mantptr,0) & ~mask) ==0)
           && !test_loop_msp(mantptr mspop 1,len-count-1)
          )
         { return x; }
      }
      // Nein -> neues Long-Float produzieren:
      var Lfloat y = allocate_lfloat(len,uexp,TheLfloat(x)->sign); // neues Long-Float
      // y_mant := NUDS mit e Bits aus x_mant mit Increment und 16n-e Nullbits:
      {var const uintD* x_mantMSDptr = LF_MSDptr(x);
       var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
       var uintD* ptr =
         copy_loop_msp(x_mantMSDptr,y_mantMSDptr,count); // count ganze Digits kopieren
       if ((mspref(ptr,0) = ((mspref(x_mantMSDptr,count) & mask) - mask)) == 0) // dann bitcount Bits kopieren und incrementieren
         { if (!( inc_loop_lsp(ptr,count) ==0)) // evtl. weiterincrementieren
             { mspref(y_mantMSDptr,0) = bit(intDsize-1); (TheLfloat(y)->expo)++; } // evtl. Exponenten erhöhen
         }
       clear_loop_msp(ptr mspop 1,len-count-1); // Rest mit Nullen füllen
      }
      return y;
#endif
}

}  // namespace cln
