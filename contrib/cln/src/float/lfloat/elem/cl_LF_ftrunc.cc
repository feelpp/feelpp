// ftruncate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_LF ftruncate (const cl_LF& x)
{
// Methode:
// x = 0.0 oder e<=0 -> Ergebnis 0.0
// 1<=e<=16n -> letzte (16n-e) Bits der Mantisse auf 0 setzen,
//              Exponent und Vorzeichen beibehalten
// e>=16n -> Ergebnis x
#if 0
      var cl_signean sign;
      var sintE exp;
      var const uintD* mantMSDptr;
      var uintC mantlen;
      LF_decode(x, { return x; }, sign=,exp=,mantMSDptr=,mantlen=,);
      if (exp<=0) { return encode_LF0(mantlen); } // e<=0 -> Ergebnis 0.0
      if ((uintE)exp >= intDsize*mantlen) // e>=16n -> x als Ergebnis
        { return x; }
        else
        // 0 < e < 16n
        // neue NUDS erzeugen mit e Bits aus mant und 16n-e Nullbits:
        { CL_ALLOCA_STACK;
          var uintD* MSDptr;
          num_stack_alloc(mantlen, MSDptr=,);
          { var uintC count = floor((uintE)exp,intDsize); // zu kopierende Digits, < mantlen
            var uintC bitcount = ((uintE)exp) % intDsize; // zu kopierende Bits danach, >=0, <intDsize
            var uintD* ptr =
              copy_loop_msp(mantMSDptr,MSDptr,count); // count ganze Digits kopieren
            msprefnext(ptr) = mspref(mantMSDptr,count) & minus_bitm(intDsize-bitcount); // dann bitcount Bits kopieren
            clear_loop_msp(ptr,mantlen-count-1); // Rest mit Nullen füllen
          }
          return encode_LF(sign,exp,MSDptr,mantlen);
        }
#else
      var uintC len = TheLfloat(x)->len;
      var uintE uexp = TheLfloat(x)->expo;
      if (uexp <= LF_exp_mid)
        { if (uexp == 0) { return x; } // x=0.0 -> Ergebnis 0.0
          return encode_LF0(len); // e<=0 -> Ergebnis 0.0
        }
      var uintE exp = uexp - LF_exp_mid;
      if (exp >= intDsize*len) // e>=16n -> x als Ergebnis
        { return x; }
      // 0 < e < 16n
      var Lfloat y = allocate_lfloat(len,uexp,TheLfloat(x)->sign); // neues Long-Float
      // y_mant := NUDS mit e Bits aus x_mant und 16n-e Nullbits:
      {var uintC count = floor(exp,intDsize); // zu kopierende Digits, < mantlen
       var uintC bitcount = exp % intDsize; // zu kopierende Bits danach, >=0, <intDsize
       var const uintD* x_mantMSDptr = LF_MSDptr(x);
       var uintD* ptr =
         copy_loop_msp(x_mantMSDptr,arrayMSDptr(TheLfloat(y)->data,len),count); // count ganze Digits kopieren
       msprefnext(ptr) = mspref(x_mantMSDptr,count) & minus_bitm(intDsize-bitcount); // dann bitcount Bits kopieren
       clear_loop_msp(ptr,len-count-1); // Rest mit Nullen füllen
      }
      return y;
#endif
}

}  // namespace cln
