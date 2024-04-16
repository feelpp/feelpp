// fround().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_LF fround (const cl_LF& x)
{
// Methode:
// x = 0.0 oder e<0 -> Ergebnis 0.0
// 0<=e<16n -> letzte (16n-e) Bits der Mantisse wegrunden,
//             Exponent und Vorzeichen beibehalten.
// e>=16n -> Ergebnis x
#if 0
      var cl_signean sign;
      var sintE exp;
      var const uintD* mantMSDptr;
      var uintC mantlen;
      LF_decode(x, { return x; }, sign=,exp=,mantMSDptr=,mantlen=,);
      if (exp<0) { return encode_LF0(mantlen); } // e<0 -> Ergebnis 0.0
      if ((uintE)exp >= intDsize*mantlen) // e>=16n -> x als Ergebnis
        { return x; }
        else
        // 0 <= e < 16n
        { // alle hinteren 16n-e Bits wegrunden:
          var uintC count = floor((uintE)exp,intDsize); // zu kopierende Digits, < mantlen
          var uintC bitcount = ((uintE)exp) % intDsize; // zu kopierende Bits danach, >=0, <intDsize
          var uintD mask = minus_bit(intDsize-bitcount-1); // Maske mit bitcount+1 Bits
          var const uintD* mantptr = mantMSDptr mspop count;
          if ((mspref(mantptr,0) & -mask) ==0) goto ab; // Bit 16n-e-1 =0 -> abrunden
          if (!((mspref(mantptr,0) & ~mask) ==0)) goto auf; // Bit 16n-e-1 =1 und Bits 16n-e-2..0 >0 -> aufrunden
          if (test_loop_msp(mantptr mspop 1,mantlen-count-1)) goto auf;
          // round-to-even, je nach Bit 16n-e :
          if (bitcount>0)
            { if ((mspref(mantptr,0) & (-2*mask)) ==0) goto ab; else goto auf; }
            elif (count>0)
              { if ((lspref(mantptr,0) & bit(0)) ==0) goto ab; else goto auf; }
              else
              // bitcount=0, count=0, also exp=0: Abrunden von +-0.5 zu 0.0
              { return encode_LF0(mantlen); }
          ab: // abrunden
          { CL_ALLOCA_STACK;
            var uintD* MSDptr;
            num_stack_alloc(mantlen, MSDptr=,);
           {var uintD* ptr =
              copy_loop_msp(mantMSDptr,MSDptr,count); // count ganze Digits kopieren
            msprefnext(ptr) = mspref(mantptr,0) & mask; // dann bitcount Bits kopieren
            clear_loop_msp(ptr,mantlen-count-1); // Rest mit Nullen füllen
            return encode_LF(sign,exp,MSDptr,mantlen);
          }}
          auf: // aufrunden
          { CL_ALLOCA_STACK;
            var uintD* MSDptr;
            num_stack_alloc(mantlen, MSDptr=,);
           {var uintD* ptr =
              copy_loop_msp(mantMSDptr,MSDptr,count); // count ganze Digits kopieren
            if ((mspref(ptr,0) = ((mspref(mantptr,0) & mask) - mask)) == 0) // dann bitcount Bits kopieren und incrementieren
              { if (!( inc_loop_lsp(ptr,count) ==0)) // evtl. weiterincrementieren
                  { mspref(MSDptr,0) = bit(intDsize-1); exp = exp+1; } // evtl. Exponenten erhöhen
              }
            clear_loop_msp(ptr mspop 1,mantlen-count-1); // Rest mit Nullen füllen
            return encode_LF(sign,exp,MSDptr,mantlen);
          }}
        }
#else
      var uintC len = TheLfloat(x)->len;
      var uintE uexp = TheLfloat(x)->expo;
      if (uexp < LF_exp_mid)
        { if (uexp == 0) { return x; } // x=0.0 -> Ergebnis 0.0
          return encode_LF0(len); // e<0 -> Ergebnis 0.0
        }
      var uintE exp = uexp - LF_exp_mid;
      if (exp >= intDsize*len) // e>=16n -> x als Ergebnis
        { return x; }
      // 0 <= e < 16n
      // alle hinteren 16n-e Bits wegrunden:
      var uintC count = floor(exp,intDsize); // zu kopierende Digits, < mantlen
      var uintC bitcount = exp % intDsize; // zu kopierende Bits danach, >=0, <intDsize
      var uintD mask = minus_bit(intDsize-bitcount-1); // Maske mit bitcount+1 Bits
      {var const uintD* mantptr = LF_MSDptr(x) mspop count;
       if ((mspref(mantptr,0) & -mask) ==0) goto ab; // Bit 16n-e-1 =0 -> abrunden
       if (!((mspref(mantptr,0) & ~mask) ==0)) goto auf; // Bit 16n-e-1 =1 und Bits 16n-e-2..0 >0 -> aufrunden
       if (test_loop_msp(mantptr mspop 1,len-count-1)) goto auf;
       // round-to-even, je nach Bit 16n-e :
       if (bitcount>0)
         { if ((mspref(mantptr,0) & (-2*mask)) ==0) goto ab; else goto auf; }
         elif (count>0)
           { if ((lspref(mantptr,0) & bit(0)) ==0) goto ab; else goto auf; }
           else
           // bitcount=0, count=0, also exp=0: Abrunden von +-0.5 zu 0.0
           { return encode_LF0(len); }
      }
      ab: // abrunden
        {var Lfloat y = allocate_lfloat(len,uexp,TheLfloat(x)->sign); // neues Long-Float
         // y_mant := NUDS mit e Bits aus x_mant und 16n-e Nullbits:
         {var const uintD* x_mantMSDptr = LF_MSDptr(x);
          var uintD* ptr =
            copy_loop_msp(x_mantMSDptr,arrayMSDptr(TheLfloat(y)->data,len),count); // count ganze Digits kopieren
          msprefnext(ptr) = mspref(x_mantMSDptr,count) & mask; // dann bitcount Bits kopieren
          clear_loop_msp(ptr,len-count-1); // Rest mit Nullen füllen
         }
         return y;
        }
      auf: // aufrunden
        {var Lfloat y = allocate_lfloat(len,uexp,TheLfloat(x)->sign); // neues Long-Float
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
        }
#endif
}

}  // namespace cln
