// float_approx().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "float/cl_F.h"

namespace cln {

float float_approx (const cl_I& x)
{
// Method: same as cl_I_to_FF().
      if (eq(x,0)) { return 0.0; }
      var cl_signean sign = -(cl_signean)minusp(x); // Vorzeichen
      var cl_I abs_x = (sign==0 ? x : -x);
      var uintC exp = integer_length(abs_x); // (integer-length x)
      // NDS zu |x|>0 bilden:
      var const uintD* MSDptr;
      var uintC len;
      I_to_NDS_nocopy(abs_x, MSDptr=,len=,,false,);
      // MSDptr/len/LSDptr ist die NDS zu x, len>0.
      // Führende Digits holen: Brauche FF_mant_len+1 Bits, dazu intDsize
      // Bits (die NDS kann mit bis zu intDsize Nullbits anfangen).
      // Dann werden diese Bits um (exp mod intDsize) nach rechts geschoben.
      var uintD msd = msprefnext(MSDptr); // erstes Digit
      #if (intDsize==64)
      var uintD msdd = 0; // weiteres Digit
      if (--len == 0) goto ok;
      msdd = msprefnext(MSDptr);
      #else // (intDsize<=32)
      var uint32 msdd = 0; // weitere min(len-1,32/intDsize) Digits
      #define NEXT_DIGIT(i)  \
        { if (--len == 0) goto ok;                                   \
          msdd |= (uint32)msprefnext(MSDptr) << (32-(i+1)*intDsize); \
        }
      DOCONSTTIMES(32/intDsize,NEXT_DIGIT);
      #undef NEXT_DIGIT
      #endif
      --len; ok:
      #if (intDsize==64)
      // Die NDS besteht aus msd, msdd, und len weiteren Digits.
      // Das höchste in 2^intDsize*msd+msdd gesetzte Bit ist Bit Nummer
      // intDsize-1 + (exp mod intDsize).
      var uintL shiftcount = exp % intDsize;
      var uint64 mant = // führende 64 Bits
        (shiftcount==0
         ? msdd
         : ((msd << (64-shiftcount)) | (msdd >> shiftcount))
        );
      // Das höchste in mant gesetzte Bit ist Bit Nummer 63.
      if ( ((mant & bit(62-FF_mant_len)) ==0) // Bit 39 =0 -> abrunden
           || ( ((mant & (bit(62-FF_mant_len)-1)) ==0) // Bit 39 =1 und Bits 38..0 =0
                && ((msdd & (bit(shiftcount)-1)) ==0) // und weitere Bits aus msdd =0
                && (!test_loop_msp(MSDptr,len)) // und alle weiteren Digits =0
                // round-to-even, je nach Bit 40 :
                && ((mant & bit(63-FF_mant_len)) ==0)
         )    )
        // abrunden
        { mant = mant >> (63-FF_mant_len); }
        else
        // aufrunden
        { mant = mant >> (63-FF_mant_len);
          mant += 1;
          if (mant >= bit(FF_mant_len+1)) // rounding overflow?
            { mant = mant>>1; exp = exp+1; }
        }
      #else
      // Die NDS besteht aus msd, msdd, und len weiteren Digits.
      // Das höchste in 2^32*msd+msdd gesetzte Bit ist Bit Nummer
      // 31 + (exp mod intDsize).
      var uintL shiftcount = exp % intDsize;
      var uint32 mant = // führende 32 Bits
        (shiftcount==0
         ? msdd
         : (((uint32)msd << (32-shiftcount)) | (msdd >> shiftcount))
        );
      // Das höchste in mant gesetzte Bit ist Bit Nummer 31.
      if ( ((mant & bit(30-FF_mant_len)) ==0) // Bit 7 =0 -> abrunden
           || ( ((mant & (bit(30-FF_mant_len)-1)) ==0) // Bit 7 =1 und Bits 6..0 =0
                && ((msdd & (bit(shiftcount)-1)) ==0) // und weitere Bits aus msdd =0
                && (!test_loop_msp(MSDptr,len)) // und alle weiteren Digits =0
                // round-to-even, je nach Bit 8 :
                && ((mant & bit(31-FF_mant_len)) ==0)
         )    )
        // abrunden
        { mant = mant >> (31-FF_mant_len); }
        else
        // aufrunden
        { mant = mant >> (31-FF_mant_len);
          mant += 1;
          if (mant >= bit(FF_mant_len+1)) // rounding overflow?
            { mant = mant>>1; exp = exp+1; }
        }
      #endif
      union { ffloat eksplicit; float machine_float; } u;
      if ((sintL)exp > (sintL)(FF_exp_high-FF_exp_mid))
        { u.eksplicit = make_FF_word(sign,bit(FF_exp_len)-1,0); } // Infinity
      else
        { u.eksplicit = make_FF_word(sign,(sintL)exp+FF_exp_mid,mant); }
      return u.machine_float;
}

}  // namespace cln
