// double_approx().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "float/cl_F.h"

namespace cln {

double double_approx (const cl_I& x)
{
// Method: same as cl_I_to_DF().
      if (eq(x,0)) { return 0.0; }
      var cl_signean sign = -(cl_signean)minusp(x); // Vorzeichen
      var cl_I abs_x = (sign==0 ? x : -x);
      var uintC exp = integer_length(abs_x); // (integer-length x)
      // NDS zu |x|>0 bilden:
      var const uintD* MSDptr;
      var uintC len;
      I_to_NDS_nocopy(abs_x, MSDptr=,len=,,false,);
      // MSDptr/len/LSDptr ist die NDS zu x, len>0.
      // Führende Digits holen: Brauche DF_mant_len+1 Bits, dazu intDsize
      // Bits (die NDS kann mit bis zu intDsize Nullbits anfangen).
      // Dann werden diese Bits um (exp mod intDsize) nach rechts geschoben.
      var uintD msd = msprefnext(MSDptr); // erstes Digit
      #if (cl_word_size==64)
      var uint64 msdd = 0; // weitere min(len-1,64/intDsize) Digits
      #define NEXT_DIGIT(i)  \
        { if (--len == 0) goto ok;                                   \
          msdd |= (uint64)msprefnext(MSDptr) << (64-(i+1)*intDsize); \
        }
      DOCONSTTIMES(64/intDsize,NEXT_DIGIT);
      #undef NEXT_DIGIT
      #else
      var uint32 msdd = 0; // weitere min(len-1,32/intDsize) Digits
      var uint32 msddf = 0; // weitere maximal 32/intDsize Digits
      #define NEXT_DIGIT(i)  \
        { if (--len == 0) goto ok;                                   \
          msdd |= (uint32)msprefnext(MSDptr) << (32-(i+1)*intDsize); \
        }
      DOCONSTTIMES(32/intDsize,NEXT_DIGIT);
      #undef NEXT_DIGIT
      #define NEXT_DIGIT(i)  \
        { if (--len == 0) goto ok;                                    \
          msddf |= (uint32)msprefnext(MSDptr) << (32-(i+1)*intDsize); \
        }
      DOCONSTTIMES(32/intDsize,NEXT_DIGIT);
      #undef NEXT_DIGIT
      #endif
      --len; ok:
      #if (cl_word_size==64)
      // Die NDS besteht aus msd, msdd und len weiteren Digits.
      // Das höchste in 2^64*msd+msdd gesetzte Bit ist Bit Nummer
      // 63 + (exp mod intDsize).
      var uintL shiftcount = exp % intDsize;
      var uint64 mant = // führende 64 Bits
        (shiftcount==0
          ? msdd
          : (((uint64)msd << (64-shiftcount)) | (msdd >> shiftcount))
        );
      // Das höchste in mant gesetzte Bit ist Bit Nummer 63.
      if ( ((mant & bit(62-DF_mant_len)) ==0) // Bit 10 =0 -> abrunden
           || ( ((mant & (bit(62-DF_mant_len)-1)) ==0) // Bit 10 =1 und Bits 9..0 =0
                && ((msdd & (bit(shiftcount)-1)) ==0) // und weitere Bits aus msddf =0
                && (!test_loop_msp(MSDptr,len)) // und alle weiteren Digits =0
                // round-to-even, je nach Bit 11 :
                && ((mant & bit(63-DF_mant_len)) ==0)
         )    )
        // abrunden
        { mant = mant >> (63-DF_mant_len); }
        else
        // aufrunden
        { mant = mant >> (63-DF_mant_len);
          mant += 1;
          if (mant >= bit(DF_mant_len+1)) // rounding overflow?
            { mant = mant>>1; exp = exp+1; }
        }
      union { dfloat eksplicit; double machine_double; } u;
      if ((sintL)exp > (sintL)(DF_exp_high-DF_exp_mid))
        { u.eksplicit =
              ((sint64)sign & bit(63))
            | ((uint64)(bit(DF_exp_len)-1) << DF_mant_len); // Infinity
        }
      else
        { u.eksplicit =
              ((sint64)sign & bit(63))                         /* Vorzeichen */
            | ((uint64)((sintL)exp+DF_exp_mid) << DF_mant_len) /* Exponent   */
            | ((uint64)mant & (bit(DF_mant_len)-1));           /* Mantisse   */
        }
      return u.machine_double;
      #else
      // Die NDS besteht aus msd, msdd, msddf und len weiteren Digits.
      // Das höchste in 2^64*msd+2^32*msdd+msddf gesetzte Bit ist Bit Nummer
      // 63 + (exp mod intDsize).
      var uintL shiftcount = exp % intDsize;
      var uint32 manthi; // führende 32 Bits
      var uint32 mantlo; // nächste 32 Bits
      if (shiftcount==0)
        { manthi = msdd; mantlo = msddf; }
        else
        { manthi = ((uint32)msd << (32-shiftcount)) | (msdd >> shiftcount);
          mantlo = (msdd << (32-shiftcount)) | (msddf >> shiftcount);
        }
      // Das höchste in mant gesetzte Bit ist Bit Nummer 63.
      if ( ((mantlo & bit(62-DF_mant_len)) ==0) // Bit 10 =0 -> abrunden
           || ( ((mantlo & (bit(62-DF_mant_len)-1)) ==0) // Bit 10 =1 und Bits 9..0 =0
                && ((msddf & (bit(shiftcount)-1)) ==0) // und weitere Bits aus msddf =0
                && (!test_loop_msp(MSDptr,len)) // und alle weiteren Digits =0
                // round-to-even, je nach Bit 11 :
                && ((mantlo & bit(63-DF_mant_len)) ==0)
         )    )
        // abrunden
        { mantlo = (mantlo >> (63-DF_mant_len)) | (manthi << (DF_mant_len-32+1));
          manthi = manthi >> (63-DF_mant_len);
        }
        else
        // aufrunden
        { mantlo = (mantlo >> (63-DF_mant_len)) | (manthi << (DF_mant_len-32+1));
          manthi = manthi >> (63-DF_mant_len);
          mantlo += 1;
          if (mantlo==0)
            { manthi += 1;
              if (manthi >= bit(DF_mant_len-32+1)) // rounding overflow?
                { manthi = manthi>>1; exp = exp+1; }
        }   }
      union { dfloat eksplicit; double machine_double; } u;
      if ((sintL)exp > (sintL)(DF_exp_high-DF_exp_mid))
        { u.eksplicit.semhi =
              ((sint32)sign & bit(31))
            | ((uint32)(bit(DF_exp_len)-1) << (DF_mant_len-32)); // Infinity
          u.eksplicit.mlo = 0;
        }
      else
        { u.eksplicit.semhi =
              ((sint32)sign & bit(31))                              /* Vorzeichen */
            | ((uint32)((sintL)exp+DF_exp_mid) << (DF_mant_len-32)) /* Exponent   */
            | ((uint32)manthi & (bit(DF_mant_len-32)-1));           /* Mantisse   */
          u.eksplicit.mlo = mantlo;
        }
      return u.machine_double;
      #endif
}

}  // namespace cln
