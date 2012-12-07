// NDS_to_I().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

CL_INLINE const cl_I CL_INLINE_DECL(NDS_to_I) (const uintD* MSDptr, uintC len)
{
      // Mehr als bn_minlength Digits -> Bignum.
      // Weniger als bn_minlength Digits -> Fixnum.
      // Genau bn_minlength Digits -> Bignum oder Fixnum.
      if (len < bn_minlength)
        { // 0..bn_minlength-1 Digits, paßt in ein Fixnum:
          if (bn_minlength>1 ? (len==0) : TRUE)
            // 0 Digits
            { return 0; }
          #if (cl_word_size <= 32)
          var sint32 wert;
          if (bn_minlength>2 ? (len==1) : TRUE)
            // 1 Digit
            len_1:
            { wert = get_sint1D_Dptr(MSDptr mspop 1); }
          elif (bn_minlength>3 ? (len==2) : TRUE)
            // 2 Digits
            len_2:
            { wert = get_sint2D_Dptr(MSDptr mspop 2); }
          elif (bn_minlength>4 ? (len==3) : TRUE)
            // 3 Digits
            len_3:
            { wert = get_sint3D_Dptr(MSDptr mspop 3); }
          elif (TRUE)
            // 4 Digits
            len_4:
            { wert = get_sint4D_Dptr(MSDptr mspop 4); }
          elif (FALSE)
            // 5 Digits
            len_5:
            { wert = get_sint4D_Dptr(MSDptr mspop 5); }
          return L_to_FN(wert);
          #else // (cl_word_size==64)
          var sint64 wert;
          #if (intDsize==32)
          if (TRUE)
            // 1 Digit
            len_1:
            { wert = (sint64)(sintD)mspref(MSDptr,0); }
          elif (FALSE)
            // 2 Digits
            len_2:
            { wert = ((sint64)(sintD)mspref(MSDptr,0) << intDsize) | (uint64)(uintD)mspref(MSDptr,1); }
          #endif
          #if (intDsize==64)
          if (FALSE)
            // 1 Digit
            len_1:
            { wert = (sintD)mspref(MSDptr,0); }
          #endif
          return cl_I_from_word(cl_combine(cl_FN_tag,wert));
          #endif
        }
      #if (cl_value_len > (bn_minlength-1)*intDsize)
      if (len == bn_minlength)
        // bn_minlength Digits, also (incl. Vorzeichen) zwischen
        // (bn_minlength-1)*intDsize+1 und bn_minlength*intDsize Bits.
        // Höchstens cl_value_len Bits -> paßt in ein Fixnum:
        { if (  (mspref(MSDptr,0) <= (uintD)(bit(cl_value_len-1-(bn_minlength-1)*intDsize)-1)) // Fixnum >=0 ?
              ||(mspref(MSDptr,0) >= (uintD)(-bit(cl_value_len-1-(bn_minlength-1)*intDsize))) // Fixnum <0 ?
             )
            #if (bn_minlength==1)
            goto len_1;
            #endif
            #if (bn_minlength==2)
            goto len_2;
            #endif
            #if (bn_minlength==3)
            goto len_3;
            #endif
            #if (bn_minlength==4)
            goto len_4;
            #endif
            #if (bn_minlength==5)
            goto len_5;
            #endif
        }
      #endif
      // mindestens bn_minlength Digits, mache ein Bignum
      { var Bignum result = allocate_bignum(len);
        // neues Bignum mit dem Inhalt der NDS füllen:
        copy_loop_msp(MSDptr,arrayMSDptr(result->data,len),len);
        return result;
      }
}

}  // namespace cln
