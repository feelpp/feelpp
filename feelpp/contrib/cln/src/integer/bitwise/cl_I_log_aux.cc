// I_to_DS_n_aux().
// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/bitwise/cl_I_log.h"


// Implementation.

#include "cln/number.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

uintD* I_to_DS_n_aux (const cl_I& obj, uintC n, uintD* destptr)
    { // Nun sind unterhalb von destptr n Digits Platz.
      // oberen Teil der DS aus obj füllen, dabei destptr erniedrigen:
      if (fixnump(obj))
        // Fixnum:
        {
          #if (intDsize==64) // && (FN_maxlength==1)
           lsprefnext(destptr) = FN_to_Q(obj);
          #else // (intDsize<=32)
           var uintV wert = FN_to_V(obj);
           #define FN_maxlength_a  (intVsize/intDsize)
           #define FN_maxlength_b  (FN_maxlength<=FN_maxlength_a ? FN_maxlength : FN_maxlength_a)
           // FN_maxlength Digits ablegen. Davon kann man FN_maxlength_b Digits aus wert nehmen.
           #if (FN_maxlength_b > 1)
           doconsttimes(FN_maxlength_b-1,
             lsprefnext(destptr) = (uintD)wert; wert = wert >> intDsize;
             );
           #endif
           lsprefnext(destptr) = (uintD)wert;
           #if (FN_maxlength > FN_maxlength_b)
           // Es ist cl_value_len-1 = intVsize, brauche
           // noch FN_maxlength-FN_maxlength_b = 1 Digit.
           lsprefnext(destptr) = (sintD)sign_of(FN_to_V(obj));
           #endif
          #endif
          n -= FN_maxlength;
        }
        else
        // Bignum:
        { var uintC len = TheBignum(obj)->length;
          n -= len;
          destptr = copy_loop_lsp(BN_LSDptr(obj),destptr,len); // DS kopieren
        }
      // unteren Teil mit Fülldigits, gebildet aus dem Vorzeichen, füllen:
      if (!(n==0))
        { destptr = fill_loop_lsp(destptr,n,sign_of_sintD(mspref(destptr,0))); }
      // destptr zeigt nun aufs untere Ende der DS.
      return destptr;
    }

}  // namespace cln

