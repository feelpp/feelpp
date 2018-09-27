// DS_to_I().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"
#include "base/digitseq/cl_DS.h"

#include "base/cl_inline.h"
#include "integer/conv/cl_I_from_NDS.cc"

namespace cln {

CL_INLINE2 const cl_I CL_INLINE2_DECL(DS_to_I) (const uintD* MSDptr, uintC len)
{
      // erst normalisieren.
      // Dabei evtl. MSDptr erhöhen und len erniedrigen:
      if (!(len==0)) // leere DS ist normalisiert
        { var uintC count = len-1;
          if ((sintD)mspref(MSDptr,0) >= 0)
            // Zahl >= 0
            { // versuche maximal len-1 führende Nullen-Digits zu streichen:
              while (!(count==0) && (mspref(MSDptr,0)==0) && ((sintD)mspref(MSDptr,1)>=0))
                { msshrink(MSDptr); len--; count--; } // Nulldigit streichen
            }
            else
            // Zahl < 0
            // versuche maximal len-1 führende Einsen-Digits zu streichen:
            { while (!(count==0) && ((sintD)mspref(MSDptr,0)==-1) && ((sintD)mspref(MSDptr,1)<0))
                { msshrink(MSDptr); len--; count--; } // Einsen-digit streichen
        }   }
      // Eventuell ist jetzt noch bei der DS 0 ausnahmsweise len=1,
      // aber NDS_to_I wird auch damit fertig.
      return NDS_to_I_inline(MSDptr,len);
}

}  // namespace cln
