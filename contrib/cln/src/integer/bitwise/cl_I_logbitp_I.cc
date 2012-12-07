// logbitp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "cln/io.h"
#include "cln/integer_io.h"
#include "cln/exception.h"
#include <sstream>

namespace cln {

bool logbitp (const cl_I& x, const cl_I& y)
{
    // Methode:
    // Falls x<0, Error.
    // Falls x>=0: Falls x>=intDsize*Länge(y), teste Vorzeichen von y.
    //             Sonst x=intDsize*k+i, Teste Bit i vom Worte Nr. k+1 (von oben herab).
      if (!minusp(x)) // x>=0 ?
        { if (fixnump(x))
            { var uintV x_ = FN_to_V(x);
              var uintC ylen;
              var const uintD* yLSDptr;
              I_to_NDS_nocopy(y, ,ylen=,yLSDptr=,true, { return false; } ); // DS zu y
              if (x_ < intDsize*ylen)
                // x ist ein Fixnum >=0, < intDsize*ylen
                { if (lspref(yLSDptr,floor(x_,intDsize)) & bit(x_%intDsize))
                    return true;
                    else
                    return false;
            }   }
          // Vorzeichen von y testen
          if (minusp(y))
            return true;
            else
            return false;
        }
        else
        // x<0
        { std::ostringstream buf;
          fprint(buf, "logbitp: Index is negative: ");
          fprint(buf, x);
          throw runtime_exception(buf.str());
        }
}

}  // namespace cln
