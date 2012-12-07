// logtest().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

bool logtest (const cl_I& x, const cl_I& y)
{
    // Methode:
    //  Fixnums separat behandeln.
    //  Sei oBdA x die kürzere der beiden Zahlen (in Digits).
    //  x echt kürzer und x<0 -> [eines der most signif. intDsize+1 Bits von y ist 1] Ja.
    //  Beide gleich lang oder x>=0 ->
    //   Kann mich auf die untersten length(x) Digits beschraenken.
    //   Mit AND durchlaufen, abbrechen (mit "Ja") falls /=0. Am Ende: Nein.
      if (fixnump(x))
        if (fixnump(y))
          // beides Fixnums
          { if ((x.word & y.word & cl_combine(0,~(cl_uint)0))==0)
              return false;
              else
              return true;
          }
          else
          // x Fixnum, y Bignum, also ist x echt kürzer
          { if (FN_V_minusp(x,FN_to_V(x))) return true; // x<0 -> ja.
            // x>=0. Kombiniere x mit den pFN_maxlength letzten Digits von y.
           {var const uintD* yLSDptr;
            var uintV x_ = FN_to_V(x);
            BN_to_NDS_nocopy(y, ,,yLSDptr=);
            #if (pFN_maxlength > 1)
            doconsttimes(pFN_maxlength-1,
              if (lsprefnext(yLSDptr) & (uintD)x_) return true;
              x_ = x_ >> intDsize;
              );
            #endif
            if (lsprefnext(yLSDptr) & (uintD)x_) return true;
            return false;
          }}
        else
        if (fixnump(y))
          // x Bignum, y Fixnum, analog wie oben, nur x und y vertauscht
          { if (FN_V_minusp(y,FN_to_V(y))) return true; // y<0 -> ja.
            // y>=0. Kombiniere y mit den pFN_maxlength letzten Digits von x.
           {var const uintD* xLSDptr;
            var uintV y_ = FN_to_V(y);
            BN_to_NDS_nocopy(x, ,,xLSDptr=);
            #if (pFN_maxlength > 1)
            doconsttimes(pFN_maxlength-1,
              if (lsprefnext(xLSDptr) & (uintD)y_) return true;
              y_ = y_ >> intDsize;
              );
            #endif
            if (lsprefnext(xLSDptr) & (uintD)y_) return true;
            return false;
          }}
          else
          // x,y Bignums
          { var const uintD* xMSDptr;
            var uintC xlen;
            var const uintD* yMSDptr;
            var uintC ylen;
            BN_to_NDS_nocopy(x, xMSDptr=,xlen=,);
            BN_to_NDS_nocopy(y, yMSDptr=,ylen=,);
            // Beachte: xlen>0, ylen>0.
            if (!(xlen==ylen))
              // beide verschieden lang
              { if (xlen<ylen)
                  { // x ist die echt kürzere DS.
                    if ((sintD)mspref(xMSDptr,0)<0) // der echt kürzere ist negativ?
                      return true;
                    // Der echt kürzere ist positiv.
                    yMSDptr = yMSDptr mspop (ylen-xlen);
                  }
                  else
                  { // y ist die echt kürzere DS.
                    if ((sintD)mspref(yMSDptr,0)<0) // der echt kürzere ist negativ?
                      return true;
                    // Der echt kürzere ist positiv.
                    xMSDptr = xMSDptr mspop (xlen-ylen);
                    xlen = ylen;
              }   }
            // Nach gemeinsamen Bits in xMSDptr/xlen/.. und yMSDptr/xlen/..
            // suchen:
            return and_test_loop_msp(xMSDptr,yMSDptr,xlen);
          }
}

}  // namespace cln
