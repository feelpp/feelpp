// partial_gcd().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/integer.h"
#include "base/digit/cl_D.h"

namespace cln {

void partial_gcd (uintD z1, uintD z2, partial_gcd_result* erg)
  { var uintD x1 = 1;
    var uintD y1 = 0;
    var uintD x2 = 0;
    var uintD y2 = 1;
    for (;;)
      {
        // Hier ist z1-y1>=z2+y2.
        // Bestimme q := floor((z1-y1)/(z2+y2)) >= 1 :
        { var uintD zaehler = z1-y1;
          var uintD nenner = z2+y2; // z2+y2 <= z1-y1 < beta !
          if (floor(zaehler,8) >= nenner) // zaehler >= 8*nenner ?
            // ja -> Dividieren lohnt sich wohl
            { var uintD q = floorD(zaehler,nenner);
              x1 += muluD_unchecked(q,x2); // x1 := x1+q*x2
              y1 += muluD_unchecked(q,y2); // y1 := y1+q*y2
              z1 -= muluD_unchecked(q,z2); // z1 := z1-q*z2
            }
            else
            // nein -> ein paarmal subtrahieren ist wohl schneller
            do { x1 += x2; y1 += y2; z1 -= z2; } // (x1,y1,z1) := (x1+x2,y1+y2,z1-z2)
               while (z1-y1 >= nenner);
        }
        if (z2-x2 <= z1+x1-1) break;
        // Hier ist z2-x2>=z1+x1.
        // Bestimme q := floor((z2-x2)/(z1+x1)) >= 1 :
        { var uintD zaehler = z2-x2;
          var uintD nenner = z1+x1; // z1+x1 <= z2-x2 < beta !
          if (floor(zaehler,8) >= nenner) // zaehler >= 8*nenner ?
            // ja -> Dividieren lohnt sich wohl
            { var uintD q = floorD(zaehler,nenner);
              x2 += muluD_unchecked(q,x1); // x2 := x2+q*x1
              y2 += muluD_unchecked(q,y1); // y2 := y2+q*y1
              z2 -= muluD_unchecked(q,z1); // z2 := z2-q*z1
            }
            else
            // nein -> ein paarmal subtrahieren ist wohl schneller
            do { x2 += x1; y2 += y1; z2 -= z1; } // (x2,y2,z2) := (x2+x1,y2+y1,z2-z1)
               while (z2-x2 >= nenner);
        }
        if (z1-y1 <= z2+y2-1) break;
      }
    // Keine Subtraktion mehr möglich.
    erg->x1 = x1; erg->y1 = y1; erg->x2 = x2; erg->y2 = y2; // Ergebnis
  }

}  // namespace cln
