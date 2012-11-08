// isqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/cl_low.h"


// Implementation.

namespace cln {

// Zieht die Ganzzahl-Wurzel aus einer 32-Bit-Zahl und
// liefert eine 16-Bit-Wurzel.
// isqrt(x)
// > uintL x : Radikand, >=0, <2^32
// < uintL ergebnis : Wurzel, >=0, <2^16
uintL isqrt (uintL x)
{
  // Methode:
  // x=0 -> y=0, fertig.
  // y := 2^k als Anfangswert, wobei k>0, k<=16 mit 2^(2k-2) <= x < 2^(2k) sei.
  // y := floor((y + floor(x/y))/2) als nächster Wert,
  // solange z := floor(x/y) < y, setze y := floor((y+z)/2).
  // y ist fertig.
  // (Beweis:
  //  1. Die Folge der y ist streng monoton fallend.
  //  2. Stets gilt y >= floor(sqrt(x)) (denn für alle y>0 ist
  //     y + x/y >= 2*sqrt(x) und daher  floor((y + floor(x/y))/2) =
  //     floor(y/2 + x/(2*y)) >= floor(sqrt(x)) ).
  //  3. Am Schluß gilt x >= y^2.
  // )
     if (x==0) return 0; // x=0 -> y=0
     { var uintC k2; integerlength32(x,k2=); // 2^(k2-1) <= x < 2^k2
      {var uintC k1 = floor(k2-1,2); // k1 = k-1, k wie oben
       if (k1 < 16-1)
         // k < 16
         { var uintL y = (x >> (k1+2)) | bit(k1); // stets 2^(k-1) <= y < 2^k
           loop
             { var uintL z;
               divu_3216_1616(x,y, z=,); // Dividiere x/y (geht, da x/y < 2^(2k)/2^(k-1) = 2^(k+1) <= 2^16)
               if (z >= y) break;
               y = floor(z+y,2); // geht, da z+y < 2*y < 2^(k+1) <= 2^16
             }
           return y;
         }
         else
         // k = 16, Vorsicht!
         { var uintL x1 = high16(x);
           var uintL y = (x >> (16+1)) | bit(16-1); // stets 2^(k-1) <= y < 2^k
           loop
             { var uintL z;
               if (x1 >= y) break; // Division x/y ergäbe Überlauf -> z > y
               divu_3216_1616(x,y, z=,); // Dividiere x/y
               if (z >= y) break;
               y = floor(z+y,2);
             }
           return y;
         }
     }}
}

#ifdef HAVE_LONGLONG
uintL isqrt (uintQ x)
{
  // As isqrt (uintL) above, but with 64-bit numbers.
     if (x==0) return 0; // x=0 -> y=0
     { var uintC k2; integerlength64(x,k2=); // 2^(k2-1) <= x < 2^k2
      {var uintC k1 = floor(k2-1,2); // k1 = k-1, k as above
       if (k1 < 32-1)
         // k < 32
         { var uintL y = (x >> (k1+2)) | bit(k1); // always 2^(k-1) <= y < 2^k
           loop
             { var uintL z;
               divu_6432_3232(high32(x),low32(x),y, z=,); // z := x/y (works, since x/y < 2^(2k)/2^(k-1) = 2^(k+1) <= 2^32)
               if (z >= y) break;
               y = floor(z+y,2); // geht, da z+y < 2*y < 2^(k+1) <= 2^32
             }
           return y;
         }
         else
         // k = 32, careful!
         { var uintL x1 = high32(x);
           var uintL y = (x >> (32+1)) | bit(32-1); // stets 2^(k-1) <= y < 2^k
           loop
             { var uintL z;
               if (x1 >= y) break; // division x/y would overflow -> z > y
               divu_6432_3232(high32(x),low32(x),y, z=,); // divide x/y
               if (z >= y) break;
               y = floor(z+y,2);
             }
           return y;
         }
     }}
}
#endif

}  // namespace cln
