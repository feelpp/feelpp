// isqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/cl_low.h"


// Implementation.

namespace cln {

// Zieht die Ganzzahl-Wurzel aus einer 64-Bit-Zahl und
// liefert eine 32-Bit-Wurzel.
// isqrt(x1,x0)
// > uintL2 x = x1*2^32+x0 : Radikand, >=0, <2^64
// < uintL ergebnis : Wurzel, >=0, <2^32
uintL isqrt (uintL x1, uintL x0)
{
  // Methode:
  // x=0 -> y=0, fertig.
  // y := 2^k als Anfangswert, wobei k>0, k<=32 mit 2^(2k-2) <= x < 2^(2k) sei.
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
     if (x1==0) { return isqrt(x0); } // x klein?
     { var uintC k2; integerlength32(x1,k2=); // 2^(k2+32-1) <= x < 2^(k2+32)
      {var uintC k = ceiling(k2+32,2); // k wie oben
       if (k < 32)
         // k < 32
         { var uintL y = ((x1 << (32-k)) | (x0 >> k) | bit(k)) >> 1; // stets 2^(k-1) <= y < 2^k
           loop
             { var uintL z;
               divu_6432_3232(x1,x0,y, z=,); // Dividiere x/y (geht, da x/y < 2^(2k)/2^(k-1) = 2^(k+1) <= 2^32)
               if (z >= y) break;
               y = floor(z+y,2); // geht, da z+y < 2*y < 2^(k+1) <= 2^32
             }
           return y;
         }
         else
         // k = 32, Vorsicht!
         { var uintL y = (x1 >> 1) | bit(32-1); // stets 2^(k-1) <= y < 2^k
           loop
             { var uintL z;
               if (x1 >= y) break; // Division x/y ergäbe Überlauf -> z > y
               divu_6432_3232(x1,x0,y, z=,); // Dividiere x/y
               if (z >= y) break;
               y = floor(z+y,2) | bit(32-1); // y muß >= 2^(k-1) bleiben
             }
           return y;
         }
     }}
}

}  // namespace cln
