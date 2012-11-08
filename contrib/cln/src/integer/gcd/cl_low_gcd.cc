// gcd().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

namespace cln {

// Liefert den ggT zweier Integers.
// gcd(a,b)
// > a,b: zwei Integers
// < ergebnis: (gcd a b), ein Integer >=0
  uintV gcd (uintV a, uintV b)
// binäre Methode:
// (gcd a b) :==
//   (prog ((j 0))
//     1 {a,b >0}
//       (when (oddp a) (if (oddp b) (go 2) (go 4)))
//       (when (oddp b) (go 3))
//       (incf j) (setq a (/ a 2)) (setq b (/ b 2))
//       (go 1)
//     2 {a,b >0, beide ungerade}
//       (cond ((> a b) (setq a (- a b)) (go 3))
//             ((= a b) (go 5))
//             ((< a b) (setq b (- b a)) (go 4))
//       )
//     3 {a,b >0, a gerade, b ungerade}
//       (repeat (setq a (/ a 2)) (until (oddp a)))
//       (go 2)
//     4 {a,b >0, a ungerade, b gerade}
//       (repeat (setq b (/ b 2)) (until (oddp b)))
//       (go 2)
//     5 {a=b>0}
//       (return (ash a j))
//   )
// Statt j zu erhöhen und immer Bit 0 von a und b abfragen,
// fragen wir stattdessen immer Bit j von a und b ab; Bits j-1..0 sind =0.
{
      #ifdef DUMMER_GGT // so macht's ein Mathematiker:
      if (a==0) { return b; }
      if (b==0) { return a; }
      var uintV bit_j = bit(0);
      loop
        { // a,b >0
          if (!((a & bit_j) ==0))
            { if (!((b & bit_j) ==0)) goto odd_odd; else goto odd_even; }
          if (!((b & bit_j) ==0)) goto even_odd;
          // a,b >0 gerade
          bit_j = bit_j<<1;
        }
      #else // Trick von B. Degel:
      var uintV bit_j = (a | b); // endet mit einer 1 und j Nullen
      bit_j = bit_j ^ (bit_j - 1); // Maske = bit(j) | bit(j-1) | ... | bit(0)
      if (!((a & bit_j) ==0))
        { if (!((b & bit_j) ==0))
            goto odd_odd;
          else
            if (b==0)
              return a;
            else
              goto odd_even;
        }
      if (!((b & bit_j) ==0))
        { if (a==0)
            return b;
          else
            goto even_odd;
        }
      return 0; // a=b=0 -> Ergebnis 0
      #endif
      loop
        { odd_odd: // a,b >0, beide ungerade
          // Vergleiche a und b:
          if (a == b) break; // a=b>0 -> fertig
          if (a > b) // a>b ?
            { a = a-b;
              even_odd: // a,b >0, a gerade, b ungerade
              do { a = a>>1; } while ((a & bit_j) ==0);
            }
            else // a<b
            { b = b-a;
              odd_even: // a,b >0, a ungerade, b gerade
              do { b = b>>1; } while ((b & bit_j) ==0);
            }
        }
      // a=b>0
      return a;
}

}  // namespace cln
