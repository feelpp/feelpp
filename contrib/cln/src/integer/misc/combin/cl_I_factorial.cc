// factorial().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "integer/misc/combin/cl_I_combin.h"

namespace cln {

  // Methode:
  // n <= 10 -> Ergebnis (Fixnum) aus Tabelle
  // Sonst:
  //   Zweierpotenzen extra am Schluß durch einen Shift um
  //   ord2(n!) = sum(k>=1, floor(n/2^k) ) = n - logcount(n)  Bits.
  //   Für k>=1 wird jede ungerade Zahl m im Intervall n/2^k < m <= n/2^(k-1)
  //   genau k mal gebraucht (als ungerader Anteil von m*2^0,...,m*2^(k-1) ).
  //   Zur Bestimmung des Produkts aller ungeraden Zahlen in einem Intervall
  //   a < m <= b verwenden wir eine rekursive Funktion, die nach Divide-and-
  //   Conquer das Produkt über die Intervalle a < m <= c und c < m <= b
  //   (c := floor((a+b)/2)) bestimmt und beide zusammenmultipliziert. Dies
  //   vermeidet, daß oft große Zahlen mit ganz kleinen Zahlen multipliziert
  //   werden.

const cl_I factorial (uintL n) // assume n >= 0 small
{
	static uintV const fakul_table [] = {
        1,
        1UL,
        1UL*2,
        #if (cl_value_len>=4)
        1UL*2*3,
        #if (cl_value_len>=6)
        1UL*2*3*4,
        #if (cl_value_len>=8)
        1UL*2*3*4*5,
        #if (cl_value_len>=11)
        1UL*2*3*4*5*6,
        #if (cl_value_len>=14)
        1UL*2*3*4*5*6*7,
        #if (cl_value_len>=17)
        1UL*2*3*4*5*6*7*8,
        #if (cl_value_len>=20)
        1UL*2*3*4*5*6*7*8*9,
        #if (cl_value_len>=23)
        1UL*2*3*4*5*6*7*8*9*10,
        #if (cl_value_len>=27)
        1UL*2*3*4*5*6*7*8*9*10*11,
        #if (cl_value_len>=30)
        1UL*2*3*4*5*6*7*8*9*10*11*12,
        #if (cl_value_len>=34)
        1UL*2*3*4*5*6*7*8*9*10*11*12*13,
        #if (cl_value_len>=38)
        1UL*2*3*4*5*6*7*8*9*10*11*12*13*14,
        #if (cl_value_len>=42)
        1UL*2*3*4*5*6*7*8*9*10*11*12*13*14*15,
        #if (cl_value_len>=46)
        1UL*2*3*4*5*6*7*8*9*10*11*12*13*14*15*16,
        #if (cl_value_len>=50)
        1UL*2*3*4*5*6*7*8*9*10*11*12*13*14*15*16*17,
        #if (cl_value_len>=54)
        1UL*2*3*4*5*6*7*8*9*10*11*12*13*14*15*16*17*18,
        #if (cl_value_len>=58)
        1UL*2*3*4*5*6*7*8*9*10*11*12*13*14*15*16*17*18*19,
        #if (cl_value_len>=63)
        1UL*2*3*4*5*6*7*8*9*10*11*12*13*14*15*16*17*18*19*20,
        #if (cl_value_len>=67)
        ...
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
        #endif
	};

      if (n < sizeof(fakul_table)/sizeof(cl_I))
        { return UV_to_I(fakul_table[n]); }
        else
        { var cl_I prod = 1; // bisheriges Produkt := 1
          var uintL k = 1;
          var uintL A = n;
          var uintL B = n; // obere Intervallgrenze floor(n/2^(k-1))
          loop
            { // 'A' enthält floor(n/2^(k-1)).
              A = A >> 1; // untere Grenze floor(n/2^k)
              // 'A' enthält floor(n/2^k).
              // Bilde Teilprodukt prod(A < i <= B & oddp(i), i)
              //       = prod(floor((A-1)/2) < i <= floor((B-1)/2), 2*i+1)
              // wobei B = floor(n/2^(k-1)), A = floor(n/2^k) = floor(B/2).
              { var uintL b = floor(B-1,2);
                if (b==0) break; // B=2 oder B=1 -> Produkt fertig
                var uintL a = floor(A-1,2);
                prod = expt_pos(cl_I_prod_ungerade(a,b),k) * prod; // aufmultiplizieren
              }
              k = k+1;
              B = A;
            }
          return prod << (n - logcount(n));
        }
}
// Bit complexity (N := n): O(log(N)^2*M(N)).

}  // namespace cln

