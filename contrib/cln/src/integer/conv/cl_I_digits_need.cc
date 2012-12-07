// cl_digits_need().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

namespace cln {

uintC cl_digits_need (const cl_I& x, uintL base)
{
  if (fixnump(x))
    { return cl_value_len; } // x < 2^cl_value_len, base >= 2, also reicht das
  else
    { var uintC len = TheBignum(x)->length;
      // 1+ceiling(len * intDsize*log(2)/log(base)) Bytes oder etwas mehr
      var uintC need = 1+floor(len,1024/intDsize); // > ceiling(len*intDsize/1024) >= 0
      switch (base) // need mit ceiling(1024*log(2)/log(base)) multiplizieren:
        { case 2: need = 1024*need; break;
          case 3: need = 647*need; break;
          case 4: need = 512*need; break;
          case 5: need = 442*need; break;
          case 6: need = 397*need; break;
          case 7: need = 365*need; break;
          case 8: need = 342*need; break;
          case 9: need = 324*need; break;
          case 10: need = 309*need; break;
          case 11: need = 297*need; break;
          case 12: need = 286*need; break;
          case 13: need = 277*need; break;
          case 14: need = 269*need; break;
          case 15: need = 263*need; break;
          case 16: need = 256*need; break;
          case 17: need = 251*need; break;
          case 18: need = 246*need; break;
          case 19: need = 242*need; break;
          case 20: need = 237*need; break;
          case 21: need = 234*need; break;
          case 22: need = 230*need; break;
          case 23: need = 227*need; break;
          case 24: need = 224*need; break;
          case 25: need = 221*need; break;
          case 26: need = 218*need; break;
          case 27: need = 216*need; break;
          case 28: need = 214*need; break;
          case 29: need = 211*need; break;
          case 30: need = 209*need; break;
          case 31: need = 207*need; break;
          case 32: need = 205*need; break;
          case 33: need = 203*need; break;
          case 34: need = 202*need; break;
          case 35: need = 200*need; break;
          case 36: need = 199*need; break;
          default: NOTREACHED
        }
      // Nun gilt need >= len*intDsize*log(2)/log(base).
      need += 1; // Platzbedarf in Bytes
      return need;
    }
}

}  // namespace cln
