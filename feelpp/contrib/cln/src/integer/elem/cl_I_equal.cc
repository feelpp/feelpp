// equal().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

bool equal (const cl_I& x, const cl_I& y)
{
// Methode:
// x und y haben gleiches Vorzeichen ->
// x Fixnum ->
//    y Fixnum -> direkt vergleichen.
//    y Bignum -> verschieden.
// x Bignum ->
//    y Fixnum -> verschieden.
//    y Bignum ->
//       falls beide gleich lang, wortweise vergleichen, sonst verschieden.
      if (fixnump(x))
        // x Fixnum
        if (fixnump(y))
          // x Fixnum, y Fixnum
          { // This assumes cl_value_shift + cl_value_len == cl_pointer_size.
            return (cl_sint)x.word == (cl_sint)y.word;
          }
          else
          // x Fixnum, y Bignum
          return false;
        else
        // x Bignum
        if (fixnump(y))
          // x Bignum, y Fixnum
          return false;
          else
          // x Bignum, y Bignum
          if (x.pointer == y.pointer)
            return true; // gleiche Pointer -> selbe Zahl
            else
            { var uintC xlen = TheBignum(x)->length;
              var uintC ylen = TheBignum(y)->length;
              if (xlen==ylen)
                // gleiche Länge -> digitweise vergleichen
                return compare_loop_msp(BN_MSDptr(x),BN_MSDptr(y),xlen) == 0;
                else
                return false;
            }
}

}  // namespace cln
