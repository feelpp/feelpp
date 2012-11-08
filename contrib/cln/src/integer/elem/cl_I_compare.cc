// compare().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

cl_signean compare (const cl_I& x, const cl_I& y)
{
      // Methode:
      // x und y haben verschiedenes Vorzeichen ->
      //    x < 0 -> x < y
      //    x >= 0 -> x > y
      // x und y haben gleiches Vorzeichen ->
      // x Fixnum ->
      //    y Fixnum -> direkt vergleichen.
      //    y Bignum ->
      //       y > 0 -> x < y
      //       y < 0 -> x > y
      // x Bignum ->
      //    y Fixnum ->
      //       x < 0 -> x < y
      //       x > 0 -> x > y
      //    y Bignum ->
      //       falls beide gleich lang -> wortweise vergleichen
      //       x kürzer als y -> bei x,y > 0 : x < y, bei x,y < 0 : x > y
      //       y kürzer als x -> bei x,y > 0 : x > y, bei x,y > 0 : x < y
      var uintC xlen;
      var uintC ylen;
      if (fixnump(x))
        // x Fixnum
        if (fixnump(y))
          // x Fixnum, y Fixnum
          { // This assumes cl_value_shift + cl_value_len == cl_pointer_size.
            if ((cl_sint)x.word == (cl_sint)y.word) return signean_null;
            else if ((cl_sint)x.word > (cl_sint)y.word) return signean_plus;
            else return signean_minus;
          }
          else
          // x Fixnum, y Bignum
          if ((sintD)mspref(BN_MSDptr(y),0) >= 0)
            // x Fixnum, y Bignum >0
            return signean_minus; // x<y
            else
            // x Fixnum, y Bignum <0
            return signean_plus; // x>y
        else
        // x Bignum
        if (fixnump(y))
          // x Bignum, y Fixnum
          if ((sintD)mspref(BN_MSDptr(x),0) >= 0)
            // x Bignum >0, y Fixnum
            return signean_plus; // x>y
            else
            // x Bignum <0, y Fixnum
            return signean_minus; // x<y
          else
          // x Bignum, y Bignum
          if ((sintD)mspref(BN_MSDptr(x),0) >= 0)
            // x Bignum >0, y Bignum
            if ((sintD)mspref(BN_MSDptr(y),0) >= 0)
              // x und y Bignums >0
              if (x.pointer == y.pointer)
                return signean_null; // gleiche Pointer -> selbe Zahl
                else
                { xlen = TheBignum(x)->length;
                  ylen = TheBignum(y)->length;
                  if (xlen==ylen)
                    samelength:
                    // gleiche Länge -> digitweise vergleichen
                    return compare_loop_msp(BN_MSDptr(x),BN_MSDptr(y),xlen);
                    else
                    return (xlen > ylen ? signean_plus : signean_minus);
                }
              else
              // x Bignum >0, y Bignum <0
              return signean_plus; // x>y
            else
            // x Bignum <0, y Bignum
            if ((sintD)mspref(BN_MSDptr(y),0) >= 0)
              // x Bignum <0, y Bignum >0
              return signean_minus; // x<y
              else
              // x und y Bignums <0
              if (x.pointer == y.pointer)
                return signean_null; // gleiche Pointer -> selbe Zahl
                else
                { xlen = TheBignum(x)->length;
                  ylen = TheBignum(y)->length;
                  if (xlen==ylen)
                    // gleiche Länge -> wortweise vergleichen
                    goto samelength; // wie oben
                    else
                    return (xlen > ylen ? signean_minus : signean_plus);
                }
}

}  // namespace cln
