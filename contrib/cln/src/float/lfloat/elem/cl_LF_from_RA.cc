// cl_RA_to_LF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "float/lfloat/cl_LF_impl.h"
#include "rational/cl_RA.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "float/cl_F.h"

namespace cln {

const cl_LF cl_RA_to_LF (const cl_RA& x, uintC len)
{
// Methode:
// x ganz -> klar.
// x = +/- a/b mit Integers a,b>0:
//   Sei k,m so gewählt, daß
//     2^(k-1) <= a < 2^k, 2^(m-1) <= b < 2^m.
//   Dann ist 2^(k-m-1) < a/b < 2^(k-m+1).
//   Ergebnis-Vorzeichen := Vorzeichen von x.
//   Berechne k=(integer-length a) und m=(integer-length b).
//   Ergebnis-Exponent := k-m.
//   Ergebnis-Mantisse:
//     Berechne floor(2^(-k+m+16n+1)*a/b) :
//       Bei k-m>=16n+1 dividiere a durch (ash b (k-m-16n-1)),
//       bei k-m<16n+1 dividiere (ash a (-k+m+16n+1)) durch b.
//     Der erste Wert ist >=2^16n, <2^(16n+2).
//     Falls er >=2^(16n+1) ist, erhöhe Exponent um 1,
//       runde 2 Bits weg und schiebe dabei um 2 Bits nach rechts;
//     falls er <2^(16n+1) ist,
//       runde 1 Bit weg und schiebe dabei um 1 Bit nach rechts.
// NB: Wenn a und b länger sind als len, ist dieser Algorithmus weniger
//     effizient, als cl_float(a,len)/cl_float(b,len) zu berechnen. Aber
//     es ist wichtig, dass cl_RA_to_LF nicht mehr als 0.5 ulp Fehler hat,
//     deswegen belassen wir es beim ineffizienten aber exakten Algorithmus.
//     Wenn es auf Rundungsfehler nicht ankommt, muss der Aufrufer im Fall
//           ceiling(integer_length(a),intDsize) >= len
//        && ceiling(integer_length(b),intDsize) >= len
//     einen anderen Algorithmus wählen.
      if (integerp(x)) {
        DeclareType(cl_I,x);
        return cl_I_to_LF(x,len);
      }
 {    // x Ratio
      DeclareType(cl_RT,x);
      var cl_I a = numerator(x); // +/- a
      var const cl_I& b = denominator(x); // b
      var cl_signean sign = -(cl_signean)minusp(a); // Vorzeichen
      if (!(sign==0)) { a = -a; } // Betrag nehmen, liefert a
      var sintC lendiff = (sintC)integer_length(a) // (integer-length a)
                          - (sintC)integer_length(b); // (integer-length b)
      // |lendiff| < intDsize*2^intCsize. Da für LF-Exponenten ein sintL zur
      // Verfügung steht, braucht man keinen Test auf Overflow oder Underflow.
      var uintC difflimit = intDsize*len + 1; // 16n+1
      var cl_I zaehler;
      var cl_I nenner;
      if (lendiff > (sintC)difflimit)
        // 0 <= k-m-16n-1 < k < intDsize*2^intCsize
        { nenner = ash(b,(uintC)(lendiff - difflimit));
          zaehler = a;
        }
        else
        // 0 < -k+m+16n+1 <= m+1 + 16n < intDsize*2^intCsize + intDsize*2^intCsize
        { zaehler = ash(a,(uintC)(difflimit - lendiff)); // (ash a -k+m+16n+1)
          nenner = b; // b
        }
      // Division zaehler/nenner durchführen:
      var cl_I_div_t q_r = cl_divide(zaehler,nenner);
      var cl_I& q = q_r.quotient;
      var cl_I& r = q_r.remainder;
      // 2^16n <= q < 2^(16n+2), also ist q Bignum mit n+1 Digits.
      var Lfloat y = allocate_lfloat(len,lendiff+LF_exp_mid,sign); // neues Long-Float
      var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
      {var uintD* q_MSDptr = arrayMSDptr(TheBignum(q)->data,len+1);
       if (mspref(q_MSDptr,0) == 1) // erstes Digit =1 oder =2,3 ?
         // 2^16n <= q < 2^(16n+1), also 2^(k-m-1) < a/b < 2^(k-m).
         { // Mantisse mit einer Schiebeschleife um 1 Bit nach rechts füllen:
           var uintD rounding_bit =
             shiftrightcopy_loop_msp(q_MSDptr mspop 1,y_mantMSDptr,len,1,1);
           if ( (rounding_bit == 0) // herausgeschobenes Bit =0 -> abrunden
                || ( eq(r,0) // =1 und Rest r > 0 -> aufrunden
                     // round-to-even
                     && ((mspref(y_mantMSDptr,len-1) & bit(0)) ==0)
              )    )
             goto ab; // abrunden
             else
             goto auf; // aufrunden
         }
         else
         // 2^(16n+1) <= q < 2^(16n+2), also 2^(k-m) < a/b < 2^(k-m+1).
         { // Mantisse mit einer Schiebeschleife um 2 Bit nach rechts füllen:
           var uintD rounding_bits =
             shiftrightcopy_loop_msp(q_MSDptr mspop 1,y_mantMSDptr,len,2,mspref(q_MSDptr,0));
           (TheLfloat(y)->expo)++; // Exponenten incrementieren auf k-m+1
           if ( ((sintD)rounding_bits >= 0) // herausgeschobenes Bit =0 -> abrunden
                || ( ((rounding_bits & bit(intDsize-2)) ==0) // =1 und nächstes Bit =1 oder Rest r > 0 -> aufrunden
                     && eq(r,0)
                     // round-to-even
                     && ((mspref(y_mantMSDptr,len-1) & bit(0)) ==0)
              )    )
             goto ab; // abrunden
             else
             goto auf; // aufrunden
         }
      }
      auf: // aufrunden
        { if ( inc_loop_lsp(y_mantMSDptr mspop len,len) )
            // Übertrag durchs Aufrunden
            { mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
              (TheLfloat(y)->expo)++; // Exponenten incrementieren
        }   }
      ab: // abrunden
      return y;
}}

// Timings on an i486 33 MHz, running Linux, in 0.01 sec.
// First timing:  cl_I_to_LF(numerator,len)/cl_I_to_LF(denominator,len)
// Second timing: cl_RA_to_LF(x,len)
// with len = 100.
//     num_length    50          70          100         200         500
// den_length
//
//         50     1.86 0.97   1.84 0.97   1.85 0.96   1.86 1.86   1.85 7.14
//
//         70     1.86 1.33   1.85 1.31   1.85 1.32   1.84 1.84   1.85 7.13
//
//        100     1.85 1.85   1.86 1.85   1.85 1.84   1.84 1.84   1.86 7.13
//
//        200     1.85 3.61   1.84 3.61   1.85 3.59   1.85 3.59   1.87 7.12
//
//        500     1.84 7.44   1.84 7.55   1.85 7.56   1.84 7.66   1.86 7.63
//
// We see that cl_RA_to_LF is faster only if
//            num_length < 2*len && den_length < len
// whereas cl_I_to_LF(numerator,len)/cl_I_to_LF(denominator,len) is faster if
//            num_length > 2*len || den_length > len

}  // namespace cln
