// UDS_to_digits().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "base/digitseq/cl_DS.h"
#include "integer/conv/cl_I_cached_power.h"

namespace cln {

// Timing für Dezimal-Umwandlung einer Zahl mit N Digits = (N*32) Bits,
// auf einem i486 33 MHz unter Linux:
//    N      standard  dnq(div)  dnq(mul)  combined
//     10     0.00031   0.00043   0.00059   0.00031
//     25     0.00103   0.00125   0.00178   0.00103
//     50     0.0030    0.0034    0.0051    0.0030
//    100     0.0100    0.0108    0.0155    0.0100
//    250     0.054     0.055     0.064     0.054
//    500     0.207     0.209     0.229     0.207
//    750     0.47      0.48      0.47      0.47
//   1000     0.81      0.81      0.86      0.81
//   1250     1.25      1.12      1.20      1.12
//   1500     1.81      1.60      1.64      1.61
//   1750     2.45      2.24      2.15      2.25
//   1940     3.01      3.03      3.12      2.80
//   2000     3.20      3.11      3.30      2.89
//   2500     5.00      4.11      4.38      3.91
//   3000     7.3       5.8       5.7       5.5
//   4000    13.0      12.4      12.9       9.7
//   5000    20.3      15.3      15.1      12.4
//  10000    81.4      57.8      56.4      32.5
//  25000                                 112
//  50000                                 265
// dnq(div) means divide-and-conquer using division by B at the topmost call,
//                threshold = 1015.
// dnq(mul) means divide-and-conquer using multiplication by 1/B at the topmost
//                call, threshold = 2050.
// combined means divide-and-conquer as long as length >= threshold.
  const unsigned int cl_digits_div_threshold = 1015;
  const int cl_digits_algo = 1;

// like I_to_digits, except that the result has exactly erg_len characters.
static inline void I_to_digits_noshrink (const cl_I& X, uintD base, uintC erg_len, cl_digits* erg)
{
  I_to_digits(X,base,erg);
  if (erg->len > erg_len) throw runtime_exception();
  var uintC count = erg_len - erg->len;
  if (count > 0)
    { var uintB* ptr = erg->MSBptr;
      do { *--ptr = '0'; } while (--count > 0);
      erg->MSBptr = ptr; erg->len = erg_len;
    }
}

void I_to_digits (const cl_I& X, uintD base, cl_digits* erg)
{
// Methode:
// Umwandlung ins Stellensystem der Basis b geht durch Umwandlung ins Stellen-
// system der Basis b^k (k>=1, b^k<2^intDsize, k maximal) vor sich.
// Aufsuchen von k und b^k aus einer Tabelle.
// Reduktion der UDS zu einer NUDS X.
// Falls X=0: die eine Ziffer 0.
// Falls X>0:
//   Dividiere X durch das Wort b^k,
//   (Single-Precision-Division, vgl. UDS_DIVIDE mit n=1:
//     r:=0, j:=m=Länge(X),
//     while j>0 do
//       j:=j-1, r:=r*beta+X[j], X[j]:=floor(r/b^k), r:=r-b^k*q[j].
//     r=Rest.)
//   zerlege den Rest (mit k-1 Divisionen durch b) in k Ziffern, wandle diese
//   Ziffern einzeln in Ascii um und lege sie an die DIGITS an.
//   Teste auf Speicherüberlauf.
//   X := Quotient.
//   Mache aus X wieder eine NUDS (maximal 1 Nulldigit streichen).
//   Dies solange bis X=0.
//   Streiche die führenden Nullen.
      // Aufsuchen von k-1 und b^k aus der Tabelle:
      var const power_table_entry* tableptr = &power_table[base-2];
      var uintC k = tableptr->k;
      var uintD b_hoch_k = tableptr->b_to_the_k; // b^k
      var uintB* erg_ptr = erg->LSBptr;
      #define next_digit(d)  { *--erg_ptr = (d<10 ? '0'+d : 'A'-10+d); }
      // Spezialfälle:
      if (zerop(X))
        { next_digit(0); goto fertig; } // 0 -> eine Ziffer '0'
      else if ((base & (base-1)) == 0)
        { // Schneller Algorithmus für Zweierpotenzen
          var const uintD* MSDptr;
          var uintC len;
          var const uintD* LSDptr;
          I_to_NDS_nocopy(X, MSDptr=,len=,LSDptr=,false,);
          var int b = (base==2 ? 1 : base==4 ? 2 : base==8 ? 3 : base==16 ? 4 : /*base==32*/ 5);
          var uintD carry = 0;
          var int carrybits = 0;
          loop
            { if (fixnump(X) && erg->LSBptr-erg_ptr>=cl_value_len)
                break;
              if (carrybits >= b)
                { var uintD d = carry & (base-1);
                  next_digit(d);
                  carry = carry >> b; carrybits -= b;
                }
                else
                { var uintD d = carry;
                  if (LSDptr != MSDptr)
                    { carry = lsprefnext(LSDptr);
                      d |= (carry << carrybits) & (base-1);
                      next_digit(d);
                      carry = carry >> (b-carrybits); carrybits = intDsize - (b-carrybits);
                    }
                    else
                    { next_digit(d); break; }
                }
            }
        }
      else if (fixnump(X) || TheBignum(X)->length < cl_digits_div_threshold
               || !cl_digits_algo)
        { // Standard-Algorithmus
          CL_ALLOCA_STACK;
          var uintD* MSDptr;
          var uintC len;
          var uintD* LSDptr;
          I_to_NDS(X, MSDptr=,len=,LSDptr=);
          // normalisiere zu einer NUDS:
          if (mspref(MSDptr,0)==0) { msshrink(MSDptr); len--; }
          loop
            { // Noch die NUDS MSDptr/len/.. mit len>0 abzuarbeiten.
              // Single-Precision-Division durch b^k:
              var uintD rest = divu_loop_msp(b_hoch_k,MSDptr,len);
              // Zerlegen des Restes in seine k Ziffern:
              var uintC count = k-1;
	      if (fixnump(X) && count>cl_value_len-1)
		  count = cl_value_len-1;
              if ((intDsize>=11) || (count>0))
                // (Bei intDsize>=11 ist wegen b<=36 zwangsläufig
                // k = ceiling(intDsize*log(2)/log(b))-1 >= 2, also count = k-1 > 0.)
                do { var uintD d;
                     #if HAVE_DD
                       divuD((uintDD)rest,base,rest=,d=);
                     #else
                       divuD(0,rest,base,rest=,d=);
                     #endif
                     next_digit(d);
                } until (--count == 0);
              next_digit(rest); // letzte der k Ziffern ablegen
              // Quotienten normalisieren (max. 1 Digit streichen):
              if (mspref(MSDptr,0)==0) { msshrink(MSDptr); len--; if (len==0) break; }
        }   }
      else
        { // Divide-and-conquer:
          // Find largest i such that B = base^(k*2^i) satisfies B <= X.
          // Divide by B: X = X1*B + X0. Convert X0 to string, fill up
          // for k*2^i characters, convert X1 to string. (Have to convert
          // X0 first because the conversion may temporarily prepend some
          // zero characters.)
          var uintC ilen_X = integer_length(X);
          var const cached_power_table_entry * p;
          var uintC ilen_B;
          var uintL i;
          for (i = 0; ; i++)
            { p = cached_power(base,i);
              ilen_B = integer_length(p->base_pow);
              if (2*ilen_B >= ilen_X) break;
              // 2*ilen_B < ilen_X, so certainly B^2 < X, let's continue with i+1.
            }
          // 2*ilen_B >= ilen_X, implies X < 2*B^2.
          // Of course also X >= B, implies ilen_X >= ilen_B.
          #ifdef MUL_REPLACES_DIV
          // Divide by B by computing
          //   q := floor((X * floor(2^ilen_X/B)) / 2^ilen_X).
          // We have q <= floor(X/B) <= q+1, so we may have to increment q.
          // Note also that
          // floor(2^ilen_X/B) = floor(floor(2^(2*ilen_B)/B)/2^(2*ilen_B-ilen_X))
          var cl_I q = (X * (p->inv_base_pow >> (2*ilen_B-ilen_X))) >> ilen_X;
          var cl_I r = X - q * p->base_pow;
          if (r < 0) throw runtime_exception();
          if (r >= p->base_pow)
            { q = q+1; r = r - p->base_pow;
              if (r >= p->base_pow) throw runtime_exception();
            }
          #else
          var cl_I_div_t q_r = floor2(X,p->base_pow);
          var const cl_I& q = q_r.quotient;
          var const cl_I& r = q_r.remainder;
          #endif
          var const cl_I& X1 = q;
          var const cl_I& X0 = r;
          var uintC B_baselen = (uintC)(k)<<i;
          I_to_digits_noshrink(X0,base,B_baselen,erg);
          erg->LSBptr -= B_baselen;
          I_to_digits(X1,base,erg);
          erg->LSBptr += B_baselen;
          erg_ptr = erg->MSBptr;
        }
      #undef next_digit
      // Streiche führende Nullen:
      while (*erg_ptr == '0') { erg_ptr++; }
      fertig:
      erg->MSBptr = erg_ptr;
      erg->len = erg->LSBptr - erg_ptr;
}
// Bit complexity (N := length(X)): O(log(N)*M(N)).

}  // namespace cln
