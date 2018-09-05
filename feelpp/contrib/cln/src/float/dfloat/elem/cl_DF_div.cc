// binary operator /

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "base/cl_N.h"
#include "float/cl_F.h"
#include "base/cl_low.h"
#include "base/digitseq/cl_DS.h"

#include "base/cl_inline.h"
#include "float/dfloat/elem/cl_DF_zerop.cc"

namespace cln {


const cl_DF operator/ (const cl_DF& x1, const cl_DF& x2)
{
// Methode:
// x2 = 0.0 -> Error
// x1 = 0.0 -> Ergebnis 0.0
// Sonst:
// Ergebnis-Vorzeichen = xor der beiden Vorzeichen von x1 und x2
// Ergebnis-Exponent = Differenz der beiden Exponenten von x1 und x2
// Ergebnis-Mantisse = Mantisse mant1 / Mantisse mant2, gerundet.
//   mant1/mant2 > 1/2, mant1/mant2 < 2;
//   nach Rundung mant1/mant2 >=1/2, <=2*mant1<2.
//   Bei mant1/mant2 >=1 brauche 52 Nachkommabits,
//   bei mant1/mant2 <1 brauche 53 Nachkommabits.
//   Fürs Runden: brauche ein Rundungsbit (Rest gibt an, ob exakt).
//   Brauche daher insgesamt 54 Nachkommabits von mant1/mant2.
//   Dividiere daher (als Unsigned Integers) 2^54*(2^53*mant1) durch (2^53*mant2).
//   Falls der Quotient >=2^54 ist, runde die letzten zwei Bits weg und
//     erhöhe den Exponenten um 1.
//   Falls der Quotient <2^54 ist, runde das letzte Bit weg. Bei rounding
//     overflow schiebe um ein weiteres Bit nach rechts, incr. Exponenten.
#if defined(FAST_DOUBLE) && !defined(__i386__)
      double_to_DF(DF_to_double(x1) / DF_to_double(x2), return ,
                   TRUE, TRUE, // Overflow und subnormale Zahl abfangen
                   !zerop_inline(x1), // ein Ergebnis +/- 0.0
                               // ist genau dann in Wirklichkeit ein Underflow
                   zerop_inline(x2), // Division durch Null abfangen
                   FALSE // kein NaN als Ergebnis möglich
                  );
#else
      // x1,x2 entpacken:
      var cl_signean sign1;
      var sintL exp1;
      #if (intDsize<=32)
      var uintL manthi1;
      var uintL mantlo1;
      #endif
      var cl_signean sign2;
      var sintL exp2;
      #if (intDsize<=32)
      var uintL manthi2;
      var uintL mantlo2;
      #endif
      #if (cl_word_size==64)
      var uint64 mantx1;
      var uint64 mantx2;
      DF_decode(x2, { throw division_by_0_exception(); }, sign2=,exp2=,mantx2=);
      DF_decode(x1, { return x1; }, sign1=,exp1=,mantx1=);
      #else
      DF_decode2(x2, { throw division_by_0_exception(); }, sign2=,exp2=,manthi2=,mantlo2=);
      DF_decode2(x1, { return x1; }, sign1=,exp1=,manthi1=,mantlo1=);
      #endif
      exp1 = exp1 - exp2; // Differenz der Exponenten
      sign1 = sign1 ^ sign2; // Ergebnis-Vorzeichen
      // Dividiere 2^54*mant1 durch mant2 oder (äquivalent)
      // 2^i*2^54*mant1 durch 2^i*mant2 für irgendein i mit 0 <= i <= 64-53 :
      // wähle i = 64-(DF_mant_len+1), also i+(DF_mant_len+2) = 65.
      #if (cl_word_size==64)
      mantx1 = mantx1 << 1;
      mantx2 = mantx2 << (64-(DF_mant_len+1));
      #if (intDsize<=32)
      manthi1 = high32(mantx1); mantlo1 = low32(mantx1);
      manthi2 = high32(mantx2); mantlo2 = low32(mantx2);
      #endif
      #else
      manthi1 = (manthi1 << 1) | (mantlo1 >> 31); mantlo1 = mantlo1 << 1;
      manthi2 = (manthi2 << (64-(DF_mant_len+1))) | (mantlo2 >> ((DF_mant_len+1)-32)); mantlo2 = mantlo2 << (64-(DF_mant_len+1));
      #endif
      var uintD mant1 [128/intDsize];
      var uintD mant2 [64/intDsize];
      #if (intDsize==64)
      arrayLSref(mant1,128/intDsize,1) = mantx1;
      arrayLSref(mant1,128/intDsize,0) = 0;
      arrayLSref(mant2,64/intDsize,0) = mantx2;
      #elif (intDsize==32) || (intDsize==16) || (intDsize==8)
      set_32_Dptr(arrayMSDptr(mant1,128/intDsize),manthi1);
      set_32_Dptr(arrayMSDptr(mant1,128/intDsize) mspop 32/intDsize,mantlo1);
      set_32_Dptr(arrayMSDptr(mant1,128/intDsize) mspop 2*32/intDsize,0);
      set_32_Dptr(arrayMSDptr(mant1,128/intDsize) mspop 3*32/intDsize,0);
      set_32_Dptr(arrayMSDptr(mant2,64/intDsize),manthi2);
      set_32_Dptr(arrayMSDptr(mant2,64/intDsize) mspop 32/intDsize,mantlo2);
      #else
      {var uintD* ptr;
       ptr = arrayLSDptr(mant1,128/intDsize);
       doconsttimes(64/intDsize, { lsprefnext(ptr) = 0; } );
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)mantlo1; mantlo1 = mantlo1>>intDsize; } );
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)manthi1; manthi1 = manthi1>>intDsize; } );
      }
      {var uintD* ptr;
       ptr = arrayLSDptr(mant2,64/intDsize);
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)mantlo2; mantlo2 = mantlo2>>intDsize; } );
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)manthi2; manthi2 = manthi2>>intDsize; } );
      }
      #endif
      #if (cl_word_size==64)
      var uint64 mantx;
      #endif
      #if (intDsize<=32)
      var uintL manthi;
      var uintL mantlo;
      #endif
      {CL_ALLOCA_STACK;
       var DS q;
       var DS r;
       UDS_divide(arrayMSDptr(mant1,128/intDsize),128/intDsize,arrayLSDptr(mant1,128/intDsize),
                  arrayMSDptr(mant2,64/intDsize),64/intDsize,arrayLSDptr(mant2,64/intDsize),
                  &q, &r
                 );
       // Es ist 2^53 <= q < 2^55, also q.len = ceiling(54/intDsize)=ceiling(55/intDsize),
       // und r=0 genau dann, wenn r.len=0.
       ASSERT(q.len==ceiling(54,intDsize))
       {var uintD* ptr = q.MSDptr;
        #if (intDsize==64)
        mantx = mspref(ptr,0);
        #else // (intDsize<=32)
        manthi = get_max32_Dptr(23,ptr);
        mantlo = get_32_Dptr(ptr mspop ceiling(23,intDsize));
        #endif
       }
       // q = 2^32*manthi+mantlo.
       #if (cl_word_size==64)
       #if (intDsize<=32)
       mantx = ((uint64)manthi<<32) | (uint64)mantlo;
       #endif
       if (mantx >= bit(DF_mant_len+2))
         // Quotient >=2^54 -> 2 Bits wegrunden
         { var uint64 rounding_bits = mantx & (bit(2)-1);
           exp1 += 1; // Exponenten incrementieren
           mantx = mantx >> 2;
           if ( (rounding_bits < bit(1)) // 00,01 werden abgerundet
                || ( (rounding_bits == bit(1)) // 10
                     && (r.len == 0) // und genau halbzahlig
                     && ((mantx & bit(0)) ==0) // -> round-to-even
              )    )
             // abrunden
             {}
             else
             // aufrunden
             { mantx += 1; }
         }
         else
         // Quotient <2^54 -> 1 Bit wegrunden
         { var uint64 rounding_bit = mantx & bit(0);
           mantx = mantx >> 1;
           if ( (rounding_bit == 0) // 0 wird abgerundet
                || ( (r.len == 0) // genau halbzahlig
                     && ((mantx & bit(0)) ==0) // -> round-to-even
              )    )
             // abrunden
             {}
             else
             // aufrunden
             { mantx += 1;
               if (mantx >= bit(DF_mant_len+1)) // rounding overflow?
                 { mantx = mantx>>1; exp1 = exp1+1; }
         }   }
       #else
       if (manthi >= bit(DF_mant_len-32+2))
         // Quotient >=2^54 -> 2 Bits wegrunden
         { var uintL rounding_bits = mantlo & (bit(2)-1);
           exp1 += 1; // Exponenten incrementieren
           mantlo = (mantlo >> 2) | (manthi << 30); manthi = manthi >> 2;
           if ( (rounding_bits < bit(1)) // 00,01 werden abgerundet
                || ( (rounding_bits == bit(1)) // 10
                     && (r.len == 0) // und genau halbzahlig
                     && ((mantlo & bit(0)) ==0) // -> round-to-even
              )    )
             // abrunden
             {}
             else
             // aufrunden
             { mantlo += 1; if (mantlo==0) { manthi += 1; } }
         }
         else
         // Quotient <2^54 -> 1 Bit wegrunden
         { var uintL rounding_bit = mantlo & bit(0);
           mantlo = (mantlo >> 1) | (manthi << 31); manthi = manthi >> 1;
           if ( (rounding_bit == 0) // 0 wird abgerundet
                || ( (r.len == 0) // genau halbzahlig
                     && ((mantlo & bit(0)) ==0) // -> round-to-even
              )    )
             // abrunden
             {}
             else
             // aufrunden
             { mantlo += 1;
               if (mantlo==0)
                 { manthi += 1;
                   if (manthi >= bit(DF_mant_len-32+1)) // rounding overflow?
                     { manthi = manthi>>1; exp1 = exp1+1; }
         }   }   }
       #endif
      }
      #if (cl_word_size==64)
      return encode_DF(sign1,exp1,mantx);
      #else
      return encode_DF(sign1,exp1,manthi,mantlo);
      #endif
#endif
}

}  // namespace cln
