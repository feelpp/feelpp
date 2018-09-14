// binary operator *

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "float/cl_F.h"
#include "base/cl_low.h"
#include "base/digitseq/cl_DS.h"

#include "base/cl_inline.h"
#include "float/dfloat/elem/cl_DF_zerop.cc"

namespace cln {


const cl_DF operator* (const cl_DF& x1, const cl_DF& x2)
{
// Methode:
// Falls x1=0.0 oder x2=0.0 -> Ergebnis 0.0
// Sonst: Ergebnis-Vorzeichen = VZ von x1 xor VZ von x2.
//        Ergebnis-Exponent = Summe der Exponenten von x1 und x2.
//        Ergebnis-Mantisse = Produkt der Mantissen von x1 und x2, gerundet:
//          2^-53 * mant1  *  2^-53 * mant2  =  2^-106 * (mant1*mant2),
//          die Klammer ist >=2^104, <=(2^53-1)^2<2^106 .
//          Falls die Klammer >=2^105 ist, um 53 Bit nach rechts schieben und
//            runden: Falls Bit 52 Null, abrunden; falls Bit 52 Eins und
//            Bits 51..0 alle Null, round-to-even; sonst aufrunden.
//          Falls die Klammer <2^105 ist, um 52 Bit nach rechts schieben und
//            runden: Falls Bit 51 Null, abrunden; falls Bit 51 Eins und
//            Bits 50..0 alle Null, round-to-even; sonst aufrunden. Nach
//            Aufrunden: Falls =2^53, um 1 Bit nach rechts schieben. Sonst
//            Exponenten um 1 erniedrigen.
#ifdef FAST_DOUBLE
      double_to_DF(DF_to_double(x1) * DF_to_double(x2), return ,
                   TRUE, TRUE, // Overflow und subnormale Zahl abfangen
                   !(zerop_inline(x1) || zerop_inline(x2)), // ein Ergebnis +/- 0.0
                               // ist genau dann in Wirklichkeit ein Underflow
                   FALSE, FALSE // keine Singularität, kein NaN als Ergebnis möglich
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
      DF_decode(x1, { return x1; }, sign1=,exp1=,mantx1=);
      #if (intDsize<=32)
      manthi1 = high32(mantx1); mantlo1 = low32(mantx1);
      #endif
      var uint64 mantx2;
      DF_decode(x2, { return x2; }, sign2=,exp2=,mantx2=);
      #if (intDsize<=32)
      manthi2 = high32(mantx2); mantlo2 = low32(mantx2);
      #endif
      #else
      DF_decode2(x1, { return x1; }, sign1=,exp1=,manthi1=,mantlo1=);
      DF_decode2(x2, { return x2; }, sign2=,exp2=,manthi2=,mantlo2=);
      #endif
      exp1 = exp1 + exp2; // Summe der Exponenten
      sign1 = sign1 ^ sign2; // Ergebnis-Vorzeichen
      // Mantissen mant1 und mant2 multiplizieren (64x64-Bit-Multiplikation):
      var uintD mant1 [64/intDsize];
      var uintD mant2 [64/intDsize];
      var uintD mant [128/intDsize];
      #if (intDsize==64)
      arrayLSref(mant1,64/intDsize,0) = mantx1;
      arrayLSref(mant2,64/intDsize,0) = mantx2;
      #elif (intDsize==32) || (intDsize==16) || (intDsize==8)
      set_32_Dptr(arrayMSDptr(mant1,64/intDsize),manthi1);
      set_32_Dptr(arrayMSDptr(mant1,64/intDsize) mspop 32/intDsize,mantlo1);
      set_32_Dptr(arrayMSDptr(mant2,64/intDsize),manthi2);
      set_32_Dptr(arrayMSDptr(mant2,64/intDsize) mspop 32/intDsize,mantlo2);
      #else
      {var uintD* ptr;
       ptr = arrayLSDptr(mant1,64/intDsize);
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)mantlo1; mantlo1 = mantlo1>>intDsize; } );
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)manthi1; manthi1 = manthi1>>intDsize; } );
      }
      {var uintD* ptr;
       ptr = arrayLSDptr(mant2,64/intDsize);
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)mantlo2; mantlo2 = mantlo2>>intDsize; } );
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)manthi2; manthi2 = manthi2>>intDsize; } );
      }
      #endif
      cl_UDS_mul(arrayLSDptr(mant1,64/intDsize),64/intDsize,
                 arrayLSDptr(mant2,64/intDsize),64/intDsize,
                 arrayLSDptr(mant,128/intDsize)
                );
      {
        #if (cl_word_size==64)
        var uint64 manterg;
        #else
        var uintL manthi;
        var uintL mantlo;
        #endif
        // Produkt mant = mant1 * mant2 ist >= 2^104, < 2^106. Bit 105 abtesten:
        #define mant_bit(k)  (arrayLSref(mant,128/intDsize,floor(k,intDsize)) & bit((k)%intDsize))
        if (mant_bit(2*DF_mant_len+1))
          // mant>=2^(2*DF_mant_len+1), um DF_mant_len+1 Bits nach rechts schieben:
          { // Bits 105..53 holen:
            #if (cl_word_size==64) && (intDsize==64)
              manterg = ((uint64)arrayLSref(mant,2,1) << 11) | ((uint64)arrayLSref(mant,2,0) >> 53); // Bits 116..53
              #define mantrest() (arrayLSref(mant,2,0) & (bit(53)-1))
            #elif (cl_word_size==64) && (intDsize==32)
              manterg = ((uint64)arrayLSref(mant,4,3) << 43) | ((uint64)arrayLSref(mant,4,2) << 11) | ((uint64)arrayLSref(mant,4,1) >> 21); // Bits 116..53
              #define mantrest() ((arrayLSref(mant,4,1) & (bit(21)-1)) || arrayLSref(mant,4,0))
            #elif (intDsize==32)
              manthi = ((uint32)arrayLSref(mant,4,3) << 11) | ((uint32)arrayLSref(mant,4,2) >> 21); // Bits 116..85
              mantlo = ((uint32)arrayLSref(mant,4,2) << 11) | ((uint32)arrayLSref(mant,4,1) >> 21); // Bits 84..53
              #define mantrest() ((arrayLSref(mant,4,1) & (bit(21)-1)) || arrayLSref(mant,4,0))
            #elif (intDsize==16)
              manthi = ((uint32)arrayLSref(mant,8,7) << 27) | ((uint32)arrayLSref(mant,8,6) << 11) | ((uint32)arrayLSref(mant,8,5) >> 5); // Bits 116..85
              mantlo = ((uint32)arrayLSref(mant,8,5) << 27) | ((uint32)arrayLSref(mant,8,4) << 11) | ((uint32)arrayLSref(mant,8,3) >> 5); // Bits 84..53
              #define mantrest() ((arrayLSref(mant,8,3) & (bit(5)-1)) || arrayLSref(mant,8,2) || arrayLSref(mant,8,1) || arrayLSref(mant,8,0))
            #elif (intDsize==8)
              manthi = ((uint32)arrayLSref(mant,16,14) << 27) | ((uint32)arrayLSref(mant,16,13) << 19) | ((uint32)arrayLSref(mant,16,12) << 11) | ((uint32)arrayLSref(mant,16,11) << 3) | ((uint32)arrayLSref(mant,16,10) >> 5); // Bits 116..85
              mantlo = ((uint32)arrayLSref(mant,16,10) << 27) | ((uint32)arrayLSref(mant,16,9) << 19) | ((uint32)arrayLSref(mant,16,8) << 11) | ((uint32)arrayLSref(mant,16,7) << 3) | ((uint32)arrayLSref(mant,16,6) >> 5); // Bits 84..53
              #define mantrest() ((arrayLSref(mant,16,6) & (bit(5)-1)) || arrayLSref(mant,16,5) || arrayLSref(mant,16,4) || arrayLSref(mant,16,3) || arrayLSref(mant,16,2) || arrayLSref(mant,16,1) || arrayLSref(mant,16,0))
            #endif
            if ( (mant_bit(DF_mant_len) ==0) // Bit DF_mant_len =0 -> abrunden
                 || ( !mantrest() // Bit DF_mant_len =1 und Bits DF_mant_len-1..0 >0 -> aufrunden
                      // round-to-even, je nach Bit DF_mant_len+1 :
                      && (mant_bit(DF_mant_len+1) ==0)
               )    )
              // abrunden
              goto ab;
              else
              // aufrunden
              goto auf;
            #undef mantrest
          }
          else
          // mant<2^(2*DF_mant_len+1), um DF_mant_len Bits nach rechts schieben:
          { exp1 = exp1-1; // Exponenten decrementieren
            // Bits 104..52 holen:
            #if (cl_word_size==64) && (intDsize==64)
              manterg = ((uint64)arrayLSref(mant,2,1) << 12) | ((uint64)arrayLSref(mant,2,0) >> 52); // Bits 115..52
              #define mantrest() (arrayLSref(mant,2,0) & (bit(52)-1))
            #elif (cl_word_size==64) && (intDsize==32)
              manterg = ((uint64)arrayLSref(mant,4,3) << 44) | ((uint64)arrayLSref(mant,4,2) << 12) | ((uint64)arrayLSref(mant,4,1) >> 20); // Bits 115..52
              #define mantrest() ((arrayLSref(mant,4,1) & (bit(20)-1)) || arrayLSref(mant,4,0))
            #elif (intDsize==32)
              manthi = ((uint32)arrayLSref(mant,4,3) << 12) | ((uint32)arrayLSref(mant,4,2) >> 20); // Bits 115..84
              mantlo = ((uint32)arrayLSref(mant,4,2) << 12) | ((uint32)arrayLSref(mant,4,1) >> 20); // Bits 83..52
              #define mantrest() ((arrayLSref(mant,4,1) & (bit(20)-1)) || arrayLSref(mant,4,0))
            #elif (intDsize==16)
              manthi = // ((uint32)arrayLSref(mant,8,7) << 28) | ((uint32)arrayLSref(mant,8,6) << 12) | ((uint32)arrayLSref(mant,8,5) >> 4); // Bits 115..84
              mantlo = // ((uint32)arrayLSref(mant,8,5) << 28) | ((uint32)arrayLSref(mant,8,4) << 12) | ((uint32)arrayLSref(mant,8,3) >> 4); // Bits 83..52
              #define mantrest() ((arrayLSref(mant,8,3) & (bit(4)-1)) || arrayLSref(mant,8,2) || arrayLSref(mant,8,1) || arrayLSref(mant,8,0))
            #elif (intDsize==8)
              manthi = ((uint32)arrayLSref(mant,16,14) << 28) | ((uint32)arrayLSref(mant,16,13) << 20) | ((uint32)arrayLSref(mant,16,12) << 12) | ((uint32)arrayLSref(mant,16,11) << 4) | ((uint32)arrayLSref(mant,16,10) >> 4); // Bits 115..84
              mantlo = ((uint32)arrayLSref(mant,16,10) << 28) | ((uint32)arrayLSref(mant,16,9) << 20) | ((uint32)arrayLSref(mant,16,8) << 12) | ((uint32)arrayLSref(mant,16,7) << 4) | ((uint32)arrayLSref(mant,16,6) >> 4); // Bits 83..52
              #define mantrest() ((arrayLSref(mant,16,6) & (bit(4)-1)) || arrayLSref(mant,16,5) || arrayLSref(mant,16,4) || arrayLSref(mant,16,3) || arrayLSref(mant,16,2) || arrayLSref(mant,16,1) || arrayLSref(mant,16,0))
            #endif
            if ( (mant_bit(DF_mant_len-1) ==0) // Bit DF_mant_len-1 =0 -> abrunden
                 || ( !mantrest() // Bit DF_mant_len-1 =1 und Bits DF_mant_len-2..0 >0 -> aufrunden
                      // round-to-even, je nach Bit DF_mant_len :
                      && (mant_bit(DF_mant_len) ==0)
               )    )
              // abrunden
              goto ab;
              else
              // aufrunden
              goto auf;
            #undef mantrest
          }
        #undef mant_bit
        auf:
        #if (cl_word_size==64)
        manterg = manterg+1;
        // Hier ist 2^DF_mant_len <= manterg <= 2^(DF_mant_len+1)
        if (manterg >= bit(DF_mant_len+1)) // rounding overflow?
          { manterg = manterg>>1; exp1 = exp1+1; } // Shift nach rechts
        #else
        mantlo = mantlo+1;
        if (mantlo==0)
          { manthi = manthi+1;
            // Hier ist 2^(DF_mant_len-32) <= manthi <= 2^(DF_mant_len-32+1)
            if (manthi >= bit(DF_mant_len-32+1)) // rounding overflow?
              { manthi = manthi>>1; exp1 = exp1+1; } // Shift nach rechts
          }
        #endif
        ab:
        // Runden fertig, 2^DF_mant_len <= manterg < 2^(DF_mant_len+1)
        #if (cl_word_size==64)
        return encode_DF(sign1,exp1,manterg);
        #else
        return encode_DF(sign1,exp1,manthi,mantlo);
        #endif
      }
#endif
}

}  // namespace cln
