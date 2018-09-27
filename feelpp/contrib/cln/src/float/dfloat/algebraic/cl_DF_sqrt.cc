// sqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "float/cl_F.h"
#include "base/cl_low.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_DF sqrt (const cl_DF& x)
{
// Methode:
// x = 0.0 -> Ergebnis 0.0
// Ergebnis-Vorzeichen := positiv,
// Ergebnis-Exponent := ceiling(e/2),
// Ergebnis-Mantisse:
//   Bilde aus [1,m51,...,m0,(55 Nullbits)] bei geradem e,
//         aus [0,1,m51,...,m0,(54 Nullbits)] bei ungeradem e
//   die Ganzzahl-Wurzel, eine 54-Bit-Zahl mit einer führenden 1.
//   Runde das letzte Bit weg:
//     Bit 0 = 0 -> abrunden,
//     Bit 0 = 1 und Wurzel exakt -> round-to-even,
//     Bit 0 = 1 und Rest >0 -> aufrunden.
//   Dabei um ein Bit nach rechts schieben.
//   Bei Aufrundung auf 2^53 (rounding overflow) Mantisse um 1 Bit nach rechts
//     schieben und Exponent incrementieren.
#if (cl_word_size==64)
      // x entpacken:
      var sintL exp;
      var uint64 mantx;
      DF_decode(x, { return x; }, ,exp=,mantx=);
      // Um die 128-Bit-Ganzzahl-Wurzel ausnutzen zu können, fügen wir beim
      // Radikanden 74 bzw. 75 statt 54 bzw. 55 Nullbits an.
      if (exp & bit(0))
        // e ungerade
        { mantx = mantx << (63-(DF_mant_len+1)); exp = exp+1; }
        else
        // e gerade
        { mantx = mantx << (64-(DF_mant_len+1)); }
      exp = exp >> 1; // exp := exp/2
      var uintD mant [128/intDsize];
      #if (intDsize==64)
      arrayLSref(mant,128/intDsize,1) = mantx;
      arrayLSref(mant,128/intDsize,0) = 0;
      #else // (intDsize<=32)
      set_32_Dptr(arrayMSDptr(mant,128/intDsize),(uint32)(mantx>>32));
      set_32_Dptr(arrayMSDptr(mant,128/intDsize) mspop 32/intDsize,(uint32)mantx);
      set_32_Dptr(arrayMSDptr(mant,128/intDsize) mspop 2*32/intDsize,0);
      set_32_Dptr(arrayMSDptr(mant,128/intDsize) mspop 3*32/intDsize,0);
      #endif
      {CL_ALLOCA_STACK;
       var DS wurzel;
       var bool exactp;
       UDS_sqrt(arrayMSDptr(mant,128/intDsize),128/intDsize,arrayLSDptr(mant,128/intDsize), &wurzel, exactp=);
       // wurzel = isqrt(2^74_75 * mant), eine 64-Bit-Zahl.
       mantx = get_64_Dptr(wurzel.MSDptr);
       // Die hinteren 63-DF_mant_len Bits wegrunden:
       if ( ((mantx & bit(62-DF_mant_len)) ==0) // Bit 10 =0 -> abrunden
            || ( ((mantx & (bit(62-DF_mant_len)-1)) ==0) // Bit 10 =1 und Bits 9..0 >0 -> aufrunden
                 && exactp                   // Bit 10 =1 und Bits 9..0 =0, aber Rest -> aufrunden
                 // round-to-even, je nach Bit 11 :
                 && ((mantx & bit(63-DF_mant_len)) ==0)
          )    )
         // abrunden
         { mantx = mantx >> (63-DF_mant_len); }
         else
         // aufrunden
         { mantx = mantx >> (63-DF_mant_len);
           mantx += 1;
           if (mantx >= bit(DF_mant_len+1)) // rounding overflow?
             { mantx = mantx>>1; exp = exp+1; }
         }
      }
      return encode_DF(0,exp,mantx);
#else
      // x entpacken:
      var sintL exp;
      var uint32 manthi;
      var uint32 mantlo;
      DF_decode2(x, { return x; }, ,exp=,manthi=,mantlo=);
      // Um die 128-Bit-Ganzzahl-Wurzel ausnutzen zu können, fügen wir beim
      // Radikanden 74 bzw. 75 statt 54 bzw. 55 Nullbits an.
      if (exp & bit(0))
        // e ungerade
        { manthi = (manthi << (63-(DF_mant_len+1))) | (mantlo >> ((DF_mant_len+1)-31));
          mantlo = mantlo << (63-(DF_mant_len+1));
          exp = exp+1;
        }
        else
        // e gerade
        { manthi = (manthi << (64-(DF_mant_len+1))) | (mantlo >> ((DF_mant_len+1)-32));
          mantlo = mantlo << (64-(DF_mant_len+1));
        }
      exp = exp >> 1; // exp := exp/2
      var uintD mant [128/intDsize];
      #if (intDsize==32) || (intDsize==16) || (intDsize==8)
      set_32_Dptr(arrayMSDptr(mant,128/intDsize),manthi);
      set_32_Dptr(arrayMSDptr(mant,128/intDsize) mspop 32/intDsize,mantlo);
      set_32_Dptr(arrayMSDptr(mant,128/intDsize) mspop 2*32/intDsize,0);
      set_32_Dptr(arrayMSDptr(mant,128/intDsize) mspop 3*32/intDsize,0);
      #else
      {var uintD* ptr;
       ptr = arrayLSDptr(mant,128/intDsize);
       doconsttimes(64/intDsize, { lsprefnext(ptr) = 0; } );
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)mantlo; mantlo = mantlo>>intDsize; } );
       doconsttimes(32/intDsize, { lsprefnext(ptr) = (uintD)manthi; manthi = manthi>>intDsize; } );
      }
      #endif
      {CL_ALLOCA_STACK;
       var DS wurzel;
       var bool exactp;
       UDS_sqrt(arrayMSDptr(mant,128/intDsize),128/intDsize,arrayLSDptr(mant,128/intDsize), &wurzel, exactp=);
       // wurzel = isqrt(2^74_75 * mant), eine 64-Bit-Zahl.
       {var uintD* ptr = wurzel.MSDptr;
        manthi = get_32_Dptr(ptr); mantlo = get_32_Dptr(ptr mspop 32/intDsize);
       }
       // Die hinteren 63-DF_mant_len Bits wegrunden:
       if ( ((mantlo & bit(62-DF_mant_len)) ==0) // Bit 10 =0 -> abrunden
            || ( ((mantlo & (bit(62-DF_mant_len)-1)) ==0) // Bit 10 =1 und Bits 9..0 >0 -> aufrunden
                 && exactp                   // Bit 10 =1 und Bits 9..0 =0, aber Rest -> aufrunden
                 // round-to-even, je nach Bit 11 :
                 && ((mantlo & bit(63-DF_mant_len)) ==0)
          )    )
         // abrunden
         { mantlo = (mantlo >> (63-DF_mant_len)) | (manthi << (DF_mant_len-32+1));
           manthi = manthi >> (63-DF_mant_len);
         }
         else
         // aufrunden
         { mantlo = (mantlo >> (63-DF_mant_len)) | (manthi << (DF_mant_len-32+1));
           manthi = manthi >> (63-DF_mant_len);
           mantlo += 1;
           if (mantlo==0)
             { manthi += 1;
               if (manthi >= bit(DF_mant_len-32+1)) // rounding overflow?
                 { manthi = manthi>>1; exp = exp+1; }
         }   }
      }
      return encode_DF(0,exp,manthi,mantlo);
#endif
}

}  // namespace cln
