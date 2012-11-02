// futruncate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/dfloat/cl_DF.h"


// Implementation.

namespace cln {

const cl_DF futruncate (const cl_DF& x)
{
// Methode:
// x = 0.0 -> Ergebnis 0.0
// e<=0 -> Ergebnis 1.0 oder -1.0, je nach Vorzeichen von x.
// 1<=e<=52 -> Greife die letzten (53-e) Bits von x heraus.
//             Sind sie alle =0 -> Ergebnis x.
//             Sonst setze sie alle und erhöhe dann die letzte Stelle um 1.
//             Kein Überlauf der 52 Bit -> fertig.
//             Sonst (Ergebnis eine Zweierpotenz): Mantisse := .1000...000,
//               e:=e+1. (Test auf Überlauf wegen e<=53 überflüssig)
// e>=53 -> Ergebnis x.
#if (cl_word_size==64)
      var dfloat x_ = TheDfloat(x)->dfloat_value;
      var uintL uexp = DF_uexp(x_); // e + DF_exp_mid
      if (uexp==0) // 0.0 ?
        { return x; }
      if (uexp <= DF_exp_mid) // e<=0 ?
        { // Exponent auf 1, Mantisse auf .1000...000 setzen.
          return ((x_ & bit(63))==0 ? cl_DF_1 : cl_DF_minus1);
        }
        else
        { if (uexp > DF_exp_mid+DF_mant_len) // e > 52 ?
            { return x; }
            else
            { var uint64 mask = // Bitmaske: Bits 52-e..0 gesetzt, alle anderen gelöscht
                bit(DF_mant_len+1+DF_exp_mid-uexp)-1;
              if ((x_ & mask)==0) // alle diese Bits =0 ?
                { return x; }
              return allocate_dfloat
                ((x_ | mask) // alle diese Bits setzen
                 + 1 // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                );
        }   }
#else
      var uint32 semhi = TheDfloat(x)->dfloat_value.semhi;
      var uint32 mlo = TheDfloat(x)->dfloat_value.mlo;
      var uintL uexp = DF_uexp(semhi); // e + DF_exp_mid
      if (uexp==0) // 0.0 ?
        { return x; }
      if (uexp <= DF_exp_mid) // e<=0 ?
        { // Exponent auf 1, Mantisse auf .1000...000 setzen.
          return ((semhi & bit(31))==0 ? cl_DF_1 : cl_DF_minus1);
        }
        else
        { if (uexp > DF_exp_mid+DF_mant_len) // e > 52 ?
            { return x; }
            else
            if (uexp > DF_exp_mid+DF_mant_len+1-32) // e > 21 ?
              { var uint32 mask = // Bitmaske: Bits 52-e..0 gesetzt, alle anderen gelöscht
                  bit(DF_mant_len+1+DF_exp_mid-uexp)-1;
                if ((mlo & mask)==0) // alle diese Bits =0 ?
                  { return x; }
                mlo = (mlo | mask) // alle diese Bits setzen
                      + 1; // letzte Stelle erhöhen,
                if (mlo==0) { semhi += 1; } // dabei evtl. Exponenten incrementieren
                return allocate_dfloat(semhi,mlo);
              }
              else
              { var uint32 mask = // Bitmaske: Bits 20-e..0 gesetzt, alle anderen gelöscht
                  bit(DF_mant_len+1+DF_exp_mid-32-uexp)-1;
                if ((mlo==0) && ((semhi & mask)==0)) // alle diese Bits und mlo =0 ?
                  { return x; }
                return allocate_dfloat
                  ((semhi | mask) // alle diese Bits setzen
                   + 1, // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                   0
                  );
        }     }
#endif
}

}  // namespace cln
