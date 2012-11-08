// futruncate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/ffloat/cl_FF.h"


// Implementation.

namespace cln {

const cl_FF futruncate (const cl_FF& x)
{
// Methode:
// x = 0.0 -> Ergebnis 0.0
// e<=0 -> Ergebnis 1.0 oder -1.0, je nach Vorzeichen von x.
// 1<=e<=23 -> Greife die letzten (24-e) Bits von x heraus.
//             Sind sie alle =0 -> Ergebnis x.
//             Sonst setze sie alle und erhöhe dann die letzte Stelle um 1.
//             Kein Überlauf der 23 Bit -> fertig.
//             Sonst (Ergebnis eine Zweierpotenz): Mantisse := .1000...000,
//               e:=e+1. (Test auf Überlauf wegen e<=24 überflüssig)
// e>=24 -> Ergebnis x.
      var ffloat x_ = cl_ffloat_value(x);
      var uintL uexp = FF_uexp(x_); // e + FF_exp_mid
      if (uexp==0) // 0.0 ?
        { return x; }
      if (uexp <= FF_exp_mid) // e<=0 ?
        { // Exponent auf 1, Mantisse auf .1000...000 setzen.
          return ((x_ & bit(31))==0 ? cl_FF_1 : cl_FF_minus1);
        }
        else
        { if (uexp > FF_exp_mid+FF_mant_len) // e > 23 ?
            { return x; }
            else
            { var uint32 mask = // Bitmaske: Bits 23-e..0 gesetzt, alle anderen gelöscht
                bit(FF_mant_len+1+FF_exp_mid-uexp)-1;
              if ((x_ & mask)==0) // alle diese Bits =0 ?
                { return x; }
              return allocate_ffloat
                ((x_ | mask) // alle diese Bits setzen
                 + 1 // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                );
        }   }
}

}  // namespace cln
