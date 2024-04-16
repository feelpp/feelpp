// futruncate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/sfloat/cl_SF.h"


// Implementation.

namespace cln {

const cl_SF futruncate (const cl_SF& x)
{
// Methode:
// x = 0.0 -> Ergebnis 0.0
// e<=0 -> Ergebnis 1.0 oder -1.0, je nach Vorzeichen von x.
// 1<=e<=16 -> Greife die letzten (17-e) Bits von x heraus.
//             Sind sie alle =0 -> Ergebnis x.
//             Sonst setze sie alle und erhöhe dann die letzte Stelle um 1.
//             Kein Überlauf der 16 Bit -> fertig.
//             Sonst (Ergebnis eine Zweierpotenz): Mantisse := .1000...000,
//               e:=e+1. (Test auf Überlauf wegen e<=17 überflüssig)
// e>=17 -> Ergebnis x.
      var uintL uexp = SF_uexp(x); // e + SF_exp_mid
      if (uexp==0) // 0.0 ?
        { return x; }
      if (uexp <= SF_exp_mid) // e<=0 ?
        { // Exponent auf 1, Mantisse auf .1000...000 setzen.
          return cl_SF_from_word(
            (x.word & ~(((bit(SF_exp_len)-1)<<SF_exp_shift)
                        + ((bit(SF_mant_len)-1)<<SF_mant_shift)))
            | ((cl_uint)(SF_exp_mid+1) << SF_exp_shift)
            | ((cl_uint)0 << SF_mant_shift)
            );
        }
        else
        { if (uexp > SF_exp_mid+SF_mant_len) // e > 16 ?
            { return x; }
            else
            { var cl_uint mask = // Bitmaske: Bits 16-e..0 gesetzt, alle anderen gelöscht
                bit(SF_mant_len+SF_mant_shift + 1+SF_exp_mid-uexp) - bit(SF_mant_shift);
              if ((x.word & mask)==0) // alle diese Bits =0 ?
                { return x; }
              return cl_SF_from_word(
                (x.word | mask) // alle diese Bits setzen
                + bit(SF_mant_shift) // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                );
        }   }
}

}  // namespace cln
