// fround().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {

const cl_DF fround (const cl_DF& x)
{
// Methode:
// x = 0.0 oder e<0 -> Ergebnis 0.0
// 0<=e<=52 -> letzte (53-e) Bits der Mantisse wegrunden,
//             Exponent und Vorzeichen beibehalten.
// e>52 -> Ergebnis x
#if (cl_word_size==64)
      var dfloat x_ = TheDfloat(x)->dfloat_value;
      var uintL uexp = DF_uexp(x_); // e + DF_exp_mid
      if (uexp < DF_exp_mid) // x = 0.0 oder e<0 ?
        { return cl_DF_0; }
        else
        { if (uexp > DF_exp_mid+DF_mant_len) // e > 52 ?
            { return x; }
            else
            if (uexp > DF_exp_mid+1) // e>1 ?
              { var uint64 bitmask = // Bitmaske: Bit 52-e gesetzt, alle anderen gelöscht
                  bit(DF_mant_len+DF_exp_mid-uexp);
                var uint64 mask = // Bitmaske: Bits 51-e..0 gesetzt, alle anderen gelöscht
                  bitmask-1;
                if ( ((x_ & bitmask) ==0) // Bit 52-e =0 -> abrunden
                     || ( ((x_ & mask) ==0) // Bit 52-e =1 und Bits 51-e..0 >0 -> aufrunden
                          // round-to-even, je nach Bit 53-e :
                          && ((x_ & (bitmask<<1)) ==0)
                   )    )
                  // abrunden
                  { mask |= bitmask; // Bitmaske: Bits 52-e..0 gesetzt, alle anderen gelöscht
                    return allocate_dfloat( x_ & ~mask );
                  }
                  else
                  // aufrunden
                  { return allocate_dfloat
                      ((x_ | mask) // alle diese Bits 51-e..0 setzen (Bit 52-e schon gesetzt)
                       + 1 // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                      );
                  }
              }
            elif (uexp == DF_exp_mid+1) // e=1 ?
              // Wie bei 1 < e <= 52, nur daß Bit 53-e stets gesetzt ist.
              { if ((x_ & bit(DF_mant_len-1)) ==0) // Bit 52-e =0 -> abrunden
                  // abrunden
                  { return allocate_dfloat( x_ & ~(bit(DF_mant_len)-1) ); }
                  else
                  // aufrunden
                  { return allocate_dfloat
                      ((x_ | (bit(DF_mant_len)-1)) // alle diese Bits 52-e..0 setzen
                       + 1 // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                      );
                  }
              }
            else // e=0 ?
              // Wie bei 1 < e <= 52, nur daß Bit 52-e stets gesetzt
              // und Bit 53-e stets gelöscht ist.
              { if ((x_ & (bit(DF_mant_len)-1)) ==0)
                  // abrunden von +-0.5 zu 0.0
                  { return cl_DF_0; }
                  else
                  // aufrunden
                  { return allocate_dfloat
                      ((x_ | (bit(DF_mant_len)-1)) // alle Bits 51-e..0 setzen
                       + 1 // letzte Stelle erhöhen, dabei Exponenten incrementieren
                      );
              }   }
        }
#else
      var uint32 semhi = TheDfloat(x)->dfloat_value.semhi;
      var uint32 mlo = TheDfloat(x)->dfloat_value.mlo;
      var uintL uexp = DF_uexp(semhi); // e + DF_exp_mid
      if (uexp < DF_exp_mid) // x = 0.0 oder e<0 ?
        { return cl_DF_0; }
        else
        { if (uexp > DF_exp_mid+DF_mant_len) // e > 52 ?
            { return x; }
            else
            if (uexp > DF_exp_mid+1) // e>1 ?
              { if (uexp > DF_exp_mid+DF_mant_len-32) // e > 20 ?
                  { var uint32 bitmask = // Bitmaske: Bit 52-e gesetzt, alle anderen gelöscht
                      bit(DF_mant_len+DF_exp_mid-uexp);
                    var uint32 mask = // Bitmaske: Bits 51-e..0 gesetzt, alle anderen gelöscht
                      bitmask-1;
                    if ( ((mlo & bitmask) ==0) // Bit 52-e =0 -> abrunden
                         || ( ((mlo & mask) ==0) // Bit 52-e =1 und Bits 51-e..0 >0 -> aufrunden
                              // round-to-even, je nach Bit 53-e :
                              && ( ((bitmask<<1) == 0) // e=21 ?
                                    ? ((semhi & bit(0)) ==0)
                                    : ((mlo & (bitmask<<1)) ==0)
                       )    )    )
                      // abrunden
                      { mask |= bitmask; // Bitmaske: Bits 52-e..0 gesetzt, alle anderen gelöscht
                        return allocate_dfloat(semhi, mlo & ~mask );
                      }
                      else
                      // aufrunden
                      { mlo = (mlo | mask) // alle diese Bits 51-e..0 setzen (Bit 52-e schon gesetzt)
                              + 1; // letzte Stelle erhöhen,
                        if (mlo==0) { semhi += 1; } // dabei evtl. Exponenten incrementieren
                        return allocate_dfloat(semhi,mlo);
                      }
                  }
                  else
                  { var uint32 bitmask = // Bitmaske: Bit 20-e gesetzt, alle anderen gelöscht
                      bit(DF_mant_len+DF_exp_mid-32-uexp);
                    var uint32 mask = // Bitmaske: Bits 19-e..0 gesetzt, alle anderen gelöscht
                      bitmask-1;
                    if ( ((semhi & bitmask) ==0) // Bit 52-e =0 -> abrunden
                         || ( (mlo==0) && ((semhi & mask) ==0) // Bit 52-e =1 und Bits 51-e..0 >0 -> aufrunden
                              // round-to-even, je nach Bit 53-e :
                              && ((semhi & (bitmask<<1)) ==0)
                       )    )
                      // abrunden
                      { mask |= bitmask; // Bitmaske: Bits 20-e..0 gesetzt, alle anderen gelöscht
                        return allocate_dfloat( semhi & ~mask, 0 );
                      }
                      else
                      // aufrunden
                      { return allocate_dfloat
                          ((semhi | mask) // alle diese Bits 19-e..0 setzen (Bit 20-e schon gesetzt)
                           + 1, // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                           0
                          );
                      }
                  }
              }
            elif (uexp == DF_exp_mid+1) // e=1 ?
              // Wie bei 1 < e <= 20, nur daß Bit 53-e stets gesetzt ist.
              { if ((semhi & bit(DF_mant_len-32-1)) ==0) // Bit 52-e =0 -> abrunden
                  // abrunden
                  { return allocate_dfloat( semhi & ~(bit(DF_mant_len-32)-1) , 0 ); }
                  else
                  // aufrunden
                  { return allocate_dfloat
                      ((semhi | (bit(DF_mant_len-32)-1)) // alle diese Bits 52-e..0 setzen
                       + 1, // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                       0
                      );
                  }
              }
            else // e=0 ?
              // Wie bei 1 < e <= 20, nur daß Bit 52-e stets gesetzt
              // und Bit 53-e stets gelöscht ist.
              { if ((mlo==0) && ((semhi & (bit(DF_mant_len-32)-1)) ==0))
                  // abrunden von +-0.5 zu 0.0
                  { return cl_DF_0; }
                  else
                  // aufrunden
                  { return allocate_dfloat
                      ((semhi | (bit(DF_mant_len-32)-1)) // alle Bits 51-e..0 setzen
                       + 1, // letzte Stelle erhöhen, dabei Exponenten incrementieren
                       0
                      );
              }   }
        }
#endif
}

}  // namespace cln
