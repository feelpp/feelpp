// binary operator +

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "base/cl_xmacros.h"

namespace cln {

const cl_SF operator+ (const cl_SF& x1, const cl_SF& x2)
{
// Methode (nach [Knuth, II, Seminumerical Algorithms, Abschnitt 4.2.1., S.200]):
// x1=0.0 -> Ergebnis x2.
// x2=0.0 -> Ergebnis x1.
// Falls e1<e2, vertausche x1 und x2.
// Also e1 >= e2.
// Falls e1 - e2 >= 16 + 3, Ergebnis x1.
// Schiebe beide Mantissen um 3 Bits nach links (Vorbereitung der Rundung:
//   Bei e1-e2=0,1 ist keine Rundung nötig, bei e1-e2>1 ist der Exponent des
//   Ergebnisses =e1-1, =e1 oder =e1+1. Brauche daher 1 Schutzbit und zwei
//   Rundungsbits: 00 exakt, 01 1.Hälfte, 10 exakte Mitte, 11 2.Hälfte.)
// Schiebe die Mantisse von x2 um e0-e1 Bits nach rechts. (Dabei die Rundung
// ausführen: Bit 0 ist das logische Oder der Bits 0,-1,-2,...)
// Falls x1,x2 selbes Vorzeichen haben: Addiere dieses zur Mantisse von x1.
// Falls x1,x2 verschiedenes Vorzeichen haben: Subtrahiere dieses von der
//   Mantisse von x1. <0 -> (Es war e1=e2) Vertausche die Vorzeichen, negiere.
//                    =0 -> Ergebnis 0.0
// Exponent ist e1.
// Normalisiere, fertig.
      // x1,x2 entpacken:
      var cl_signean sign1;
      var sintL exp1;
      var uintL mant1;
      var cl_signean sign2;
      var sintL exp2;
      var uintL mant2;
      SF_decode(x1, { return x2; }, sign1=,exp1=,mant1=);
      SF_decode(x2, { return x1; }, sign2=,exp2=,mant2=);
      var cl_uint max_x1_x2 = x1.word;
      if (exp1 < exp2)
        { max_x1_x2 = x2.word;
          swap(cl_signean, sign1,sign2);
          swap(sintL,      exp1 ,exp2 );
          swap(uintL,      mant1,mant2);
        }
      // Nun ist exp1>=exp2.
     {var uintL expdiff = exp1 - exp2; // Exponentendifferenz
      if (expdiff >= SF_mant_len+3) // >= 16+3 ?
        { return cl_SF_from_word(max_x1_x2); }
      mant1 = mant1 << 3; mant2 = mant2 << 3;
      // Nun 2^(SF_mant_len+3) <= mant1,mant2 < 2^(SF_mant_len+4).
      {var uintL mant2_last = mant2 & (bit(expdiff)-1); // letzte expdiff Bits von mant2
       mant2 = mant2 >> expdiff; if (!(mant2_last==0)) { mant2 |= bit(0); }
      }
      // mant2 = um expdiff Bits nach rechts geschobene und gerundete Mantisse
      // von x2.
      if ((x1.word ^ x2.word) & bit(SF_sign_shift))
        // verschiedene Vorzeichen -> Mantissen subtrahieren
        { if (mant1 > mant2) { mant1 = mant1 - mant2; goto norm_2; }
          if (mant1 == mant2) // Ergebnis 0 ?
            { return SF_0; }
          // negatives Subtraktionsergebnis
          mant1 = mant2 - mant1; sign1 = sign2; goto norm_2;
        }
        else
        // gleiche Vorzeichen -> Mantissen addieren
        { mant1 = mant1 + mant2; }
      // mant1 = Ergebnis-Mantisse >0, sign1 = Ergebnis-Vorzeichen,
      // exp1 = Ergebnis-Exponent.
      // Außerdem: Bei expdiff=0,1 sind die zwei letzten Bits von mant1 Null,
      // bei expdiff>=2 ist mant1 >= 2^(SF_mant_len+2).
      // Stets ist mant1 < 2^(SF_mant_len+5). (Daher werden die 2 Rundungsbits
      // nachher um höchstens eine Position nach links geschoben werden.)
      // [Knuth, S.201, leicht modifiziert:
      //   N1. m>=1 -> goto N4.
      //   N2. [Hier m<1] m>=1/2 -> goto N5.
      //       N3. m:=2*m, e:=e-1, goto N2.
      //   N4. [Hier 1<=m<2] m:=m/2, e:=e+1.
      //   N5. [Hier 1/2<=m<1] Runde m auf 17 Bits hinterm Komma.
      //       Falls hierdurch m=1 geworden, setze m:=m/2, e:=e+1.
      // ]
      // Bei uns ist m=mant1/2^(SF_mant_len+4),
      // ab Schritt N5 ist m=mant1/2^(SF_mant_len+1).
      norm_1: // [Knuth, S.201, Schritt N1]
      if (mant1 >= bit(SF_mant_len+4)) goto norm_4;
      norm_2: // [Knuth, S.201, Schritt N2]
              // Hier ist mant1 < 2^(SF_mant_len+4)
      if (mant1 >= bit(SF_mant_len+3)) goto norm_5;
      // [Knuth, S.201, Schritt N3]
      mant1 = mant1 << 1; exp1 = exp1-1; // Mantisse links schieben
      goto norm_2;
      norm_4: // [Knuth, S.201, Schritt N4]
              // Hier ist 2^(SF_mant_len+4) <= mant1 < 2^(SF_mant_len+5)
      exp1 = exp1+1;
      mant1 = (mant1>>1) | (mant1 & bit(0)); // Mantisse rechts schieben
      norm_5: // [Knuth, S.201, Schritt N5]
              // Hier ist 2^(SF_mant_len+3) <= mant1 < 2^(SF_mant_len+4)
      // Auf SF_mant_len echte Mantissenbits runden, d.h. rechte 3 Bits
      // wegrunden, und dabei mant1 um 3 Bits nach rechts schieben:
      {var uintL rounding_bits = mant1 & (bit(3)-1);
       mant1 = mant1 >> 3;
       if ( (rounding_bits < bit(2)) // 000,001,010,011 werden abgerundet
            || ( (rounding_bits == bit(2)) // 100 (genau halbzahlig)
                 && ((mant1 & bit(0)) ==0) // -> round-to-even
          )    )
         // abrunden
         {}
         else
         // aufrunden
         { mant1 = mant1+1;
           if (mant1 >= bit(SF_mant_len+1))
             // Bei Überlauf während der Rundung nochmals rechts schieben
             // (Runden ist hier überflüssig):
             { mant1 = mant1>>1; exp1 = exp1+1; } // Mantisse rechts schieben
         }
      }// Runden fertig
      return encode_SF(sign1,exp1,mant1);
     }
}

}  // namespace cln
