// LF_LF_plus_LF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"
#include "float/cl_F.h"
#include "base/cl_xmacros.h"

namespace cln {

const cl_LF LF_LF_plus_LF (const cl_LF& arg1, const cl_LF& arg2)
{
// Methode (nach [Knuth, II, Seminumerical Algorithms, Abschnitt 4.2.1., S.200]):
// Falls e1<e2, vertausche x1 und x2.
// Also e1 >= e2.
// Falls e2=0, also x2=0.0, Ergebnis x1.
// Falls e1 - e2 >= 16n+2, Ergebnis x1.
// Erweitere die Mantissen rechts um 3 Bits (Bit -1 als Schutzbit, Bits -2,-3
//   als Rundungsbits: 00 exakt, 01 1.Hälfte, 10 exakte Mitte, 11 2.Hälfte.)
// Schiebe die Mantisse von x2 um e0-e1 Bits nach rechts. (Dabei die Rundung
// ausführen: Bit -3 ist das logische Oder der Bits -3,-4,-5,...)
// Falls x1,x2 selbes Vorzeichen haben: Addiere dieses zur Mantisse von x1.
// Falls x1,x2 verschiedenes Vorzeichen haben: Subtrahiere dieses von der
//   Mantisse von x1. <0 -> (Es war e1=e2) Vertausche die Vorzeichen, negiere.
//                    =0 -> Ergebnis 0.0
// Exponent ist e1.
// Normalisiere, fertig.
      var cl_LF x1 = arg1;
      var cl_LF x2 = arg2;
      var uintE uexp1 = TheLfloat(arg1)->expo;
      var uintE uexp2 = TheLfloat(arg2)->expo;
      if (uexp1 < uexp2)
        // x1 und x2 vertauschen
        { x1 = arg2; x2 = arg1; swap(uintE, uexp1,uexp2); }
      // uexp1 >= uexp2
      if (uexp2==0) { return x1; } // x2=0.0 -> x1 als Ergebnis
      var uintC len = TheLfloat(x1)->len; // Länge n von x1 und x2
      var uintE expdiff = uexp1-uexp2; // e1-e2
      if ((expdiff == 0) && (TheLfloat(x1)->sign != TheLfloat(x2)->sign))
        // verschiedene Vorzeichen, aber gleicher Exponent
        { // Vorzeichen des Ergebnisses festlegen:
          var cl_signean erg = // Mantissen (je len Digits) vergleichen
            compare_loop_msp(arrayMSDptr(TheLfloat(x1)->data,len),arrayMSDptr(TheLfloat(x2)->data,len),len);
          if (erg==0) // Mantissen gleich
            { return encode_LF0(len); } // Ergebnis 0.0
          if (erg<0) // |x1| < |x2|
            // x1 und x2 vertauschen, expdiff bleibt =0
            { x1.pointer = arg2.pointer; x2.pointer = arg1.pointer;
              swap(uintE, uexp1,uexp2);
            }
        }
      if (expdiff >= intDsize * len + 2) // e1-e2 >= 16n+2 ?
        { return x1; } // ja -> x1 als Ergebnis
      // neues Long-Float allozieren:
      var Lfloat y = allocate_lfloat(len,uexp1,TheLfloat(x1)->sign);
      var uintL i = floor(expdiff,intDsize); // e1-e2 div 16 (>=0, <=n)
      var uintL j = expdiff % intDsize; // e1-e2 mod 16 (>=0, <16)
      // Mantisse von x2 muß um intDsize*i+j Bits nach rechts geschoben werden.
      var uintC x2_len = len - i; // n-i Digits von x2 gebraucht
      // x2_len Digits um j Bits nach rechts schieben und dabei kopieren:
      CL_ALLOCA_STACK;
      var uintD* x2_MSDptr;
      var uintD* x2_LSDptr;
      var uintD rounding_bits;
      num_stack_alloc(x2_len, x2_MSDptr=,x2_LSDptr=); // x2_len Digits Platz
      if (j==0)
        { copy_loop_msp(arrayMSDptr(TheLfloat(x2)->data,len),x2_MSDptr,x2_len); rounding_bits = 0; }
        else
        { rounding_bits = shiftrightcopy_loop_msp(arrayMSDptr(TheLfloat(x2)->data,len),x2_MSDptr,x2_len,j,0); }
      // x2_MSDptr/x2_len/x2_LSDptr sind die essentiellen Digits von x2.
      // rounding_bits enthält die letzten j herausgeschobenen Bits.
      // Aus rounding_bits und den nächsten i Digits die 3 Rundungsbits
      // (als Bits intDsize-1..intDsize-3 von rounding_bits) aufbauen:
      if (j>=2)
        // j>=2 -> Bits -1,-2 sind OK, Bit -3 bestimmen:
        { if ((rounding_bits & (bit(intDsize-3)-1)) ==0)
            { if (test_loop_msp(arrayMSDptr(TheLfloat(x2)->data,len) mspop x2_len,i))
                { rounding_bits |= bit(intDsize-3); } // Rundungsbit -3 setzen
            }
            else
            { rounding_bits |= bit(intDsize-3); // Rundungsbit -3 setzen
              rounding_bits &= bitm(intDsize)-bit(intDsize-3); // andere Bits löschen
        }   }
        else
        // j<=3 -> Bits intDsize-4..0 von rounding_bits sind bereits Null.
        // nächstes und weitere i-1 Digits heranziehen:
        { if (i > 0) // i=0 -> Bits -1,-2,-3 sind OK.
            { var uintD* ptr = arrayMSDptr(TheLfloat(x2)->data,len) mspop x2_len;
              rounding_bits |= (mspref(ptr,0) >> j); // weitere relevante Bits des nächsten Digit dazu
              if ((rounding_bits & (bit(intDsize-3)-1)) ==0) // Alle Bits -3,-4,... =0 ?
                { if (   (!((mspref(ptr,0) & (bit(3)-1)) ==0)) // j (<=3) untere Bits von ptr[0] alle =0 ?
                      || test_loop_msp(ptr mspop 1,i-1)
                     )
                    { rounding_bits |= bit(intDsize-3); } // Rundungsbit -3 setzen
                }
                else
                { rounding_bits |= bit(intDsize-3); // Rundungsbit -3 setzen
                  rounding_bits &= bitm(intDsize)-bit(intDsize-3); // andere Bits löschen
        }   }   }
      // x2 liegt in verschobener Form in der UDS x2_MSDptr/x2_len/x2_LSDptr
      // vor, mit Rundungsbits in Bit intDsize-1..intDsize-3 von rounding_bits.
      {var uintD* y_mantMSDptr = arrayMSDptr(TheLfloat(y)->data,len);
       var uintD* y_mantLSDptr = arrayLSDptr(TheLfloat(y)->data,len);
       if (TheLfloat(x1)->sign == TheLfloat(x2)->sign)
         // gleiche Vorzeichen -> Mantissen addieren
         { // erst rechten Mantissenteil (x2_len Digits) durch Addition:
           var uintD carry =
             add_loop_lsp(arrayLSDptr(TheLfloat(x1)->data,len),x2_LSDptr,
                          y_mantLSDptr, x2_len
                         );
           // dann linken Mantissenteil (i Digits) direkt kopieren:
           var uintD* ptr =
             copy_loop_msp(arrayMSDptr(TheLfloat(x1)->data,len),y_mantMSDptr,i);
           // dann Übertrag vom rechten zum linken Mantissenteil addieren:
           if (!(carry==0))
             { if ( inc_loop_lsp(ptr,i) )
                 // Übertrag über das erste Digit hinaus
                 { // Exponent von y incrementieren:
                   if ( ++(TheLfloat(y)->expo) == LF_exp_high+1 ) { throw floating_point_overflow_exception(); }
                   // normalisiere durch Schieben um 1 Bit nach rechts:
                  {var uintD carry_rechts =
                     shift1right_loop_msp(y_mantMSDptr,len,~(uintD)0);
                   rounding_bits = rounding_bits>>1; // Rundungsbits mitschieben
                   if (!(carry_rechts==0)) { rounding_bits |= bit(intDsize-1); }
             }   }}
         }
         else
         // verschiedene Vorzeichen -> Mantissen subtrahieren
         { // erst rechten Mantissenteil (x2_len Digits) durch Subtraktion:
           rounding_bits = -rounding_bits;
           {var uintD carry =
              subx_loop_lsp(arrayLSDptr(TheLfloat(x1)->data,len),x2_LSDptr,
                            y_mantLSDptr, x2_len,
                            (rounding_bits==0 ? 0 : ~(uintD)0)
                           );
            // dann linken Mantissenteil (i Digits) direkt kopieren:
            var uintD* ptr =
              copy_loop_msp(arrayMSDptr(TheLfloat(x1)->data,len),y_mantMSDptr,i);
            // dann Übertrag des rechten vom linken Mantissenteil subtrahieren:
            if (!(carry==0))
              { if ( dec_loop_lsp(ptr,i) )
                  // Übertrag über das erste Digit hinaus, also e1=e2
                  { NOTREACHED } // diesen Fall haben wir schon behandelt
              }
           }
           // UDS y_mantMSDptr/len/y_mantLSDptr/rounding_bits normalisieren:
           {var uintD* ptr = y_mantMSDptr;
            var uintC k = 0;
            var uintC count;
            dotimesC(count,len,
              { if (!(mspref(ptr,0)==0)) goto nonzero_found;
                ptr = ptr mspop 1; k++;
              });
            if (!(rounding_bits==0)) goto nonzero_found;
            // Die UDS ist ganz Null. Also war e1=e2, keine Rundungsbits.
            { NOTREACHED } // diesen Fall haben wir schon behandelt
            nonzero_found: // Digit /=0 gefunden
            // UDS von ptr nach y_mantMSDptr um k Digits nach unten kopieren:
            if (k>0)
              // mindestens ein führendes Nulldigit. Also war e1-e2 = 0 oder 1.
              { ptr = copy_loop_msp(ptr,y_mantMSDptr,len-k); // len-k Digits verschieben
                msprefnext(ptr) = rounding_bits; // Rundungsbits als weiteres Digit
                clear_loop_msp(ptr,k-1); // dann k-1 Nulldigits
                rounding_bits = 0; // und keine weiteren Rundungsbits
                // Exponenten um intDsize*k erniedrigen:
                k = intDsize*k;
               {var uintE uexp = TheLfloat(y)->expo;
                #if !(LF_exp_low==1)
                if (uexp < k+LF_exp_low)
                #else
                if (uexp <= k)
                #endif
                  { if (underflow_allowed())
                      { throw floating_point_underflow_exception(); }
                      else
                      { return encode_LF0(len); } // Ergebnis 0.0
                  }
                TheLfloat(y)->expo = uexp - k;
              }}
           }
           // NUDS y_mantMSDptr/len/y_mantLSDptr/rounding_bits normalisieren:
           {var uintL s;
            integerlengthD(mspref(y_mantMSDptr,0), s = intDsize - );
            // s = Anzahl der führenden Nullbits im ersten Word (>=0, <intDsize)
            if (s > 0)
              { // Muß die NUDS y_mantMSDptr/len/y_mantLSDptr/rounding_bits
                // um s Bits nach links schieben.
                // (Bei e1-e2>1 ist dabei zwangsläufig s=1.)
                if (s==1)
                  { shift1left_loop_lsp(y_mantLSDptr,len);
                    if (rounding_bits & bit(intDsize-1))
                      { lspref(y_mantLSDptr,0) |= bit(0); }
                    rounding_bits = rounding_bits << 1;
                  }
                  else // s>1, also e1-e2 <= 1 <= s.
                  { shiftleft_loop_lsp(y_mantLSDptr,len,s,rounding_bits>>(intDsize-s));
                    rounding_bits = 0; // = rounding_bits << s;
                  }
                // Exponenten um s erniedrigen:
               {var uintE uexp = TheLfloat(y)->expo;
                #if !(LF_exp_low==1)
                if (uexp < s+LF_exp_low)
                #else
                if (uexp <= s)
                #endif
                  { if (underflow_allowed())
                      { throw floating_point_underflow_exception(); }
                      else
                      { return encode_LF0(len); } // Ergebnis 0.0
                  }
                TheLfloat(y)->expo = uexp - s;
              }}
         } }
       // Hier enthält rounding_bits Bit -1 als Bit intDsize-1, Bit -2 als
       // Bit intDsize-2, Bit -3 als Oder(Bits intDsize-3..0) !
       // Runden. Dazu rounding_bits inspizieren:
       if ((rounding_bits & bit(intDsize-1)) ==0) goto ab; // Bit -1 gelöscht -> abrunden
       rounding_bits = rounding_bits<<1; // Bits -2,-3
       if (!(rounding_bits==0)) goto auf; // Bit -2 oder Bit -3 gesetzt -> aufrunden
       // round-to-even:
       if ((lspref(y_mantLSDptr,0) & bit(0)) ==0) goto ab;
       auf: // aufrunden
         if ( inc_loop_lsp(y_mantLSDptr,len) )
           { // Übertrag durchs Aufrunden
             mspref(y_mantMSDptr,0) = bit(intDsize-1); // Mantisse := 10...0
             // Exponent erhöhen:
             if (++(TheLfloat(y)->expo) == LF_exp_high+1) { throw floating_point_overflow_exception(); }
           }
       ab: // abrunden
         ;
      }
      // y fertig.
      return y;
}

}  // namespace cln
