// gcd().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "base/digit/cl_D.h"
#include "base/cl_xmacros.h"

namespace cln {

#define GCD_ALGO 3  // 1: binär, 2: Schulmethode, 3: Lehmer


#if (GCD_ALGO == 1)

// binäre Methode:
// (gcd a b) :==
// b=0 --> (abs a)
// a=0 --> (abs b)
// sonst:
//   (abs a) und (abs b) in zwei Buffer packen, als Unsigned Digit Sequences.
//   [Schreibe oBdA wieder a,b]
//   (prog ((j 0))
//     1 {a,b >0}
//       (if (evenp a)
//         (if (evenp b)
//           (progn (incf j) (setq a (/ a 2)) (setq b (/ b 2)) (go 1))
//           (go 4)
//         )
//         (while (evenp b) (setq b (/ b 2)))
//       )
//     2 {a,b >0, beide ungerade}
//       (cond ((> a b))
//             ((= a b) (go 5))
//             ((< a b) (rotatef a b))
//       )
//     3 {a,b >0, beide ungerade, a>b}
//       (setq a (- a b))
//     4 {a,b >0, a gerade, b ungerade}
//       (repeat (setq a (/ a 2)) (until (oddp a)))
//       (go 2)
//     5 {a=b>0}
//       (return (ash a j))
//   )
// weil es oft auftritt (insbesondere bei GCD's mehrerer Zahlen):
// a=1 oder b=1 --> 1

  const cl_I gcd (const cl_I& a, const cl_I& b)
    { if (eq(a,1)) { return 1; } // a=1 -> 1
      if (eq(b,1)) { return 1; } // b=1 -> 1
      if (eq(b,0)) { return abs(a); } // b=0 -> (abs a)
      if (eq(a,0)) { return abs(b); } // a=0 -> (abs b)
      CL_ALLOCA_STACK;
      var uintD* a_MSDptr;
      var uintC a_len;
      var uintD* a_LSDptr;
      var uintD* b_MSDptr;
      var uintC b_len;
      var uintD* b_LSDptr;
      // Macro: erzeugt die NUDS zu (abs x), erniedrigt num_stack
      #define I_abs_to_NUDS(x)  \
        I_to_NDS_1(x, x##_MSDptr = , x##_len = , x##_LSDptr = ); /* (nichtleere) NDS holen */\
        if ((sintD)mspref(x##_MSDptr,0) < 0) /* falls <0, negieren:                        */\
          { neg_loop_lsp(x##_LSDptr,x##_len); }                                              \
        if (mspref(x##_MSDptr,0) == 0) /* normalisieren (max. 1 Nulldigit entfernen)       */\
          { msshrink(x##_MSDptr); x##_len--; }
      I_abs_to_NUDS(a); // (abs a) als NUDS erzeugen
      I_abs_to_NUDS(b); // (abs b) als NUDS erzeugen
      // Jetzt ist a = a_MSDptr/a_len/a_LSDptr, b = b_MSDptr/b_len/b_LSDptr,
      // beides NUDS, und a_len>0, b_len>0.
      // Macro: Halbiere x.
      #define halb(x)  \
        { shift1right_loop_msp(x##_MSDptr,x##_len,0); /* um 1 Bit rechts schieben */             \
          if (mspref(x##_MSDptr,0) == 0) { msshrink(x##_MSDptr); x##_len--; } /* normalisieren */\
        }
      // Macro: Ob x gerade ist.
      #define evenp(x)  \
        ((lspref(x##_LSDptr,0) & bit(0)) ==0)
      { var uintL j = 0;
        label_1: // a,b >0
          if (evenp(a))
            { if (evenp(b))
                { j++; halb(a); halb(b); goto label_1; }
                else
                goto label_4;
            }
          while (evenp(b)) { halb(b); }
        label_2: // a,b >0, beide ungerade
          // Vergleiche a und b:
          if (a_len > b_len) goto label_3; // a>b ?
          if (a_len == b_len)
            { var cl_signean vergleich = compare_loop_msp(a_MSDptr,b_MSDptr,a_len);
              if (vergleich > 0) goto label_3; // a>b ?
              if (vergleich == 0) goto label_5; // a=b ?
            }
          // a<b -> a,b vertauschen:
          swap(uintD*, a_MSDptr,b_MSDptr);
          swap(uintC, a_len,b_len);
          swap(uintD*, a_LSDptr,b_LSDptr);
        label_3: // a,b >0, beide ungerade, a>b
          // subtrahiere a := a - b
          if (!( subfrom_loop_lsp(b_LSDptr,a_LSDptr,b_len) ==0))
            { dec_loop_lsp(a_LSDptr lspop b_len,a_len-b_len); }
          while (mspref(a_MSDptr,0) == 0) { msshrink(a_MSDptr); a_len--; } // normalisieren
        label_4: // a,b >0, a gerade, b ungerade
          do { halb(a); } while (evenp(a));
          goto label_2;
        label_5: // a=b>0
          // a zu einer NDS machen:
          return ash(NUDS_to_I(a_MSDptr,a_len),j); // ggT der ungeraden Anteile als Integer, mal 2^j
      }
      #undef evenp
      #undef halb
      #undef I_abs_to_NUDS
    }

#endif /* GCD_ALGO == 1 */


#if (GCD_ALGO == 2)

// Schulmethode:
//   (gcd a b) :==
//   [a:=(abs a), b:=(abs b), while b>0 do (a,b) := (b,(mod a b)), -> a]
// verbessert:
// a=0 -> (abs b)
// b=0 -> (abs a)
// a=1 -> 1
// b=1 -> 1
// a:=(abs a), b:=(abs b)
// Falls a=b: return a; falls a<b: vertausche a und b.
// (*) {Hier a>b>0}
// Falls b=1, return 1. {spart eine Division durch 1}
// Sonst dividieren (divide a b), a:=b, b:=Rest.
//       Falls b=0, return a, sonst goto (*).

  const cl_I gcd (const cl_I& a, const cl_I& b)
    { if (eq(a,1)) { return 1; } // a=1 -> 1
      if (eq(b,1)) { return 1; } // b=1 -> 1
      if (eq(b,0)) { return abs(a); } // b=0 -> (abs a)
      if (eq(a,0)) { return abs(b); } // a=0 -> (abs b)
      // Beträge nehmen:
     {var cl_I abs_a = abs(a);
      var cl_I abs_b = abs(b);
      var cl_I& a = abs_a;
      var cl_I& b = abs_b;
      if (fixnump(a) && fixnump(b)) // ggT zweier Fixnums >0
        { // bleibt Fixnum, da (gcd a b) <= (min a b)
          return V_to_FN(gcd(FN_to_UV(a),FN_to_UV(b)));
        }
      { var cl_signean vergleich = compare(a,b);
        if (vergleich == 0) { return a; } // a=b -> fertig
        if (vergleich < 0) { var cl_I tmp = a; a = b; b = a; } // a<b -> a,b vertauschen
      }
      loop // Hier a > b > 0
        { if (eq(b,1)) { return 1; } // b=1 -> Ergebnis 1
          { var cl_I_div_t div = cl_divide(a,b); a = b; b = div.remainder; }
          if (eq(b,0)) { return a; }
        }
    }}

#endif /* GCD_ALGO == 2 */


#if (GCD_ALGO == 3)

// Lehmer-Methode:
// vgl. [ D. E. Knuth: The Art of Computer Programming, Vol. 2: Seminumerical
//        Algorithms, Sect. 4.5.2., Algorithm L ]
// und [ Collins, Loos: SAC-2, Algorithms IGCD, DPCC ].
// (gcd a b) :==
// a=0 -> (abs b)
// b=0 -> (abs a)
// a=1 -> 1
// b=1 -> 1
// a:=(abs a), b:=(abs b)
// (*) {Hier a,b>0}
// Falls a=b: return a; falls a<b: vertausche a und b.
// {Hier a>b>0}
// Falls (- (integer-length a) (integer-length b)) >= intDsize/2,
//   lohnt sich eine Division: (a,b) := (b , a mod b). Falls b=0: return a.
// Falls dagegen 0 <= (- (integer-length a) (integer-length b)) < intDsize/2,
//   seien a' die führenden intDsize Bits von a
//   (2^(intDsize-1) <= a' < 2^intDsize) und b' die entsprechenden Bits von b
//   (2^(intDsize/2) <= b' <= a' < 2^intDsize).
//   Rechne den Euklid-Algorithmus mit Beifaktoren für ALLE Zahlen (a,b) aus,
//   die mit a' bzw. b' anfangen; das liefert x1,y1,x2,y2, so daß
//   ggT(a,b) = ggT(x1*a-y1*b,-x2*a+y2*b) und x1*a-y1*b>=0,-x2*a+y2*b>=0.
//   Genauer: Mit offensichtlicher Skalierung betrachten wir
//            a als beliebiges Element des Intervalls [a',a'+1) und
//            b als beliebiges Element des Intervalls [b',b'+1) und
//            führen den Euklid-Algorithmus schrittweise durch:
//            (x1,y1,z1) := (1,0,a'), (x2,y2,z2) := (0,1,b'),
//            Schleife:
//            {Hier x1*a'-y1*b'=z1, x1*a-y1*b in [z1-y1,z1+x1), z1-y1>=0, z1>0,
//             und -x2*a'+y2*b'=z2, -x2*a+y2*b in [z2-x2,z2+y2), z2-x2>=0, z2>0,
//             x1*y2-x2*y1=1, x1*z2+x2*z1=b', y1*z2+y2*z1=a'.}
//            Falls z1-y1>=z2+y2:
//              (x1,y1,z1) := (x1+x2,y1+y2,z1-z2), goto Schleife.
//            Falls z2-x2>=z1+x1:
//              (x2,y2,z2) := (x2+x1,y2+y1,z2-z1), goto Schleife.
//            Sonst muß man abbrechen.
//            {Zu den Schleifeninvarianten:
//             1. Die Gleichungen x1*a'-y1*b'=z1, -x2*a'+y2*b'=z2,
//                x1*y2-x2*y1=1, x1*z2+x2*z1=b', y1*z2+y2*z1=a' mit Induktion.
//             2. Die Ungleichungen x1>0, y1>=0, x2>=0, y2>0 mit Induktion.
//             3. Die Ungleichungen z1>=0, z2>=0 nach Fallauswahl.
//             4. Die Ungleichung x1+x2>0 aus x1*z2+x2*z1=b'>0,
//                die Ungleichung y1+y2>0 aus y1*z2+y2*z1=a'>0.
//             5. Die Ungleichung z1>0 wegen Fallauswahl und y1+y2>0,
//                Die Ungleichung z2>0 wegen Fallauswahl und x1+x2>0.
//             6. Die Ungleichungen z1-y1>=0, z2-x2>=0 wegen Fallauswahl.
//             7. Die Ungleichung max(z1,z2) <= a' mit Induktion.
//             8. Die Ungleichung x1+x2 <= x1*z2+x2*z1 = b',
//                die Ungleichung y1+y2 <= y1*z2+y2*z1 = a'.
//             Damit bleiben alle Größen im Intervall [0,beta), kein Überlauf.
//             9. Die Ungleichungen z1+x1<=beta, z2+y2<=beta mit Induktion.
//             10. x1*a-y1*b in (z1-y1,z1+x1) (bzw. [z1,z1+x1) bei y1=0),
//                -x2*a+y2*b in (z2-x2,z2+y2) (bzw. [z2,z2+y2) bei x2=0),
//                da a in a'+[0,1) und b in b'+[0,1).
//                Jedenfalls 0 < x1*a-y1*b < z1+x1 <= x2*z1+x1*z2 = b' falls x2>0,
//                und        0 < -x2*a+y2*b < z2+y2 <= y1*z2+y2*z1 = a' falls y1>0.}
//            Man kann natürlich auch mehrere Subtraktionsschritte auf einmal
//            durchführen:
//            Falls q := floor((z1-y1)/(z2+y2)) > 0 :
//              (x1,y1,z1) := (x1+q*x2,y1+q*y2,z1-q*z2), goto Schleife.
//            Falls q := floor((z2-x2)/(z1+x1)) > 0 :
//              (x2,y2,z2) := (x2+q*x1,y2+q*y1,z2-q*z1), goto Schleife.
//            {Am Schluß gilt -(x1+x2) < z1-z2 < y1+y2 und daher
//             z2-x2 <= b'/(x1+x2) < z1+x1, z1-y1 <= a'/(y1+y2) < z2+y2,
//             und - unter Berücksichtigung von x1*y2-x2*y1=1 -
//             z1-y1 <= b'/(x1+x2) < z2+y2, z2-x2 <= a'/(y1+y2) < z1+x1,
//             also  max(z1-y1,z2-x2) <= min(b'/(x1+x2),a'/(y1+y2))
//                          <= max(b'/(x1+x2),a'/(y1+y2)) < min(z1+x1,z2+y2).}
//   Im Fall y1=x2=0 => x1=y2=1 (der nur bei a'=z1=z2=b' eintreten kann)
//     ersetze (a,b) := (a-b,b). {Beide >0, da a>b>0 war.}
//   Der Fall y1=0,x2>0 => x1=y2=1 => a' = z1 < z2+x2*z1 = b'
//     kann nicht eintreten.
//   Im Fall x2=0,y1>0 => x1=y2=1 ersetze (a,b) := (a-y1*b,b).
//     {Das ist OK, da 0 <= z1-y1 = a'-y1*(b'+1) < a-y1*b < a.}
//   Sonst (y1>0,x2>0) ersetze (a,b) := (x1*a-y1*b,-x2*a+y2*b).
//     {Das ist OK, da 0 <= z1-y1 = x1*a'-y1*(b'+1) < x1*a-y1*b
//                  und 0 <= z2-x2 = -x2*(a'+1)+y2*b' < -x2*a+y2*b
//      und x1*a-y1*b < x1*(a'+1)-y1*b' = z1+x1 <= x2*z1+x1*z2 = b' <= b
//      und -x2*a+y2*b < -x2*a'+y2*(b'+1) = z2+y2 <= y1*z2+y2*z1 = a' <= a.}
// goto (*).

  // Define this to 1 in order to use double-word sized a' and b'.
  // This gives better x1,y1,x2,y2, because normally the values x1,y1,x2,y2
  // have only about intDsize/2 bits and so half of the multiplication work
  // is lost. Actually, this flag multiplies the gcd speed by 1.5, not 2.0.
  #define DOUBLE_SPEED 1
  // Speed increases only for large integers. For small ones, computing
  // with double-word sized a' and b' is too costly. The threshold is
  // between 12 and 20, around 15.
  #define cl_gcd_double_threshold  16

  // gcd of two single-word numbers >0
  static uintD gcdD (uintD a, uintD b)
    { var uintD bit_j = (a | b); // endet mit einer 1 und j Nullen
      bit_j = bit_j ^ (bit_j - 1); // Maske = bit(j) | bit(j-1) | ... | bit(0)
      if (!((a & bit_j) ==0))
        { if (!((b & bit_j) ==0)) goto odd_odd; else goto odd_even; }
      if (!((b & bit_j) ==0)) goto even_odd;
      NOTREACHED;
      loop
        { odd_odd: // a,b >0, beide ungerade
          // Vergleiche a und b:
          if (a == b) break; // a=b>0 -> fertig
          if (a > b) // a>b ?
            { a = a-b;
              even_odd: // a,b >0, a gerade, b ungerade
              do { a = a>>1; } while ((a & bit_j) ==0);
            }
            else // a<b
            { b = b-a;
              odd_even: // a,b >0, a ungerade, b gerade
              do { b = b>>1; } while ((b & bit_j) ==0);
            }
        }
      // a=b>0
      return a;
    }

  const cl_I gcd (const cl_I& a, const cl_I& b)
    { if (eq(a,1)) { return 1; } // a=1 -> 1
      if (eq(b,1)) { return 1; } // b=1 -> 1
      if (eq(b,0)) { return abs(a); } // b=0 -> (abs a)
      if (eq(a,0)) { return abs(b); } // a=0 -> (abs b)
      if (fixnump(a) && fixnump(b)) // ggT zweier Fixnums /=0
        { var sintV a_ = FN_to_V(a);
          if (a_ < 0) { a_ = -a_; }
          var sintV b_ = FN_to_V(b);
          if (b_ < 0) { b_ = -b_; }
          return UV_to_I(gcd((uintV)a_,(uintV)b_));
        }
      CL_ALLOCA_STACK;
      var uintD* a_MSDptr;
      var uintC a_len;
      var uintD* a_LSDptr;
      var uintD* b_MSDptr;
      var uintC b_len;
      var uintD* b_LSDptr;
      // Macro: erzeugt die NUDS zu (abs x), erniedrigt num_stack
      #define I_abs_to_NUDS(x)  \
        I_to_NDS_1(x, x##_MSDptr = , x##_len = , x##_LSDptr = ); /* (nichtleere) NDS holen */\
        if ((sintD)mspref(x##_MSDptr,0) < 0) /* falls <0, negieren:                        */\
          { neg_loop_lsp(x##_LSDptr,x##_len); }                                              \
        if (mspref(x##_MSDptr,0) == 0) /* normalisieren (max. 1 Nulldigit entfernen)       */\
          { msshrink(x##_MSDptr); x##_len--; }
      I_abs_to_NUDS(a); // (abs a) als NUDS erzeugen
      I_abs_to_NUDS(b); // (abs b) als NUDS erzeugen
      // Jetzt ist a = a_MSDptr/a_len/a_LSDptr, b = b_MSDptr/b_len/b_LSDptr,
      // beides NUDS, und a_len>0, b_len>0.
      // Platz für zwei Rechenregister besorgen, mit je max(a_len,b_len)+1 Digits:
      {var uintD* divroomptr; // Platz für Divisionsergebnis
       var uintD* c_LSDptr;
       var uintD* d_LSDptr;
       {var uintC c_len = (a_len>=b_len ? a_len : b_len) + 1;
        num_stack_alloc(c_len,divroomptr=,c_LSDptr=);
        num_stack_alloc(c_len,,d_LSDptr=);
        // Jetzt ist ../c_len/c_LSDptr, ../c_len/d_LSDptr frei.
       }
       loop
         { // Hier a,b>0, beides NUDS.
           // Vergleiche a und b:
           if (a_len > b_len) goto a_greater_b; // a>b ?
           if (a_len == b_len)
             { var cl_signean vergleich = compare_loop_msp(a_MSDptr,b_MSDptr,a_len);
               if (vergleich > 0) goto a_greater_b; // a>b ?
               if (vergleich == 0) break; // a=b ?
             }
           // a<b -> a,b vertauschen:
           swap(uintD*, a_MSDptr,b_MSDptr);
           swap(uintC, a_len,b_len);
           swap(uintD*, a_LSDptr,b_LSDptr);
           a_greater_b:
           // Hier a>b>0, beides NUDS.
           if (b_len==1) // Beschleunigung eines häufigen Falles
             { var uintD b0 = mspref(b_MSDptr,0);
               if (b0==1)
                  // a>b=1 -> Ergebnis 1.
                 { return 1; }
               // a>b>1 -> evtl. Division durch b
               var uintD a0;
               if (a_len==1)
                 { a0 = mspref(a_MSDptr,0); }
                 else
                 { a0 = divu_loop_msp(b0,a_MSDptr,a_len);
                   if (a0==0)
                     { return UD_to_I(b0); }
                 }
               return UD_to_I(gcdD(a0,b0));
             }
           // Entscheidung, ob Division oder Linearkombination:
           { var uintD a_msd; // führende intDsize Bits von a
             var uintD b_msd; // entsprechende Bits von b
             #if DOUBLE_SPEED
             var uintD a_nsd; // nächste intDsize Bits von a
             var uintD b_nsd; // entsprechende Bits von b
             #endif
             { var uintC len_diff = a_len-b_len; // Längendifferenz
               if (len_diff > 1) goto divide; // >=2 -> Bitlängendifferenz>intDsize -> dividieren
               #define bitlendiff_limit  (intDsize/2) // sollte >0,<intDsize sein
              {var uintC a_msd_size;
               a_msd = mspref(a_MSDptr,0); // führendes Digit von a
               integerlengthD(a_msd,a_msd_size=); // dessen Bit-Länge (>0,<=intDsize) berechnen
               b_msd = mspref(b_MSDptr,0);
               #if HAVE_DD
               {var uintDD b_msdd = // 2 führende Digits von b
                  (len_diff==0
                   ? highlowDD(b_msd, mspref(b_MSDptr,1))
                   : (uintDD)b_msd
                  );
                // a_msd_size+intDsize - b_msdd_size >= bitlendiff_limit -> dividieren:
                b_msd = lowD(b_msdd >> a_msd_size);
                if (b_msd < (uintD)bit(intDsize-bitlendiff_limit)) goto divide;
                #if DOUBLE_SPEED
                b_nsd = lowD(highlowDD(lowD(b_msdd), (b_len<=2-len_diff ? 0 : mspref(b_MSDptr,2-len_diff))) >> a_msd_size);
                #endif
               }
               {var uintDD a_msdd = // 2 führende Digits von a
                  highlowDD(a_msd, mspref(a_MSDptr,1));
                a_msd = lowD(a_msdd >> a_msd_size);
                #if DOUBLE_SPEED
                a_nsd = lowD(highlowDD(lowD(a_msdd), (a_len<=2 ? 0 : mspref(a_MSDptr,2))) >> a_msd_size);
                #endif
               }
               if (a_msd == b_msd) goto subtract;
               #else
               if (len_diff==0)
                 { // a_msd_size - b_msd_size >= bitlendiff_limit -> dividieren:
                   if ((a_msd_size > bitlendiff_limit)
                       && (b_msd < (uintD)bit(a_msd_size-bitlendiff_limit))
                      )
                     goto divide;
                   // Entscheidung für Linearkombination ist gefallen.
                   // a_msd und b_msd so erweitern, daß a_msd die führenden
                   // intDsize Bits von a enthält:
                  {var uintC shiftcount = intDsize-a_msd_size; // Shiftcount nach links (>=0, <intDsize)
                   if (shiftcount>0)
                     { a_msd = a_msd << shiftcount;
                       b_msd = b_msd << shiftcount;
                       a_msd |= mspref(a_MSDptr,1) >> a_msd_size;
                       b_msd |= mspref(b_MSDptr,1) >> a_msd_size;
                     }
                   if (a_msd == b_msd) goto subtract;
                   #if DOUBLE_SPEED
                   a_nsd = mspref(a_MSDptr,1);
                   b_nsd = mspref(b_MSDptr,1);
                   if (shiftcount>0)
                     { a_nsd = a_nsd << shiftcount;
                       b_nsd = b_nsd << shiftcount;
                       if (a_len>2)
                         { a_nsd |= mspref(a_MSDptr,2) >> a_msd_size;
                           b_nsd |= mspref(b_MSDptr,2) >> a_msd_size;
                     }   }
                   #endif
                 }}
                 else
                 // len_diff=1
                 { // a_msd_size+intDsize - b_msd_size >= bitlendiff_limit -> dividieren:
                   if ((a_msd_size >= bitlendiff_limit)
                       || (b_msd < (uintD)bit(a_msd_size+intDsize-bitlendiff_limit))
                      )
                     goto divide;
                   // Entscheidung für Linearkombination ist gefallen.
                   // a_msd und b_msd so erweitern, daß a_msd die führenden
                   // intDsize Bits von a enthält:
                   // 0 < a_msd_size < b_msd_size + bitlendiff_limit - intDsize <= bitlendiff_limit < intDsize.
                   a_msd = (a_msd << (intDsize-a_msd_size)) | (mspref(a_MSDptr,1) >> a_msd_size);
                   #if DOUBLE_SPEED
                   a_nsd = mspref(a_MSDptr,1) << (intDsize-a_msd_size);
                   b_nsd = b_msd << (intDsize-a_msd_size);
                   a_nsd |= mspref(a_MSDptr,2) >> a_msd_size;
                   b_nsd |= mspref(b_MSDptr,1) >> a_msd_size;
                   #endif
                   b_msd = b_msd >> a_msd_size;
                 }
               #endif
               #undef bitlendiff_limit
             }}
             // Nun ist a_msd = a' > b' = b_msd.
             { // Euklid-Algorithmus auf den führenden Digits durchführen:
               var partial_gcd_result likobi;
               #if DOUBLE_SPEED
               if (a_len >= cl_gcd_double_threshold)
                 {
                   #if HAVE_DD
                   partial_gcd(highlowDD(a_msd,a_nsd),highlowDD(b_msd,b_nsd),&likobi); // liefert x1,y1,x2,y2
                   #else
                   partial_gcd(a_msd,a_nsd,b_msd,b_nsd,&likobi); // liefert x1,y1,x2,y2
                   #endif
                 }
               else
               #endif
                 { partial_gcd(a_msd,b_msd,&likobi); } // liefert x1,y1,x2,y2, aber nur halb so gut
               // Hier y1>0.
               if (likobi.x2==0)
                 { // Ersetze (a,b) := (a-y1*b,b).
                   if (likobi.y1==1) goto subtract; // einfacherer Fall
                   // Dazu evtl. a um 1 Digit erweitern, so daß a_len=b_len+1:
                   if (a_len == b_len) { lsprefnext(a_MSDptr) = 0; a_len++; }
                   // und y1*b von a subtrahieren:
                   mspref(a_MSDptr,0) -= mulusub_loop_lsp(likobi.y1,b_LSDptr,a_LSDptr,b_len);
                 }
                 else
                 { // Ersetze (a,b) := (x1*a-y1*b,-x2*a+y2*b).
                   // Dazu evtl. b um 1 Digit erweitern, so daß a_len=b_len:
                   if (!(a_len==b_len)) { lsprefnext(b_MSDptr) = 0; b_len++; }
                   // c := x1*a-y1*b bilden:
                   mulu_loop_lsp(likobi.x1,a_LSDptr,c_LSDptr,a_len);
                   /* lspref(c_LSDptr,a_len) -= */
                     mulusub_loop_lsp(likobi.y1,b_LSDptr,c_LSDptr,a_len);
                   // d := -x2*a+y2*b bilden:
                   mulu_loop_lsp(likobi.y2,b_LSDptr,d_LSDptr,a_len);
                   /* lspref(d_LSDptr,a_len) -= */
                     mulusub_loop_lsp(likobi.x2,a_LSDptr,d_LSDptr,a_len);
                   // Wir wissen, daß 0 < c < b und 0 < d < a. Daher müßten
                   // lspref(c_LSDptr,a_len) und lspref(d_LSDptr,a_len) =0 sein.
                   // a := c und b := d kopieren:
                   copy_loop_lsp(c_LSDptr,a_LSDptr,a_len);
                   copy_loop_lsp(d_LSDptr,b_LSDptr,a_len);
                   // b normalisieren:
                   while (mspref(b_MSDptr,0)==0) { msshrink(b_MSDptr); b_len--; }
             }   }
             if (false)
               { subtract: // Ersetze (a,b) := (a-b,b).
                 if (!( subfrom_loop_lsp(b_LSDptr,a_LSDptr,b_len) ==0))
                   // Übertrag nach b_len Stellen, muß also a_len=b_len+1 sein.
                   { mspref(a_MSDptr,0) -= 1; }
               }
             // a normalisieren:
             while (mspref(a_MSDptr,0)==0) { msshrink(a_MSDptr); a_len--; }
           }
           if (false)
             { divide: // Ersetze (a,b) := (b , a mod b).
              {var uintD* old_a_LSDptr = a_LSDptr;
               var DS q;
               var DS r;
               cl_UDS_divide(a_MSDptr,a_len,a_LSDptr,b_MSDptr,b_len,b_LSDptr, divroomptr, &q,&r);
               a_MSDptr = b_MSDptr; a_len = b_len; a_LSDptr = b_LSDptr; // a := b
               b_len = r.len; if (b_len==0) break; // b=0 -> fertig
               b_LSDptr = old_a_LSDptr; // b übernimmt den vorherigen Platz von a
               b_MSDptr = copy_loop_lsp(r.LSDptr,b_LSDptr,b_len); // b := r kopieren
               goto a_greater_b; // Nun ist a>b>0
             }}
      }  }
      return NUDS_to_I(a_MSDptr,a_len); // NUDS a als Ergebnis
      #undef I_abs_to_NUDS
    }

#endif /* GCD_ALGO == 3 */

}  // namespace cln
