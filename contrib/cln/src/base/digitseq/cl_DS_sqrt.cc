// cl_UDS_sqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/digitseq/cl_DS.h"


// Implementation.

#include "base/cl_low.h"
#include "cln/exception.h"

namespace cln {

// We observe the following timings:
// Time for square root of a_len = 2*N by b_len = N digits,
// OS: Linux 2.2, intDsize==32,        OS: TRU64/4.0, intDsize==64,
// Machine: P-III/450MHz               Machine: EV5/300MHz:
//      N   standard  Newton           standard  Newton
//      30   0.00002   0.00009          0.00011   0.00027
//     100   0.00012   0.00052          0.00057   0.0017
//     300   0.00087   0.0031           0.0037    0.0091
//    1000   0.0089    0.020            0.037     0.069
//    3000   0.087     0.11  <-(~3200)  0.30      0.28  <- (~2750)
//   10000   1.27      0.55             3.5       1.3
//   30000  12.7       1.35            31.1       3.4
// Newton faster for 3200<N            Newton faster for 2750<N
// When in doubt, prefer to choose the standard algorithm.
#if CL_USE_GMP
  static inline bool cl_recipsqrt_suitable (uintC n)
  { return n >= 3200; }
#else
// Use the old default values from CLN version <= 1.0.3 as a crude estimate.
// Time for square root of a_len = 2*N by b_len = N digits,
// on a i486 33 MHz running Linux:
//      N   standard  Newton
//      10    0.00022 0.00132
//      25    0.00082 0.0047
//      50    0.0026  0.0130
//     100    0.0095  0.038
//     250    0.057   0.154
//     500    0.22    0.46
//    1000    0.90    1.39
//    2500    6.0     4.6
//    5000   24.1    10.7
//   10000   98      23.2
//   -----> Newton faster for 1570 <= N <= 1790 and for N >= 2100.
  static inline bool cl_recipsqrt_suitable (uintC n)
  { return n >= 2100; }
#endif

// Bildet zu einer Unsigned Digit sequence a die Wurzel
// (genauer: Gaußklammer aus Wurzel aus a).
// squarep = cl_UDS_sqrt(a_MSDptr,a_len,a_LSDptr, &b);
// > a_MSDptr/a_len/a_LSDptr: eine UDS
// < NUDS b: Gaußklammer der Wurzel aus a
// < squarep: true falls a = b^2, false falls b^2 < a < (b+1)^2.
// Methode:
// erst A normalisieren. A=0 --> B=0, fertig.
// Wähle n so, daß beta^(2n-2) <= A < beta^(2n).
// Wähle s (0<=s<16) so, daß beta^(2n)/4 <= A*2^(2s) < beta^(2n).
// Setze A:=A*2^(2s) und kopiere dabei A. Suche B=floor(sqrt(A)).
// Mache Platz für B=[0,b[n-1],...,b[0]], (mit einem Nulldigit Platz davor,
// da dort nicht B, sondern 2*B abgespeichert werden wird).
// Auf den Plätzen [a[2n-1],...,a[2n-2j]] wird die Differenz
// [a[2n-1],...,a[2n-2j]] - [b[n-1],...,b[n-j]] ^ 2 abgespeichert.
// Bestimme b[n-1] = floor(sqrt(a[2n-1]*beta+a[2n-2])) mit Heron/Newton:
//   {x:=beta als vorheriger Anfangswert, dann:}
//   x := floor((beta+a[2n-1])/2)
//   wiederhole: d:=floor((a[2n-1]*beta+a[2n-2])/x).
//               Falls d<beta (kein Überlauf) und d<x,
//                 setze x:=floor((x+d)/2), nochmals.
//   b[n-1]:=x. In B um ein Bit nach links verschoben abspeichern.
// {Wegen a[2n-1]>=beta/4 ist b[n-1]>=beta/2.}
// Erniedrige [a[2n-1],a[2n-2]] um b[n-1]^2.
// Für j=1,...,n:
//   {Hier [b[n-1],...,b[n-j]] = floor(sqrt(altes [a[2n-1],...,a[2n-2j]])),
//     in [a[2n-1],...,a[2n-2j]] steht jetzt der Rest
//     [a[2n-1],...,a[2n-2j]] - [b[n-1],...,b[n-j]]^2, er ist >=0 und
//     und <= 2 * [b[n-1],...,b[n-j]], belegt daher höchstens j Digits und 1 Bit.
//     Daher sind nur [a[2n-j],...,a[2n-2j]] von Belang.}
//   Für j<n: Bestimme die nächste Ziffer:
//     b* := min(beta-1,floor([a[2n-j],...,a[2n-2j-1]]/(2*[b[n-1],...,b[n-j]]))).
//     und [a[2n-j],...,a[2n-2j-1]] :=
//         [a[2n-j],...,a[2n-2j-1]] - b* * 2 * [b[n-1],...,b[n-j]] (>= 0).
//     Im einzelnen:
//       b* := min(beta-1,floor([a[2n-j],a[2n-j-1],a[2n-j-2]]/(2*b[n-1]))),
//       [a[2n-j],...,a[2n-2j-1]] wie angegeben erniedigen.
//       Solange die Differenz <0 ist, setze b* := b* - 1 und
//         erhöhe [a[2n-j],...,a[2n-2j-1]] um 2 * [b[n-1],...,b[n-j]].
//     Erniedrige [a[2n-j],...,a[2n-2j-2]] um b* ^ 2.
//     Tritt dabei ein negativer Carry auf,
//       so setze b* := b* - 1,
//          setze b[n-j-1] := b* (im Speicher um 1 Bit nach links verschoben),
//          erhöhe [a[2n-j],...,a[2n-2j-2]] um 2*[b[n-1],...,b[n-j-1]]+1.
//       Sonst setze b[n-j-1] := b* (im Speicher um 1 Bit nach links verschoben).
//     Nächstes j.
//   Für j=n:
//     Falls [a[n],...,a[0]] = [0,...,0], ist die Wurzel exakt, sonst nicht.
//     Ergebnis ist [b[n-1],...,b[0]] * 2^(-s), schiebe also im Speicher
//       [b[n],...,b[0]] um s+1 Bits nach rechts.
//     Das Ergebnis ist eine NUDS der Länge n.
bool cl_UDS_sqrt (const uintD* a_MSDptr, uintC a_len, const uintD* a_LSDptr, DS* b_)
{
      // A normalisieren:
      while ((a_len>0) && (mspref(a_MSDptr,0)==0)) { msshrink(a_MSDptr); a_len--; }
      if (a_len==0) // A=0 -> B := NUDS 0
        { b_->LSDptr = b_->MSDptr; b_->len = 0; return true; }
      CL_ALLOCA_STACK;
      // n und s bestimmen:
      var uintC n = ceiling(a_len,2); // a_len = 2n oder 2n-1, n>0.
      var uintL s;
      { var uintD msd = mspref(a_MSDptr,0); // a[2n] bzw. a[2n-1]
        #if 0
        s = 0;
        while /* ((msd & (bit(intDsize-1)|bit(intDsize-2))) ==0) */
              (((sintD)msd >= 0) && ((sintD)(msd<<1) >= 0))
          { msd = msd<<2; s++; }
        #else
        integerlengthD(msd, s = intDsize - ); s = s>>1;
        #endif
      }
      // Noch ist s nur modulo intDsize/2 bestimmt.
      // A um 2s Bits nach links verschoben kopieren:
      var uintD* new_a_MSDptr;
      { var uintD* new_a_LSDptr;
        num_stack_alloc(2*n,new_a_MSDptr=,new_a_LSDptr=); // 2n Digits Platz belegen
       {var uintL shiftcount = 2*s;
        if (!((a_len & bit(0)) ==0)) // a_len ungerade?
          { s += intDsize/2; lsprefnext(new_a_LSDptr) = 0; } // ja -> ein Nulldigit einschieben
        if (shiftcount==0)
          { copy_loop_lsp(a_LSDptr,new_a_LSDptr,a_len); }
          else
          { shiftleftcopy_loop_lsp(a_LSDptr,new_a_LSDptr,a_len,shiftcount); }
      }}
      #define a_MSDptr  new_a_MSDptr
      // Nun ist A = a_MSDptr/2n/..
      if (cl_recipsqrt_suitable(n))
        { // C := 1/sqrt(A) und dann D := A*C näherungsweise errechnen.
          // D evtl. korrigieren, liefert B.
          var uintD* c_MSDptr;
          var uintD* c_LSDptr;
          var uintD* d_MSDptr;
          var uintD* d_LSDptr;
          var uintD* d2_MSDptr;
          num_stack_alloc(n+2, c_MSDptr=,c_LSDptr=);
          num_stack_alloc(2*n+3, d_MSDptr=,d_LSDptr=);
          num_stack_alloc(2*n, d2_MSDptr=,);
          // 1/4 <= a < 1.
          cl_UDS_recipsqrt(a_MSDptr,2*n,c_MSDptr,n);
          // 1 <= c <= 2, | 1/sqrt(a) - c | < 1/2*beta^-n.
          cl_UDS_mul(a_MSDptr mspop (n+1),n+1,c_LSDptr,n+2,d_LSDptr);
          // 1/4 <= d < 2, | sqrt(a) - d | < beta^-n.
          if (mspref(d_MSDptr,0) > 0)
            { dec_loop_lsp(d_MSDptr mspop (n+1),n+1);
              if (mspref(d_MSDptr,0) > 0) throw runtime_exception();
            }
          // D is our guess for B. Square to see how much we have to correct.
          cl_UDS_mul_square(d_MSDptr mspop (1+n),n,d2_MSDptr mspop 2*n);
          // Store D.
          b_->LSDptr = copy_loop_msp(d_MSDptr mspop 1,b_->MSDptr,n);
          b_->len = n;
          // Store 2*D in place of D.
          if (shift1left_loop_lsp(d_MSDptr mspop (1+n),n))
            mspref(d_MSDptr,0) = 1;
          // Compare D^2 against A.
          if (subfrom_loop_lsp(d2_MSDptr mspop 2*n,a_MSDptr mspop 2*n,2*n))
            // guessed too high, decrement D
            { dec_loop_lsp(b_->LSDptr,n);
              dec_loop_lsp(d_MSDptr mspop (1+n),1+n); // store 2*D+1
              if (!addto_loop_lsp(d_MSDptr mspop (1+n),a_MSDptr mspop 2*n,1+n))
                throw runtime_exception();
              if (!inc_loop_lsp(a_MSDptr mspop (n-1),n-1))
                throw runtime_exception();
            }
          else if (test_loop_msp(a_MSDptr,n-1))
            // guessed way too low
            throw runtime_exception();
          else if (compare_loop_msp(a_MSDptr mspop (n-1),d_MSDptr,1+n) > 0)
            // guessed too low, increment D
            { inc_loop_lsp(b_->LSDptr,n);
              mspref(d_MSDptr,n) |= bit(0); // store 2*D-1
              subfrom_loop_lsp(d_MSDptr mspop (1+n),a_MSDptr mspop 2*n,1+n);
              inc_loop_lsp(d_MSDptr mspop (1+n),1+n); // store 2*D
              if (compare_loop_msp(a_MSDptr mspop (n-1),d_MSDptr,1+n) > 0)
                throw runtime_exception();
            }
          else
            // guessed ok
            {}
          // Schiebe b um s Bits nach rechts:
          if (s > 0)
            shiftright_loop_msp(b_->MSDptr,n,s);
          // Teste, ob alle a[n],...,a[0]=0 sind:
          if (test_loop_msp(a_MSDptr mspop (n-1),n+1))
            return false;
          else
            return true; // ja -> Wurzel exakt
        }
      // Platz für B belegen:
      { var uintD* b_MSDptr = b_->MSDptr mspop -1; // ab hier n+1 Digits Platz
        var uintD b_msd;
        // B = [0,b[n-1],...,b[0]] = b_MSDptr/n+1/..
        // Bestimmung von b[n-1]:
        { var uintD a_msd = mspref(a_MSDptr,0); // a[2n-1]
          var uintD a_2msd = mspref(a_MSDptr,1); // a[2n-2]
          #if HAVE_DD
          var uintDD a_msdd = highlowDD(a_msd,a_2msd); // a[2n-1]*beta+a[2n-2]
          #endif
          // Anfangswert: x := floor((beta + a[2n-1])/2)
          var uintD x = floor(a_msd,2) | bit(intDsize-1);
          loop // Heron-Iterationsschleife
            { var uintD d;
              // Dividiere d := floor((a[2n-1]*beta+a[2n-2])/x) :
              if (a_msd>=x) break; // Überlauf -> d>=beta -> fertig
              #if HAVE_DD
                divuD(a_msdd,x, d=,);
              #else
                divuD(a_msd,a_2msd,x, d=,);
              #endif
              if (d >= x) break; // d>=x -> fertig
              // Nächste Iteration: x := floor((x+d)/2)
              // (Da die Folge der x bekanntlich monoton fallend ist
              // und bei b[n-1] >= beta/2 endet, muß x >= beta/2 werden,
              // d.h. x+d>=beta.)
              #if HAVE_DD
                x = (uintD)(floor((uintDD)x + (uintDD)d, 2));
              #else
                x = floor((uintD)(x+d),2) | bit(intDsize-1);
              #endif
            }
          // x = b[n-1] fertig berechnet.
          b_msd = x;
          // Quadrieren und von [a[2n-1],a[2n-2]] abziehen:
          #if HAVE_DD
            a_msdd -= muluD(x,x);
            mspref(a_MSDptr,0) = highD(a_msdd); mspref(a_MSDptr,1) = lowD(a_msdd);
          #else
            {var uintD x2hi;
             var uintD x2lo;
             muluD(x,x, x2hi=,x2lo=);
             mspref(a_MSDptr,0) = a_msd - x2hi;
             if (a_2msd < x2lo) { mspref(a_MSDptr,0) -= 1; }
             mspref(a_MSDptr,1) = a_2msd - x2lo;
            }
          #endif
          mspref(b_MSDptr,0) = 1; mspref(b_MSDptr,1) = x<<1; // b[n-1] ablegen
        }
       {var uintC j = 0;
        var uintD* a_mptr = a_MSDptr mspop 0;
        var uintD* a_lptr = a_MSDptr mspop 2;
        var uintD* b_ptr = b_MSDptr mspop 2;
        // Wurzel-Hauptschleife
        until (++j == n) // j=1,...,n
          { // b_MSDptr = Pointer auf b[n], b_ptr = Pointer hinter b[n-j].
            // a_mptr = Pointer auf a[2n-j], a_lptr = Pointer hinter a[2n-2j].
            // Bestimme b* :
            var uintD b_stern;
            { var uintD a_1d = mspref(a_mptr,0); // a[2n-j], =0 oder =1
              var uintD a_2d = mspref(a_mptr,1); // a[2n-j-1]
              var uintD a_3d = mspref(a_mptr,2); // a[2n-j-2]
              // a[2n-j]*beta^2+a[2n-j-1]*beta+a[2n-j-2] durch 2 dividieren,
              // dann durch b_msd = b[n-1] dividieren:
              #if HAVE_DD
                var uintDD a_123dd = highlowDD(a_2d,a_3d);
                a_123dd = a_123dd>>1; if (!(a_1d==0)) { a_123dd |= bit(2*intDsize-1); }
                if (highD(a_123dd) >= b_msd)
                  { b_stern = bitm(intDsize)-1; } // bei Überlauf: beta-1
                  else
                  { divuD(a_123dd,b_msd, b_stern=,); }
              #else
                a_3d = a_3d>>1; if (!((a_2d & bit(0)) ==0)) { a_3d |= bit(intDsize-1); }
                a_2d = a_2d>>1; if (!(a_1d==0)) { a_2d |= bit(intDsize-1); }
                if (a_2d >= b_msd)
                  { b_stern = bitm(intDsize)-1; } // bei Überlauf: beta-1
                  else
                  { divuD(a_2d,a_3d,b_msd, b_stern=,); }
              #endif
            }
            // b_stern = b* in der ersten Schätzung.
            a_lptr = a_lptr mspop 1; // Pointer hinter a[2n-2j-1]
            // Subtraktion [a[2n-j],...,a[2n-2j-1]] -= b* * [b[n],b[n-1],...,b[n-j]] :
            { var uintD carry = mulusub_loop_lsp(b_stern,b_ptr,a_lptr,j+1);
              if (mspref(a_mptr,0) >= carry)
                { mspref(a_mptr,0) -= carry; }
                else
                { mspref(a_mptr,0) -= carry; // a[2n-j] wird <0
                  // negativer Übertrag -> b* nach unten korrigieren:
                  loop
                    { b_stern = b_stern-1; // b* := b* - 1
                      // erhöhe [a[2n-j],...,a[2n-2j-1]] um [b[n],...,b[n-j]]:
                      if (!(( addto_loop_lsp(b_ptr,a_lptr,j+1) ==0)))
                        if ((mspref(a_mptr,0) += 1) ==0) // Übertrag zu a[2n-j]
                          break; // macht a[2n-j] wieder >=0 -> Subtraktionsergebnis >=0
            }   }   }
            // b_stern = b* in der zweiten Schätzung.
            a_mptr = a_mptr mspop 1; // Pointer auf a[2n-j-1]
            a_lptr = a_lptr mspop 1; // Pointer hinter a[2n-2j-2]
            // Ziehe b* ^ 2 von [a[2n-j],...,a[2n-2j-2]] ab:
            #if HAVE_DD
            { var uintDD b_stern_2 = muluD(b_stern,b_stern);
              var uintDD a_12dd = highlowDD(lspref(a_lptr,1),lspref(a_lptr,0)); // a[2n-2j-1]*beta+a[2n-2j-2]
              var uintDD a_12dd_new = a_12dd - b_stern_2;
              lspref(a_lptr,1) = highD(a_12dd_new); lspref(a_lptr,0) = lowD(a_12dd_new);
              if (a_12dd >= b_stern_2) goto b_stern_ok;
            }
            #else
            { var uintD b_stern_2_hi;
              var uintD b_stern_2_lo;
              muluD(b_stern,b_stern, b_stern_2_hi=,b_stern_2_lo=);
             {var uintD a_1d = lspref(a_lptr,1); // a[2n-2j-1]
              var uintD a_2d = lspref(a_lptr,0); // a[2n-2j-2]
              var uintD a_1d_new = a_1d - b_stern_2_hi;
              var uintD a_2d_new = a_2d - b_stern_2_lo;
              if (a_2d < b_stern_2_lo) { a_1d_new -= 1; }
              lspref(a_lptr,1) = a_1d_new; lspref(a_lptr,0) = a_2d_new;
              if ((a_1d > b_stern_2_hi)
                  || ((a_1d == b_stern_2_hi) && (a_2d >= b_stern_2_lo))
                 )
                goto b_stern_ok;
            }}
            #endif
            if (TRUE)
              { // muß noch [a[2n-j],...,a[2n-2j]] um 1 erniedrigen:
                if ( dec_loop_lsp(a_lptr lspop 2,j+1) ==0) goto b_stern_ok;
                // Subtraktion von b*^2 lieferte negativen Carry
                b_stern = b_stern-1; // b* := b* - 1
                // erhöhe [a[2n-j-1],...,a[2n-2j-2]] um [b[n],...,b[n-j],0] + 2 * b* + 1
                if ((sintD)b_stern < 0) { mspref(b_ptr,-1) |= bit(0); } // höchstes Bit von b* in b[n-j] ablegen
                mspref(b_ptr,0) = (uintD)(b_stern<<1)+1; // niedrige Bits von b* und eine 1 als b[n-j-1] ablegen
                addto_loop_lsp(b_ptr mspop 1,a_lptr,j+2);
                // (a[2n-j] wird nicht mehr gebraucht.)
                mspref(b_ptr,0) -= 1; // niedrige Bits von b* in b[n-j-1] ablegen
                b_ptr = b_ptr mspop 1;
              }
              else
              b_stern_ok:
              { // b* als b[n-j-1] ablegen:
                if ((sintD)b_stern < 0) { mspref(b_ptr,-1) |= bit(0); } // höchstes Bit von b* in b[n-j] ablegen
                mspref(b_ptr,0) = (uintD)(b_stern<<1); // niedrige Bits von b* als b[n-j-1] ablegen
                b_ptr = b_ptr mspop 1;
              }
          }
        // b_MSDptr = Pointer auf b[n], b_ptr = Pointer hinter b[0].
        // a_mptr = Pointer auf a[n].
        // Schiebe [b[n],...,b[0]] um s+1 Bits nach rechts:
        if (s == intDsize-1)
          { lsshrink(b_ptr); }
          else
          { shiftright_loop_msp(b_MSDptr,n+1,s+1); msshrink(b_MSDptr); }
        // b = b_MSDptr/n/b_ptr ist fertig, eine NUDS.
        b_->MSDptr = b_MSDptr; b_->len = n; b_->LSDptr = b_ptr;
        // Teste, ob alle a[n],...,a[0]=0 sind:
        if (test_loop_msp(a_mptr,n+1))
          { return false; }
          else
          { return true; } // ja -> Wurzel exakt
      }}
}
// Bit complexity (N := a_len): O(M(N)).

}  // namespace cln
