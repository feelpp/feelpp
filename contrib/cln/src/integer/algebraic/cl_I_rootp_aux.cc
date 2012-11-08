// cl_rootp_aux().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "base/digit/cl_2D.h"
#include "base/digitseq/cl_2DS.h"

namespace cln {

// Stellt fest, ob ein Integer >=0 eine n-te Potenz ist.
// rootp(x,n,&w)
// > x: ein Integer >=0
// > n: ein Integer >0
// > Annahme: x > 1 und n < (integer-length x).
// < w: Integer (expt x (/ n)) falls x eine n-te Potenz
// < ergebnis: true         ........................, false sonst

bool cl_rootp_aux (cl_I x, uintL n, cl_I* w)
{

// Methode:
// Falls x=0 oder x=1: x = x^n -> JA, x als Ergebnis.
// Hier also x>1. Suche ein Integer y > 1 mit x=y^n.
// Falls n >= integer_length(x): NEIN. (Da y>=2, müßte x>=2^n gelten.)
// Hier also n>0 klein.
// Solange n gerade ist: x := (sqrt x), n := (/ n 2). x nicht ganz -> NEIN.
// Hier also n>0 ungerade.
// Falls n=1: x = x^n -> JA, x als Ergebnis.
// Falls o := ord2(x) nicht durch n teilbar ist: NEIN.
// Sonst dividiere x durch 2^o, am Schluß y mit 2^(o/n) multiplizieren.
// Hier also n>0 ungerade, x ungerade.
// beta := 2^intDsize, m := ceiling(integer_length(x)/(intDsize*n)).
// Suche ein y mit y>=0, y<beta^m mit  x == y^n mod beta^m :
//   Mit dem Hensel-Lemma. Das Polynom f(X) = X^n-x hat die Diskriminante
//   (-1)^((n-1)*n/2) * Res(X^n-x,n*X^(n-1)) = +- n^n * x^(n-1), und diese ist
//   nicht durch p=2 teilbar. Daher ist das Hensel-Lemma mit p=2 anwendbar.
//   Verwende quadratisches Hensel-Lifting, bei linearem Hensel-Lifting der
//   der Verwaltungsaufwand vergleichsweise größer ist und die schnelle
//   Multiplikation nicht zum Zuge kommt.
//   Sei  y0 mod beta^k  mit  y0^n == x mod beta^k  bekannt. k=m -> fertig.
//   Setze  y == y0 + beta^k*y1 mod beta^2k  an, wobei 2k := min(2*k,m).
//   Dabei wird y1 mod beta^(2k-k) so gewählt, daß mod beta^2k
//   x == y^n == y0^n + n * y0^(n-1) * beta^k*y1. Dazu wird
//   (x - y0^n) mod beta^2k errechnet, durch beta^k dividiert (die Division
//   muß nach Voraussetzung an y0 aufgehen) und
//   y1 := ((x-y0^n)/beta^k) / (n*y0^(n-1)) mod beta^(2k-k) gebildet.
//   Damit hat man  (y0 + beta^k*y1)^n == x mod beta^2k . 2k=m -> fertig.
//   Den Anfang (k=1) bekommt man analog, mit beta:=2 und k=1,k=2,k=4,...
// Dann testet man, ob wirklich x = y^n, und ist fertig.

      while ((n % 2) == 0) // n gerade?
        { if (!sqrtp(x,&x)) // Quadratwurzel ziehen versuchen
            { return false; } // nicht ganzzahlig -> fertig
          n = n >> 1; // n := (/ n 2)
        }
      // Nun ist n ungerade.
      if (n==1) { *w = x; return true; } // n=1 -> x als Ergebnis
      var uintC oq = 0; // Shift von y am Schluß
      {var uintC o = ord2(x);
       if (!(o==0))
         {var uintL o_r; divu_3232_3232(o,n, oq=,o_r=); // o_r = o mod n
          if (!(o_r==0)) { return false; } // o nicht durch n teilbar -> fertig
          // oq = o/n.
          // dividiere x durch 2^o:
          x = ash(x,-(sintC)o);
      }  }
      // Nun ist n ungerade, x ungerade.
      CL_ALLOCA_STACK;
      var uintC n_len;
      var uintD* n_LSDptr;
      var uintC x_len;
      var const uintD* x_LSDptr;
      // UDS zu n bilden, 0<n_len<=ceiling(32/intDsize):
      var uintD n_UDS[ceiling(32,intDsize)];
      #if (intDsize==64)
      arrayLSref(n_UDS,1,0) = n; n_LSDptr = arrayLSDptr(n_UDS,1); n_len = 1;
      #else // (intDsize<=32)
      {var uintD* n_MSDptr = arrayMSDptr(n_UDS,32/intDsize);
       set_32_Dptr(n_MSDptr,n); n_LSDptr = arrayLSDptr(n_UDS,32/intDsize);
       n_len = 32/intDsize; // und (zwecks Effizienz) normieren:
       doconsttimes(32/intDsize-1,
         { if (!(msprefnext(n_MSDptr) == 0)) goto n_UDS_ok; n_len--; }
         );
       n_UDS_ok: ; // n_MSDptr/n_len/n_LSDptr ist NUDS zu n.
      }
      #endif
      I_to_NDS_nocopy(x, ,x_len=,x_LSDptr=,false,); // UDS zu x bilden, x_len>0
      var uintD x_lsd = lspref(x_LSDptr,0); // letztes Digit von x
      var uintD y_lsd; // n-te Wurzel von x_lsd mod 2^intDsize
      y_lsd = 1; // Wurzel mod 2^1
      // Für k=1,2,4,...:
      // y_lsd := y_lsd + 2^k * (x_lsd-y_lsd^n)/2^k / (n*y_lsd^(n-1))
      //        = y_lsd + (x_lsd-y_lsd^n) / (n*y_lsd^(n-1))
      doconsttimes(log2_intDsize, // log2(intDsize) Iterationen reichen aus
        { var uintD y_lsd_n1 = expt_pos(y_lsd,n-1); // y_lsd^(n-1)
          var uintD y_lsd_n = mul2adic(y_lsd_n1,y_lsd); // y_lsd^n
          var uintD delta = x_lsd-y_lsd_n; // x_lsd - y_lsd^n
          if (delta==0) goto y_lsd_ok;
          y_lsd = y_lsd + div2adic(delta,mul2adic((uintD)n,y_lsd_n1));
        });
      y_lsd_ok:
      ASSERT(expt_pos(y_lsd,n)==x_lsd);
      // Nun ist y_lsd^n == x_lsd mod beta=2^intDsize.
      { var uintC m = ceiling(x_len,n); // für y nötige Länge, >0, <=x_len
        var uintD* y_LSDptr;
        { var uintD* z1_LSDptr;
          var uintD* z2_LSDptr;
          var uintD* z3_LSDptr;
          num_stack_alloc_1(m, ,y_LSDptr=); // Platz für y
          {var uintC need = 2*m+(32/intDsize-1); // >= max(2*m,m+32/intDsize)
           num_stack_alloc(need, ,z1_LSDptr=); // Platz für Rechenregister 1
           num_stack_alloc(need, ,z2_LSDptr=); // Platz für Rechenregister 2
           num_stack_alloc(need, ,z3_LSDptr=); // Platz für Rechenregister 3
          }
         {var uintC k = 1; // y ist bisher mod beta^k bekannt
          lspref(y_LSDptr,0) = y_lsd; // Startwert von y
          until (k==m)
            { var uintC k2 = 2*k; if (k2>m) { k2=m; } // k2 = min(2*k,m) > k
              // bisheriges y mod beta^k2 mit n-1 potenzieren:
              // Methode für z := y^(n-1) :
              //   zz:=y, e:=n-1.
              //   Solange e gerade, setze zz:=zz*zz, e:=e/2.
              //   z:=zz.
              //   Solange (e:=floor(e/2)) >0,
              //     setze zz:=zz*zz, und falls e ungerade, setze z:=z*zz.
              var uintL e = n-1; // e:=n-1
              var uintD* free_LSDptr = z1_LSDptr;
              var uintD* zz_LSDptr = z2_LSDptr;
              var uintD* z_LSDptr;
              // Ab jetzt {zz_LSDptr,free_LSDptr} = {z1_LSDptr,z2_LSDptr}.
              clear_loop_lsp(y_LSDptr lspop k,k2-k); // y auf k2 Digits erweitern
              copy_loop_lsp(y_LSDptr,zz_LSDptr,k2); // zz:=y
              do { var uintD* new_zz_LSDptr = free_LSDptr;
                   cl_UDS_mul(zz_LSDptr,k2,zz_LSDptr,k2,new_zz_LSDptr); // zz:=zz*zz
                   free_LSDptr = zz_LSDptr; zz_LSDptr = new_zz_LSDptr;
                   e = e>>1; // e:=e/2
                 }
                 while ((e & bit(0)) ==0); // solange e gerade
              z_LSDptr = zz_LSDptr; // z:=zz
              // (zz nicht kopieren; ab der nächsten Veränderung von zz wird
              // {zz_LSDptr,z_LSDptr,free_LSDptr} = {z1_LSDptr,z2_LSDptr,z3_LSDptr}
              // gelten.)
              until ((e = e>>1) == 0)
                { {var uintD* new_zz_LSDptr = free_LSDptr;
                   cl_UDS_mul(zz_LSDptr,k2,zz_LSDptr,k2,new_zz_LSDptr); // zz:=zz*zz
                   free_LSDptr = (z_LSDptr==zz_LSDptr ? z3_LSDptr : zz_LSDptr);
                   zz_LSDptr = new_zz_LSDptr;
                  }
                  if (!((e & bit(0)) == 0))
                    {var uintD* new_z_LSDptr = free_LSDptr;
                     cl_UDS_mul(z_LSDptr,k2,zz_LSDptr,k2,new_z_LSDptr); // z:=z*zz
                     free_LSDptr = z_LSDptr; z_LSDptr = new_z_LSDptr;
                }   }
              // z = y^(n-1) mod beta^k2 ist fertig.
              if (z_LSDptr==zz_LSDptr) { zz_LSDptr = z3_LSDptr; } // zz ist jetzt auch frei
              cl_UDS_mul(z_LSDptr,k2,y_LSDptr,k2,free_LSDptr); // y^n
              sub_loop_lsp(x_LSDptr,free_LSDptr,zz_LSDptr,k2); // zz:=x-y^n
              ASSERT(!DS_test_loop(zz_LSDptr lspop k,k,zz_LSDptr)); // zz == 0 mod beta^k
              cl_UDS_mul(z_LSDptr,k2-k,n_LSDptr,n_len,free_LSDptr); // n*y^(n-1)
              // Quotienten mod beta^(k2-k) bilden und an y mod beta^k ankleben:
              div2adic(k2-k,zz_LSDptr lspop k,free_LSDptr,y_LSDptr lspop k);
              k = k2; // jetzt gilt y^n == x sogar mod beta^k2.
        }}  }
        // y mit y^n == x mod beta^m ist gefunden.
        var cl_I y = UDS_to_I(y_LSDptr lspop m,m); // y als Integer >=0
        // y^n (mit n ungerade) bilden:
        //   c:=a:=y, b:=n.
        //   Solange b:=floor(b/2) >0 ist,
        //     setze a:=a*a, und falls b ungerade, setze c:=a*c.
        //   Liefere c.
        { var cl_I c = y;
          var cl_I a = y;
          until ((n = n>>1) == 0)
            { a = square(a); if (!((n & bit(0)) == 0)) { c = a * c; } }
          // c = y^n
          // mit x vergleichen:
          if (!(x == c))
            // Die ganze Rechnung war umsonst.
            { return false; }
        }
        // y ist tatsächlich n-te Wurzel von x.
        // Noch mit 2^oq multiplizieren:
        if (oq==0) // kein Shift nötig?
          { *w = y; }
          else
          { *w = ash(y,oq); }
        return true;
      }
}

}  // namespace cln
