// xgcd().

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


#if (GCD_ALGO == 2)

// Schulmethode:
//   (gcd A B) :==
//   [a:=(abs A), b:=(abs B), while b>0 do (a,b) := (b,(mod a b)), -> a]
// verbessert:
// A=1 -> return g=1, (u,v)=(1,0)
// B=1 -> return g=1, (u,v)=(0,1)
// a:=(abs A), ua:=(signum A), va:=0
// b:=(abs B), ub:=0, vb:=(signum B)
// A=0 -> return g=b, (u,v) = (ub,vb)
// B=0 -> return g=a, (u,v) = (ua,va)
// {Stets ua*A+va*B=a, ub*A+vb*B=b, ua*vb-ub*va = +/- 1.}
// Falls a=b: return a,ua,va;
// falls a<b: vertausche a und b, ua und ub, va und vb.
// (*) {Hier a>b>0}
// Falls b=1, return 1,ub,vb. {spart eine Division durch 1}
// Sonst dividieren (divide a b) -> q,r.
//       Falls r=0, return b,ub,vb.
//       a:=b, b := Rest r = a-q*b, (ua,va,ub,vb) := (ub,vb,ua-q*ub,va-q*vb).
//       goto (*).

  const cl_I xgcd (const cl_I& a, const cl_I& b, cl_I* u, cl_I* v)
    { if (eq(a,1)) // a=1 -> g=1, (u,v)=(1,0)
        { *u = 1; *v = 0; return a; }
      if (eq(b,1)) // b=1 -> g=1, (u,v)=(0,1)
        { *u = 0; *v = 1; return b; }
      // Vorzeichen nehmen:
      var cl_I ua = (minusp(a) ? cl_I(-1) : cl_I(1)); // ua := +/- 1
      var cl_I va = 0;
      var cl_I ub = 0;
      var cl_I vb = (minusp(b) ? cl_I(-1) : cl_I(1)); // vb := +/- 1
      // Beträge nehmen:
     {var cl_I abs_a = abs(a);
      var cl_I abs_b = abs(b);
      var cl_I& a = abs_a;
      var cl_I& b = abs_b;
      if (eq(b,0)) // b=0 -> g=a, (u,v) = (ua,va)
        { *u = ua; *v = va; return a; }
      if (eq(a,0)) // a=0 -> g=b, (u,v) = (ub,vb)
        { *u = ub; *v = vb; return b; }
      { var cl_signean vergleich = compare(a,b);
        if (vergleich == 0) // a=b -> fertig
          { *u = ua; *v = va; return a; }
        if (vergleich < 0) // a<b -> a,b vertauschen
          { swap(cl_I,a,b); swap(cl_I,ua,ub); swap(cl_I,va,vb); }
      }
      loop // Hier a>b>0
        { if (eq(b,1)) // b=1 -> g=b, (u,v) = (ub,vb)
            { *u = ub; *v = vb; return b; }
          var cl_I_div_t div = cl_divide(a,b); // Division a / b
          var cl_I& q = div.quotient;
          var cl_I& r = div.remainder;
          if (eq(r,0)) // r=0 -> fertig
            { *u = ub; *v = vb; return b; }
          { var cl_I x = ua-q*ub; ua = ub; ub = x; }
          { var cl_I x = va-q*vb; va = vb; vb = x; }
          a = b; b = r;
        }
    }}

#endif /* GCD_ALGO == 2 */


#if (GCD_ALGO == 3)
// (xgcd A B) :==
// wie oben bei (gcd A B).
// Zusätzlich werden Variablen sA,sB,sk,uAa,uBa,uAb,uBb geführt,
// wobei sA,sB,sk Vorzeichen (+/- 1) und uAa,uBa,uAb,uBb Integers >=0 sind mit
//     uAa * sA*A - uBa * sB*B = a,
//   - uAb * sA*A + uBb * sB*B = b,
// ferner  uAa * uBb - uAb * uBa = sk  und daher (Cramersche Regel)
//   uBb * a + uBa * b = sk*sA*A, uAb * a + uAa * b = sk*sB*B.
// Zu Beginn (a,b) := (|A|,|B|), (sA,sB) := ((signum A), (signumB)),
//           (uAa,uBa,uAb,uBb) := (1,0,0,1).
// Beim Ersetzen (a,b) := (a-b,b)
//   ersetzt man (uAa,uBa,uAb,uBb) := (uAa+uAb,uBa+uBb,uAb,uBb).
// Beim Ersetzen (a,b) := (a-y1*b,b)
//   ersetzt man (uAa,uBa,uAb,uBb) := (uAa+y1*uAb,uBa+y1*uBb,uAb,uBb).
// Beim Ersetzen (a,b) := (x1*a-y1*b,-x2*a+y2*b) mit x1*y2-x2*y1=1
//   ersetzt man (uAa,uBa,uAb,uBb) :=
//               (x1*uAa+y1*uAb,x1*uBa+y1*uBb,x2*uAa+y2*uAb,x2*uBa+y2*uBb).
// Beim Ersetzen (a,b) := (b,a)
//   ersetzt man (uAa,uBa,uAb,uBb) := (uAb,uBb,uAa,uBa),
//               sk := -sk, (sA,sB) := (-sA,-sB).
// Beim Ersetzen (a,b) := (b,a-q*b)
//   ersetzt man (uAa,uBa,uAb,uBb) := (uAb,uBb,uAa+q*uAb,uBa+q*uBb),
//               sk := -sk, (sA,sB) := (-sA,-sB).
// Zum Schluß ist a der ggT und a = uAa*sA * A + -uBa*sB * B
// die gewünschte Linearkombination.
// Da stets gilt sk*sA*A = |A|, sk*sB*B = |B|, a>=1, b>=1,
// folgt 0 <= uAa <= |B|, 0 <= uAb <= |B|, 0 <= uBa <= |A|, 0 <= uBb <= |A|.
// Ferner wird sk nie benutzt, braucht also nicht mitgeführt zu werden.

  // Define this to 1 in order to use double-word sized a' and b'.
  // This gives better x1,y1,x2,y2, because normally the values x1,y1,x2,y2
  // have only about intDsize/2 bits and so half of the multiplication work
  // is lost. Actually, this flag multiplies the gcd speed by 1.5, not 2.0.
  #define DOUBLE_SPEED 1

  // Bildet u := u + v, wobei für u genügend Platz sei:
  // (Benutzt v.MSDptr nicht.)
  static void NUDS_likobi0_NUDS (DS* u, DS* v)
    { var uintC u_len = u->len;
      var uintC v_len = v->len;
      if (u_len >= v_len)
        { if (!( addto_loop_lsp(v->LSDptr,u->LSDptr,v_len) ==0))
            { if (!( inc_loop_lsp(u->LSDptr lspop v_len,u_len-v_len) ==0))
                { lsprefnext(u->MSDptr) = 1; u->len++; }
        }   }
        else // u_len <= v_len
        { u->MSDptr = copy_loop_lsp(v->LSDptr lspop u_len,u->LSDptr lspop u_len,v_len-u_len);
          u->len = v_len;
          if (!( addto_loop_lsp(v->LSDptr,u->LSDptr,u_len) ==0))
            { if (!( inc_loop_lsp(u->LSDptr lspop u_len,v_len-u_len) ==0))
                { lsprefnext(u->MSDptr) = 1; u->len++; }
        }   }
    }

  // Bildet u := u + q*v, wobei für u genügend Platz sei:
  // (Dabei sei nachher u>0.)
  static void NUDS_likobi1_NUDS (DS* u, DS* v, uintD q)
    { var uintC v_len = v->len;
      if (v_len>0) // nur nötig, falls v /=0
        { var uintC u_len = u->len;
          var uintD carry;
          if (u_len <= v_len) // evtl. u vergrößern
            { u->MSDptr = clear_loop_lsp(u->MSDptr,v_len-u_len+1);
              u->len = u_len = v_len+1;
            } // Nun ist u_len > v_len.
          carry = muluadd_loop_lsp(q,v->LSDptr,u->LSDptr,v_len);
          if (!(carry==0))
            { var uintD* ptr = u->LSDptr lspop v_len;
              if ((lspref(ptr,0) += carry) < carry)
                { if (!( inc_loop_lsp(ptr lspop 1,u_len-v_len-1) ==0))
                    { lsprefnext(u->MSDptr) = 1; u->len++; }
            }   }
          while (mspref(u->MSDptr,0)==0) { msshrink(u->MSDptr); u->len--; } // normalisieren
    }   }

  // Bildet (u,v) := (x1*u+y1*v,x2*u+y2*v), wobei für u,v genügend Platz sei:
  // (Dabei sei u>0 oder v>0, nachher u>0 und v>0.)
  static void NUDS_likobi2_NUDS (DS* u, DS* v, partial_gcd_result* q, uintD* c_LSDptr, uintD* d_LSDptr)
    { var uintC u_len = u->len;
      var uintC v_len = v->len;
      var uintC c_len;
      var uintC d_len;
      if (u_len >= v_len)
        { mulu_loop_lsp(q->x1,u->LSDptr,c_LSDptr,u_len); c_len = u_len+1;
          mulu_loop_lsp(q->x2,u->LSDptr,d_LSDptr,u_len); d_len = u_len+1;
          if (!(v_len==0))
            {{var uintD carry =
                muluadd_loop_lsp(q->y1,v->LSDptr,c_LSDptr,v_len);
              if (!(carry==0))
                { var uintD* ptr = c_LSDptr lspop v_len;
                  if ((lspref(ptr,0) += carry) < carry)
                    { if (!( inc_loop_lsp(ptr lspop 1,u_len-v_len) ==0))
                        { lspref(c_LSDptr,c_len) = 1; c_len++; }
             }  }   }
             {var uintD carry =
                muluadd_loop_lsp(q->y2,v->LSDptr,d_LSDptr,v_len);
              if (!(carry==0))
                { var uintD* ptr = d_LSDptr lspop v_len;
                  if ((lspref(ptr,0) += carry) < carry)
                    { if (!(inc_loop_lsp(ptr lspop 1,u_len-v_len) ==0))
                        { lspref(d_LSDptr,d_len) = 1; d_len++; }
            }}  }   }
        }
        else
        { mulu_loop_lsp(q->y1,v->LSDptr,c_LSDptr,v_len); c_len = v_len+1;
          mulu_loop_lsp(q->y2,v->LSDptr,d_LSDptr,v_len); d_len = v_len+1;
          if (!(u_len==0))
            {{var uintD carry =
                muluadd_loop_lsp(q->x1,u->LSDptr,c_LSDptr,u_len);
              if (!(carry==0))
                { var uintD* ptr = c_LSDptr lspop u_len;
                  if ((lspref(ptr,0) += carry) < carry)
                    { if (!( inc_loop_lsp(ptr lspop 1,v_len-u_len) ==0))
                        { lspref(c_LSDptr,c_len) = 1; c_len++; }
             }  }   }
             {var uintD carry =
                muluadd_loop_lsp(q->x2,u->LSDptr,d_LSDptr,u_len);
              if (!(carry==0))
                { var uintD* ptr = d_LSDptr lspop u_len;
                  if ((lspref(ptr,0) += carry) < carry)
                    { if (!( inc_loop_lsp(ptr lspop 1,v_len-u_len) ==0))
                        { lspref(d_LSDptr,d_len) = 1; d_len++; }
            }}  }   }
        }
      u->MSDptr = copy_loop_lsp(c_LSDptr,u->LSDptr,c_len);
      while (mspref(u->MSDptr,0)==0) { msshrink(u->MSDptr); c_len--; }
      u->len = c_len;
      v->MSDptr = copy_loop_lsp(d_LSDptr,v->LSDptr,d_len);
      while (mspref(v->MSDptr,0)==0) { msshrink(v->MSDptr); d_len--; }
      v->len = d_len;
    }

  // Los geht's:
  const cl_I xgcd (const cl_I& a, const cl_I& b, cl_I* u, cl_I* v)
    { if (eq(a,1)) // a=1 -> g=1, (u,v)=(1,0)
        { *u = 1; *v = 0; return a; }
      if (eq(b,1)) // b=1 -> g=1, (u,v)=(0,1)
        { *u = 0; *v = 1; return b; }
      var sintL sA = (minusp(a) ? ~0 : 0); // Vorzeichen von A
      var sintL sB = (minusp(b) ? ~0 : 0); // Vorzeichen von B
      CL_ALLOCA_STACK;
      var uintD* a_MSDptr;
      var uintC a_len;
      var uintD* a_LSDptr;
      var uintD* b_MSDptr;
      var uintC b_len;
      var uintD* b_LSDptr;
      // Macro: erzeugt die NUDS zu (abs x), erniedrigt num_stack
      #define I_abs_to_NUDS(x,zero_statement)  \
        I_to_NDS_1(x, x##_MSDptr = , x##_len = , x##_LSDptr = ); /* (nichtleere) NDS holen */\
        if (x##_len == 0) { zero_statement } /* falls =0, fertig      */\
        if ((sintD)mspref(x##_MSDptr,0) < 0) /* falls <0, negieren:   */\
          { neg_loop_lsp(x##_LSDptr,x##_len); }				\
        if (mspref(x##_MSDptr,0) == 0) /* normalisieren (max. 1 Nulldigit entfernen) */\
          { msshrink(x##_MSDptr); x##_len--; }
      I_abs_to_NUDS(a, // (abs A) als NUDS erzeugen
                    // A=0 -> g=|B|, (u,v) = (0,sB)
                    { *u = 0; *v = (sB==0 ? cl_I(1) : cl_I(-1));
                      return abs(b);
                    });
      I_abs_to_NUDS(b, // (abs B) als NUDS erzeugen
                    // B=0 -> g=|A|, (u,v) = (sA,0)
                    { *u = (sA==0 ? cl_I(1) : cl_I(-1)); *v = 0;
                      return abs(a);
                    });
      // Jetzt ist a = a_MSDptr/a_len/a_LSDptr, b = b_MSDptr/b_len/b_LSDptr,
      // beides NUDS, und a_len>0, b_len>0.
      {// Beifaktoren:
       var DS uAa;
       var DS uBa;
       var DS uAb;
       var DS uBb;
       // Rechenregister:
       var uintD* divroomptr; // Platz für Divisionsergebnis
       var uintD* c_LSDptr;
       var uintD* d_LSDptr;
       // Platz für uAa,uBa,uAb,uBb besorgen:
       {var uintC u_len = b_len+1;
        num_stack_alloc(u_len,,uAa.LSDptr=); uAa.MSDptr = uAa.LSDptr;
        num_stack_alloc(u_len,,uAb.LSDptr=); uAb.MSDptr = uAb.LSDptr;
       }
       {var uintC u_len = a_len+1;
        num_stack_alloc(u_len,,uBa.LSDptr=); uBa.MSDptr = uBa.LSDptr;
        num_stack_alloc(u_len,,uBb.LSDptr=); uBb.MSDptr = uBb.LSDptr;
       }
       lsprefnext(uAa.MSDptr) = 1; uAa.len = 1; // uAa := 1
       uBa.len = 0; // uBa := 0
       uAb.len = 0; // uAb := 0
       lsprefnext(uBb.MSDptr) = 1; uBb.len = 1; // uBb := 1
       // Jetzt ist uAa = uAa.MSDptr/uAa.len/uAa.LSDptr,
       //           uBa = uBa.MSDptr/uBa.len/uBa.LSDptr,
       //           uAb = uAb.MSDptr/uAb.len/uAb.LSDptr,
       //           uBb = uBb.MSDptr/uBb.len/uBb.LSDptr,
       // alles NUDS.
       // Platz für zwei Rechenregister besorgen, mit je max(a_len,b_len)+1 Digits:
       {var uintC c_len = (a_len>=b_len ? a_len : b_len) + 1;
        num_stack_alloc(c_len,,c_LSDptr=);
        num_stack_alloc(c_len,divroomptr=,d_LSDptr=);
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
           a_greater_b_swap:
           swap(DS, uAa,uAb); // und uAa und uAb vertauschen
           swap(DS, uBa,uBb); // und uBa und uBb vertauschen
           sA = ~sA; sB = ~sB; // und sA und sB umdrehen
           a_greater_b:
           // Hier a>b>0, beides NUDS.
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
                   ? highlowDD(b_msd, (b_len==1 ? 0 : mspref(b_MSDptr,1)))
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
                  highlowDD(a_msd, (a_len==1 ? 0 : mspref(a_MSDptr,1)));
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
                       if (a_len>1)
                         { a_msd |= mspref(a_MSDptr,1) >> a_msd_size;
                           b_msd |= mspref(b_MSDptr,1) >> a_msd_size;
                     }   }
                   if (a_msd == b_msd) goto subtract;
                   #if DOUBLE_SPEED
                   if (a_len>1)
                     { a_nsd = mspref(a_MSDptr,1);
                       b_nsd = mspref(b_MSDptr,1);
                       if (shiftcount>0)
                         { a_nsd = a_nsd << shiftcount;
                           b_nsd = b_nsd << shiftcount;
                           if (a_len>2)
                             { a_nsd |= mspref(a_MSDptr,2) >> a_msd_size;
                               b_nsd |= mspref(b_MSDptr,2) >> a_msd_size;
                     }   }   }
                   else
                     { a_nsd = 0; b_nsd = 0; }
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
                   if (a_len>2)
                     { a_nsd |= mspref(a_MSDptr,2) >> a_msd_size;
                       b_nsd |= mspref(b_MSDptr,1) >> a_msd_size;
                     }
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
               #if HAVE_DD
               partial_gcd(highlowDD(a_msd,a_nsd),highlowDD(b_msd,b_nsd),&likobi); // liefert x1,y1,x2,y2
               #else
               partial_gcd(a_msd,a_nsd,b_msd,b_nsd,&likobi); // liefert x1,y1,x2,y2
               #endif
               #else
               partial_gcd(a_msd,b_msd,&likobi); // liefert x1,y1,x2,y2, aber nur halb so gut
               #endif
               // Hier y1>0.
               if (likobi.x2==0)
                 { // Ersetze (a,b) := (a-y1*b,b).
                   if (likobi.y1==1) goto subtract; // einfacherer Fall
                   // Dazu evtl. a um 1 Digit erweitern, so daß a_len=b_len+1:
                   if (a_len == b_len) { lsprefnext(a_MSDptr) = 0; a_len++; }
                   // und y1*b von a subtrahieren:
                   mspref(a_MSDptr,0) -= mulusub_loop_lsp(likobi.y1,b_LSDptr,a_LSDptr,b_len);
                   NUDS_likobi1_NUDS(&uAa,&uAb,likobi.y1); // uAa := uAa + y1 * uAb
                   NUDS_likobi1_NUDS(&uBa,&uBb,likobi.y1); // uBa := uBa + y1 * uBb
                 }
                 else
                 { // Ersetze (uAa,uAb) := (x1*uAa+y1*uAb,x2*uAa+y2*uAb) :
                   NUDS_likobi2_NUDS(&uAa,&uAb,&likobi,c_LSDptr,d_LSDptr);
                   // Ersetze (uBa,uBb) := (x1*uBa+y1*uBb,x2*uBa+y2*uBb) :
                   NUDS_likobi2_NUDS(&uBa,&uBb,&likobi,c_LSDptr,d_LSDptr);
                   // Ersetze (a,b) := (x1*a-y1*b,-x2*a+y2*b).
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
                 NUDS_likobi0_NUDS(&uAa,&uAb); // uAa := uAa + uAb
                 NUDS_likobi0_NUDS(&uBa,&uBb); // uBa := uBa + uBb
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
               b_len = r.len; if (b_len==0) goto return_a_coeffsb; // b=0 -> fertig
               b_LSDptr = old_a_LSDptr; // b übernimmt den vorherigen Platz von a
               b_MSDptr = copy_loop_lsp(r.LSDptr,b_LSDptr,b_len); // b := r kopieren
               // (uAa,uAb) := (uAb,uAa+q*uAb) :
               if (!(uAb.len==0))
                 { cl_UDS_mul(q.LSDptr,q.len,uAb.LSDptr,uAb.len,c_LSDptr); // q * uAb
                   var DS c;
                   c.LSDptr = c_LSDptr; c.len = q.len + uAb.len;
                   if (lspref(c_LSDptr,c.len-1)==0) { c.len--; } // normalisieren
                   NUDS_likobi0_NUDS(&uAa,&c); // zu uAa addieren
                 } // noch uAa,uAb vertauschen (später)
               // (uBa,uBb) := (uBb,uBa+q*uBb) :
               if (!(uBb.len==0))
                 { cl_UDS_mul(q.LSDptr,q.len,uBb.LSDptr,uBb.len,c_LSDptr); // q * uBb
                   var DS c;
                   c.LSDptr = c_LSDptr; c.len = q.len + uBb.len;
                   if (lspref(c_LSDptr,c.len-1)==0) { c.len--; } // normalisieren
                   NUDS_likobi0_NUDS(&uBa,&c); // zu uBa addieren
                 } // noch uBa,uBb vertauschen (später)
               goto a_greater_b_swap; // Nun ist a>b>0
             }}
         }
       // Nun ist a = b. Wähle diejenige der beiden Linearkombinationen
       //   a =  uAa*sA * A + -uBa*sB * B
       //   b = -uAb*sA * A +  uBb*sB * B
       // die die betragsmäßig kleinsten Koeffizienten hat.
       // Teste auf uBa < uBb. (Das kann auftreten, z.B. bei
       // A=560014183, B=312839871 wird a=b=1, uAa < uAb, uBa < uBb.)
       // Falls uBa = uBb, teste auf uAa < uAb. (Das kann auftreten, z.B. bei
       // A=2, B=3 wird a=b=1, uAa < uAb, uBa = uBb.)
       if (uBb.len > uBa.len) goto return_a_coeffsa;
       if (uBb.len < uBa.len) goto return_a_coeffsb;
       // (uBb.len == uBa.len)
       { var cl_signean vergleich = compare_loop_msp(uBb.MSDptr,uBa.MSDptr,uBb.len);
         if (vergleich > 0) goto return_a_coeffsa;
         if (vergleich < 0) goto return_a_coeffsb;
       }
       if (uAb.len > uAa.len) goto return_a_coeffsa;
       if (uAb.len < uAa.len) goto return_a_coeffsb;
       // (uAb.len == uAa.len)
       if (compare_loop_msp(uAb.MSDptr,uAa.MSDptr,uAb.len) > 0)
         return_a_coeffsa:
         { // uAa mit Vorfaktor sA versehen:
           lsprefnext(uAa.MSDptr) = 0; uAa.len++;
           if (!(sA==0)) { neg_loop_lsp(uAa.LSDptr,uAa.len); }
           // uBa mit Vorfaktor -sB versehen:
           lsprefnext(uBa.MSDptr) = 0; uBa.len++;
           if (sB==0) { neg_loop_lsp(uBa.LSDptr,uBa.len); }
           *u = DS_to_I(uAa.MSDptr,uAa.len); // DS uAa als Vorfaktor von A
           *v = DS_to_I(uBa.MSDptr,uBa.len); // DS uBa als Vorfaktor von B
         }
         else
         return_a_coeffsb:
         { // uAb mit Vorfaktor -sA versehen:
           lsprefnext(uAb.MSDptr) = 0; uAb.len++;
           if (sA==0) { neg_loop_lsp(uAb.LSDptr,uAb.len); }
           // uBb mit Vorfaktor sB versehen:
           lsprefnext(uBb.MSDptr) = 0; uBb.len++;
           if (!(sB==0)) { neg_loop_lsp(uBb.LSDptr,uBb.len); }
           *u = DS_to_I(uAb.MSDptr,uAb.len); // DS uAb als Vorfaktor von A
           *v = DS_to_I(uBb.MSDptr,uBb.len); // DS uBb als Vorfaktor von B
         }
      }
      return NUDS_to_I(a_MSDptr,a_len); // NUDS a als ggT
      #undef I_abs_to_NUDS
    }

#endif /* GCD_ALGO == 3 */

}  // namespace cln
