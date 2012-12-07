  static void mulu_karatsuba (const uintD* sourceptr1, uintC len1,
                              const uintD* sourceptr2, uintC len2,
                              uintD* destptr)
    // Karatsuba-Multiplikation
    // Prinzip: (x1*b^k+x0) * (y1*b^k+y0)
    //        = x1*y1 * b^2k + ((x1+x0)*(y1+y0)-x1*y1-x0*y0) * b^k + x0*y0
    // Methode 1 (Collins/Loos, Degel):
    // source2 wird in floor(len2/len1) einzelne UDS mit je einer
    // Länge len3 (len1 <= len3 < 2*len1) unterteilt,
    // jeweils k=floor(len3/2).
    // Methode 2 (Haible):
    // source2 wird in ceiling(len2/len1) einzelne UDS mit je einer
    // Länge len3 (0 < len3 <= len1) unterteilt, jeweils k=floor(len1/2).
    // Aufwand für die hinteren Einzelteile:
    // bei beiden Methoden jeweils 3*len1^2.
    // Aufwand für das vorderste Teil (alles, falls len1 <= len2 < 2*len1)
    // mit r = len1, s = (len2 mod len1) + len1 (>= len1, < 2*len1):
    // bei Methode 1:
    //                       |   :       |  r
    //                    |      :       |  s
    //      (r-s/2)*s/2 + s/2*s/2 + s/2*s/2 = r*s/2 + s^2/4 .
    // bei Methode 2:
    //                       |     :     |  r
    //                    |  |     :     |  s
    //      (s-r)*r + r/2*r/2 + r/2*r/2 + r/2*r/2 = r*s - r^2/4 .
    // Wegen (r*s/2 + s^2/4) - (r*s - r^2/4) = (r-s)^2/4 >= 0
    // ist Methode 2 günstiger.
    // Denkfehler! Dies gilt - wenn überhaupt - nur knapp oberhalb des
    // Break-Even-Points.
    // Im allgemeinen ist der Multiplikationsaufwand für zwei Zahlen der
    // Längen u bzw. v nämlich gegeben durch  min(u,v)^c * max(u,v),
    // wobei  c = log3/log2 - 1 = 0.585...
    // Dadurch wird der Aufwand in Abhängigkeit des Parameters t = k,
    // r/2 <= t <= s/2 (der einzig sinnvolle Bereich), zu
    // (r-t)^c*(s-t) + t^c*(s-t) + t^(1+c).
    // Dessen Optimum liegt (im Bereich r <= s <= 2*r)
    // - im klassischen Fall c=1 tatsächlich stets bei t=r/2 [Methode 2],
    // - im Karatsuba-Fall c=0.6 aber offenbar bei t=s/2 [Methode 1]
    //   oder ganz knapp darunter.
    // Auch erweist sich Methode 1 im Experiment als effizienter.
    // Daher implementieren wir Methode 1 :
    { // Es ist 2 <= len1 <= len2.
      // Spezialfall Quadrieren abfangen (häufig genug, daß sich das lohnt):
      if (sourceptr1 == sourceptr2)
        if (len1 == len2)
          { mulu_karatsuba_square(sourceptr1,len1,destptr); return; }
      var bool first_part = true; // Flag, ob jetzt das erste Teilprodukt berechnet wird
      if (len2 >= 2*len1)
        { CL_SMALL_ALLOCA_STACK;
          // Teilprodukte von jeweils len1 mal len1 Digits bilden:
          var uintC k_lo = floor(len1,2); // Länge der Low-Teile: floor(len1/2) >0
          var uintC k_hi = len1 - k_lo; // Länge der High-Teile: ceiling(len1/2) >0
          // Es gilt k_lo <= k_hi <= len1, k_lo + k_hi = len1.
          // Summe x1+x0 berechnen:
          var uintD* sum1_MSDptr;
          var uintC sum1_len = k_hi; // = max(k_lo,k_hi)
          var uintD* sum1_LSDptr;
          num_stack_small_alloc_1(sum1_len,sum1_MSDptr=,sum1_LSDptr=);
          {var uintD carry = // Hauptteile von x1 und x0 addieren:
             add_loop_lsp(sourceptr1 lspop k_lo,sourceptr1,sum1_LSDptr,k_lo);
           if (!(k_lo==k_hi))
             // noch k_hi-k_lo = 1 Digits abzulegen
             { mspref(sum1_MSDptr,0) = lspref(sourceptr1,len1-1); // = lspref(sourceptr1,2*k_lo)
               if (!(carry==0)) { if (++(mspref(sum1_MSDptr,0)) == 0) carry=1; else carry=0; }
             }
           if (carry) { lsprefnext(sum1_MSDptr) = 1; sum1_len++; }
          }
         {  // Platz für Summe y1+y0 belegen:
            var uintC sum2_maxlen = k_hi+1;
            var uintD* sum2_LSDptr;
            num_stack_small_alloc(sum2_maxlen,,sum2_LSDptr=);
            // Platz für Produkte x0*y0, x1*y1 belegen:
          { var uintD* prod_MSDptr;
            var uintD* prod_LSDptr;
            var uintD* prodhi_LSDptr;
            num_stack_small_alloc(2*len1,prod_MSDptr=,prod_LSDptr=);
            prodhi_LSDptr = prod_LSDptr lspop 2*k_lo;
            // prod_MSDptr/2*len1/prod_LSDptr wird zuerst die beiden
            // Produkte x1*y1 in prod_MSDptr/2*k_hi/prodhi_LSDptr
            //      und x0*y0 in prodhi_LSDptr/2*k_lo/prod_LSDptr,
            // dann das Produkt (b^k*x1+x0)*(b^k*y1+y0) enthalten.
            // Platz fürs Produkt (x1+x0)*(y1+y0) belegen:
           {var uintD* prodmid_MSDptr;
            var uintD* prodmid_LSDptr;
            num_stack_small_alloc(sum1_len+sum2_maxlen,prodmid_MSDptr=,prodmid_LSDptr=);
            // Schleife über die hinteren Einzelteile:
            do { // Produkt x0*y0 berechnen:
                 cl_UDS_mul(sourceptr1,k_lo,sourceptr2,k_lo,prod_LSDptr);
                 // Produkt x1*y1 berechnen:
                 cl_UDS_mul(sourceptr1 lspop k_lo,k_hi,sourceptr2 lspop k_lo,k_hi,prodhi_LSDptr);
                 // Summe y1+y0 berechnen:
                {var uintC sum2_len = k_hi; // = max(k_lo,k_hi)
                 var uintD* sum2_MSDptr = sum2_LSDptr lspop sum2_len;
                 {var uintD carry = // Hauptteile von y1 und y0 addieren:
                    add_loop_lsp(sourceptr2 lspop k_lo,sourceptr2,sum2_LSDptr,k_lo);
                  if (!(k_lo==k_hi))
                    // noch k_hi-k_lo = 1 Digits abzulegen
                    { mspref(sum2_MSDptr,0) = lspref(sourceptr2,len1-1); // = lspref(sourceptr2,2*k_lo)
                      if (!(carry==0)) { if (++(mspref(sum2_MSDptr,0)) == 0) carry=1; else carry=0; }
                    }
                  if (carry) { lsprefnext(sum2_MSDptr) = 1; sum2_len++; }
                 }
                 // Produkt (x1+x0)*(y1+y0) berechnen:
                 cl_UDS_mul(sum1_LSDptr,sum1_len,sum2_LSDptr,sum2_len,prodmid_LSDptr);
                 // Das Produkt beansprucht  2*k_hi + (0 oder 1) <= sum1_len + sum2_len  Digits.
                 {var uintC prodmid_len = sum1_len+sum2_len;
                  // Davon x1*y1 abziehen:
                  {var uintD carry =
                     subfrom_loop_lsp(prodhi_LSDptr,prodmid_LSDptr,2*k_hi);
                   // Falls Carry: Produkt beansprucht 2*k_hi+1 Digits.
                   // Carry um maximal 1 Digit weitertragen:
                   if (!(carry==0)) { lspref(prodmid_LSDptr,2*k_hi) -= 1; }
                  }
                  // Und x0*y0 abziehen:
                  {var uintD carry =
                     subfrom_loop_lsp(prod_LSDptr,prodmid_LSDptr,2*k_lo);
                   // Carry um maximal prodmid_len-2*k_lo Digits weitertragen:
                   if (!(carry==0))
                     { dec_loop_lsp(prodmid_LSDptr lspop 2*k_lo,prodmid_len-2*k_lo); }
                  }
                  // prodmid_LSDptr[-prodmid_len..-1] enthält nun x0*y1+x1*y0.
                  // Dies wird zu prod = x1*y1*b^(2*k) + x0*y0 addiert:
                  {var uintD carry =
                     addto_loop_lsp(prodmid_LSDptr,prod_LSDptr lspop k_lo,prodmid_len);
                     // (Benutze dabei k_lo+prodmid_len <= k_lo+2*(k_hi+1) = 2*len1-k_lo+2 <= 2*len1 .)
                   if (!(carry==0))
                     { inc_loop_lsp(prod_LSDptr lspop (k_lo+prodmid_len),2*len1-(k_lo+prodmid_len)); }
                }}}
                 // Das Teilprodukt zum Gesamtprodukt addieren:
                 if (first_part)
                   { copy_loop_lsp(prod_LSDptr,destptr,2*len1);
                     destptr = destptr lspop len1;
                     first_part = false;
                   }
                   else
                   { var uintD carry =
                       addto_loop_lsp(prod_LSDptr,destptr,len1);
                     destptr = destptr lspop len1;
                     copy_loop_lsp(prod_LSDptr lspop len1,destptr,len1);
                     if (!(carry==0)) { inc_loop_lsp(destptr,len1); }
                   }
                 sourceptr2 = sourceptr2 lspop len1; len2 -= len1;
               }
               while (len2 >= 2*len1);
         }}}
        }
      // Nun ist len1 <= len2 < 2*len1.
      // letztes Teilprodukt von len1 mal len2 Digits bilden:
     {CL_SMALL_ALLOCA_STACK;
      var uintD* prod_MSDptr;
      var uintC prod_len = len1+len2;
      var uintD* prod_LSDptr;
      num_stack_small_alloc(prod_len,prod_MSDptr=,prod_LSDptr=);
      { var uintC k_hi = floor(len2,2); // Länge der High-Teile: floor(len2/2) >0
        var uintC k_lo = len2 - k_hi; // Länge der Low-Teile: ceiling(len2/2) >0
        // Es gilt k_hi <= k_lo <= len1 <= len2, k_lo + k_hi = len2.
        var uintC x1_len = len1-k_lo; // <= len2-k_lo = k_hi <= k_lo
        // Summe x1+x0 berechnen:
        var uintD* sum1_MSDptr;
        var uintC sum1_len = k_lo; // = max(k_lo,k_hi)
        var uintD* sum1_LSDptr;
        num_stack_small_alloc_1(sum1_len,sum1_MSDptr=,sum1_LSDptr=);
        {var uintD carry = // x1 und unteren Teil von x0 addieren:
           add_loop_lsp(sourceptr1 lspop k_lo,sourceptr1,sum1_LSDptr,x1_len);
         // und den oberen Teil von x0 dazu:
         copy_loop_lsp(sourceptr1 lspop x1_len,sum1_LSDptr lspop x1_len,k_lo-x1_len);
         if (!(carry==0))
           { carry = inc_loop_lsp(sum1_LSDptr lspop x1_len,k_lo-x1_len);
             if (carry) { lsprefnext(sum1_MSDptr) = 1; sum1_len++; }
           }
        }
       {// Summe y1+y0 berechnen:
        var uintD* sum2_MSDptr;
        var uintC sum2_len = k_lo; // = max(k_lo,k_hi)
        var uintD* sum2_LSDptr;
        num_stack_small_alloc_1(sum2_len,sum2_MSDptr=,sum2_LSDptr=);
        {var uintD carry = // Hauptteile von y1 und y0 addieren:
           add_loop_lsp(sourceptr2 lspop k_lo,sourceptr2,sum2_LSDptr,k_hi);
         if (!(k_lo==k_hi))
           // noch k_lo-k_hi = 1 Digits abzulegen
           { mspref(sum2_MSDptr,0) = lspref(sourceptr2,k_lo-1); // = lspref(sourceptr2,k_hi)
             if (!(carry==0)) { if (++(mspref(sum2_MSDptr,0)) == 0) carry=1; else carry=0; }
           }
         if (carry) { lsprefnext(sum2_MSDptr) = 1; sum2_len++; }
        }
        // Platz für Produkte x0*y0, x1*y1:
        { var uintC prodhi_len = x1_len+k_hi;
          var uintD* prodhi_LSDptr = prod_LSDptr lspop 2*k_lo;
          // prod_MSDptr/len1+len2/prod_LSDptr wird zuerst die beiden
          // Produkte x1*y1 in prod_MSDptr/x1_len+k_hi/prodhi_LSDptr
          //      und x0*y0 in prodhi_LSDptr/2*k_lo/prod_LSDptr,
          // dann das Produkt (b^k*x1+x0)*(b^k*y1+y0) enthalten.
          // Platz fürs Produkt (x1+x0)*(y1+y0) belegen:
         {var uintD* prodmid_MSDptr;
          var uintC prodmid_len = sum1_len+sum2_len;
          var uintD* prodmid_LSDptr;
          num_stack_small_alloc(prodmid_len,prodmid_MSDptr=,prodmid_LSDptr=);
          // Produkt (x1+x0)*(y1+y0) berechnen:
          cl_UDS_mul(sum1_LSDptr,sum1_len,sum2_LSDptr,sum2_len,prodmid_LSDptr);
          // Das Produkt beansprucht  2*k_lo + (0 oder 1) <= sum1_len + sum2_len = prodmid_len  Digits.
          // Produkt x0*y0 berechnen:
          cl_UDS_mul(sourceptr1,k_lo,sourceptr2,k_lo,prod_LSDptr);
          // Produkt x1*y1 berechnen:
          if (!(x1_len==0))
            { cl_UDS_mul(sourceptr1 lspop k_lo,x1_len,sourceptr2 lspop k_lo,k_hi,prodhi_LSDptr);
             // Und x1*y1 abziehen:
             {var uintD carry =
                subfrom_loop_lsp(prodhi_LSDptr,prodmid_LSDptr,prodhi_len);
              // Carry um maximal prodmid_len-prodhi_len Digits weitertragen:
              if (!(carry==0))
                { dec_loop_lsp(prodmid_LSDptr lspop prodhi_len,prodmid_len-prodhi_len); }
            }}
            else
            // Produkt x1*y1=0, nichts abzuziehen
            { clear_loop_lsp(prodhi_LSDptr,prodhi_len); }
          // Und x0*y0 abziehen:
          {var uintD carry =
             subfrom_loop_lsp(prod_LSDptr,prodmid_LSDptr,2*k_lo);
           // Falls Carry: Produkt beansprucht 2*k_lo+1 Digits.
           // Carry um maximal 1 Digit weitertragen:
           if (!(carry==0)) { lspref(prodmid_LSDptr,2*k_lo) -= 1; }
          }
          // prodmid_LSDptr[-prodmid_len..-1] enthält nun x0*y1+x1*y0.
          // Dies ist < b^k_lo * b^k_hi + b^x1_len * b^k_lo
          //          = b^len2 + b^len1 <= 2 * b^len2,
          // paßt also in len2+1 Digits.
          // Im Fall x1_len=0 ist es sogar < b^k_lo * b^k_hi = b^len2,
          // es paßt also in len2 Digits.
          // prodmid_len, wenn möglich, um maximal 2 verkleinern:
          // (benutzt prodmid_len >= 2*k_lo >= len2 >= 2)
          if (mspref(prodmid_MSDptr,0)==0)
            { prodmid_len--;
              if (mspref(prodmid_MSDptr,1)==0) { prodmid_len--; }
            }
          // Nun ist k_lo+prodmid_len <= len1+len2 .
          // (Denn es war prodmid_len = sum1_len+sum2_len <= 2*(k_lo+1)
          //  <= len2+3, und nach 2-maliger Verkleinerung jedenfalls
          //  prodmid_len <= len2+1. Im Falle k_lo < len1 also
          //  k_lo + prodmid_len <= (len1-1)+(len2+1) = len1+len2.
          //  Im Falle k_lo = len1 aber ist x1_len=0, sum1_len = k_lo, also
          //  war prodmid_len = sum1_len+sum2_len <= 2*k_lo+1 <= len2+2,
          //  nach 2-maliger Verkleinerung jedenfalls prodmid_len <= len2.)
          // prodmid*b^k = (x0*y1+x1*y0)*b^k zu prod = x1*y1*b^(2*k) + x0*y0 addieren:
          {var uintD carry =
             addto_loop_lsp(prodmid_LSDptr,prod_LSDptr lspop k_lo,prodmid_len);
           if (!(carry==0))
             { inc_loop_lsp(prod_LSDptr lspop (k_lo+prodmid_len),prod_len-(k_lo+prodmid_len)); }
      }}}}}
      // Das Teilprodukt zum Gesamtprodukt addieren:
      if (first_part)
        { copy_loop_lsp(prod_LSDptr,destptr,prod_len); }
        else
        { var uintD carry =
            addto_loop_lsp(prod_LSDptr,destptr,len1);
          destptr = destptr lspop len1;
          copy_loop_lsp(prod_LSDptr lspop len1,destptr,len2);
          if (!(carry==0)) { inc_loop_lsp(destptr,len2); }
        }
    }}
