  // Eine vereinfachte Version von mulu_karatsuba für den Fall
  // sourceptr1 == sourceptr2 && len1 == len2.
  // Weniger Variablen, eine Additionsschleife weniger, eine Kopierschleife
  // weniger, und bei rekursiven Aufrufen ist wieder
  // sourceptr1 == sourceptr2 && len1 == len2.
  static void mulu_karatsuba_square (const uintD* sourceptr, uintC len,
                                     uintD* destptr)
    { // Es ist 2 <= len.
      CL_SMALL_ALLOCA_STACK;
      var uintC prod_len = 2*len;
      var uintD* prod_LSDptr = destptr;
      var uintC k_hi = floor(len,2); // Länge der High-Teile: floor(len/2) >0
      var uintC k_lo = len - k_hi; // Länge der Low-Teile: ceiling(len/2) >0
      // Es gilt k_hi <= k_lo <= len, k_lo + k_hi = len.
      // Summe x1+x0 berechnen:
      var uintD* sum_MSDptr;
      var uintC sum_len = k_lo; // = max(k_lo,k_hi)
      var uintD* sum_LSDptr;
      num_stack_small_alloc_1(sum_len,sum_MSDptr=,sum_LSDptr=);
      {var uintD carry = // Hauptteile von x1 und x0 addieren:
         add_loop_lsp(sourceptr lspop k_lo,sourceptr,sum_LSDptr,k_hi);
       if (!(k_lo==k_hi))
         // noch k_lo-k_hi = 1 Digits abzulegen
         { mspref(sum_MSDptr,0) = lspref(sourceptr,k_lo-1); // = lspref(sourceptr,k_hi)
           if (!(carry==0)) { if (++(mspref(sum_MSDptr,0)) == 0) carry=1; else carry=0; }
         }
       if (carry) { lsprefnext(sum_MSDptr) = 1; sum_len++; }
      }
      // Platz für Produkte x0*x0, x1*x1:
      { var uintC prodhi_len = 2*k_hi;
        var uintD* prodhi_LSDptr = prod_LSDptr lspop 2*k_lo;
        // prod_MSDptr/2*len/prod_LSDptr wird zuerst die beiden
        // Produkte x1*x1 in prod_MSDptr/2*k_hi/prodhi_LSDptr
        //      und x0*x0 in prodhi_LSDptr/2*k_lo/prod_LSDptr,
        // dann das Produkt (b^k*x1+x0)*(b^k*x1+x0) enthalten.
        // Platz fürs Produkt (x1+x0)*(x1+x0) belegen:
       {var uintD* prodmid_MSDptr;
        var uintC prodmid_len = 2*sum_len;
        var uintD* prodmid_LSDptr;
        num_stack_small_alloc(prodmid_len,prodmid_MSDptr=,prodmid_LSDptr=);
        // Produkt (x1+x0)*(x1+x0) berechnen:
        cl_UDS_mul_square(sum_LSDptr,sum_len,prodmid_LSDptr);
        // Das Produkt beansprucht  2*k_lo + (0 oder 1) <= 2*sum_len = prodmid_len  Digits.
        // Produkt x0*x0 berechnen:
        cl_UDS_mul_square(sourceptr,k_lo,prod_LSDptr);
        // Produkt x1*x1 berechnen:
        cl_UDS_mul_square(sourceptr lspop k_lo,k_hi,prodhi_LSDptr);
        // Und x1*x1 abziehen:
        {var uintD carry =
           subfrom_loop_lsp(prodhi_LSDptr,prodmid_LSDptr,prodhi_len);
         // Carry um maximal prodmid_len-prodhi_len Digits weitertragen:
         if (!(carry==0))
           { dec_loop_lsp(prodmid_LSDptr lspop prodhi_len,prodmid_len-prodhi_len); }
        }
        // Und x0*x0 abziehen:
        {var uintD carry =
           subfrom_loop_lsp(prod_LSDptr,prodmid_LSDptr,2*k_lo);
         // Falls Carry: Produkt beansprucht 2*k_lo+1 Digits.
         // Carry um maximal 1 Digit weitertragen:
         if (!(carry==0)) { lspref(prodmid_LSDptr,2*k_lo) -= 1; }
        }
        // prodmid_LSDptr[-prodmid_len..-1] enthält nun 2*x0*x1.
        // Dies ist < 2 * b^k_lo * b^k_hi = 2 * b^len,
        // paßt also in len+1 Digits.
        // prodmid_len, wenn möglich, um maximal 2 verkleinern:
        // (benutzt prodmid_len >= 2*k_lo >= len >= 2)
        if (mspref(prodmid_MSDptr,0)==0)
          { prodmid_len--;
            if (mspref(prodmid_MSDptr,1)==0) { prodmid_len--; }
          }
        // Nun ist k_lo+prodmid_len <= 2*len .
        // (Denn es war prodmid_len = 2*sum_len <= 2*(k_lo+1)
        //  <= len+3, und nach 2-maliger Verkleinerung jedenfalls
        //  prodmid_len <= len+1. Wegen k_lo < len also
        //  k_lo + prodmid_len <= (len-1)+(len+1) = 2*len.)
        // prodmid*b^k = 2*x0*x1*b^k zu prod = x1*x1*b^(2*k) + x0*x0 addieren:
        {var uintD carry =
           addto_loop_lsp(prodmid_LSDptr,prod_LSDptr lspop k_lo,prodmid_len);
         if (!(carry==0))
           { inc_loop_lsp(prod_LSDptr lspop (k_lo+prodmid_len),prod_len-(k_lo+prodmid_len)); }
    } }}}
