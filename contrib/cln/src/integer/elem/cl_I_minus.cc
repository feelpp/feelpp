// binary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_I operator- (const cl_I& x, const cl_I& y)
{
    // Methode:
    // x Fixnum ->
    //   y Fixnum -> beide direkt subtrahieren, mit L_to_I beenden
    //   y Bignum -> falls x=0, (- y); sonst beide zu DS machen, subtrahieren.
    // x Bignum ->
    //   y Fixnum -> falls y=0, x; sonst beide zu DS machen, subtrahieren.
    //   y Bignum -> beide zu DS machen, subtrahieren.
      var uintD* MSDptr;
      var uintC len;
      var uintD* LSDptr;
      // MSDptr/len/LSDptr bilden die DS des Ergebnisses.
      if (fixnump(x))
        { // x ist Fixnum
          if (fixnump(y))
            { // x,y sind Fixnums
              #if (cl_value_len < intVsize)
              return V_to_I( FN_to_V(x) - FN_to_V(y) ); // als intVsize-Bit-Zahlen subtrahieren
              #elif (cl_word_size==64)
              return Q_to_I( FN_to_Q(x) - FN_to_Q(y) ); // als 64-Bit-Zahlen subtrahieren
              #elif (intVsize==32) // && (cl_value_len == intVsize)
              var sint32 xhi = sign_of(FN_to_V(x));
              var uint32 xlo = FN_to_V(x);
              var sint32 yhi = sign_of(FN_to_V(y));
              var uint32 ylo = FN_to_V(y);
              xhi -= yhi;
              if (xlo < ylo) { xhi -= 1; }
              xlo -= ylo;
              return L2_to_I(xhi,xlo);
              #endif
            }
            else
            { // x ist Fixnum, y ist Bignum, also y länger
              #if (intDsize==64)
              var sint64 x_ = FN_to_V(x); // Wert von x
              #else
              var sintV x_ = FN_to_V(x); // Wert von x
              #endif
              if (FN_V_zerop(x,x_)) { return -y; } // bei x=0 Ergebnis (- y)
              CL_ALLOCA_STACK;
              BN_to_NDS_1(y, MSDptr=,len=,LSDptr=); // NDS zu y bilden.
              // vorsorglich 1 Digit mehr belegen:
              { var sintD sign = sign_of_sintD(mspref(MSDptr,0));
                lsprefnext(MSDptr) = sign; len++;
              }
              // Negierschleife:
              neg_loop_lsp(LSDptr,len);
              // MSDigit ist nun = 0x0000 oder = 0xFFFF
              // x_ zu den oberen pFN_maxlength Digits von -y addieren:
              {
                #if (intDsize==64)
                var uint64 y_ = lspref(LSDptr,0);
                var uint64 y_new = y_+(uint64)x_;
                lspref(LSDptr,0) = y_new;
                #else
                var uintV y_ = pFN_maxlength_digits_at(LSDptr);
                var uintV y_new = y_+(uintV)x_;
                set_pFN_maxlength_digits_at(LSDptr,y_new);
                #endif
                var uintD* midptr = LSDptr lspop pFN_maxlength;
                if (y_new < y_)
                  { // Carry.
                    if (!FN_V_minusp(x,x_)) // kürzerer Summand war positiv
                      // Dann ist ein positiver Übertrag weiterzutragen
                      // (Beispiel: 0002FFFC + 0007 = 00030003)
                      { DS_1_plus(midptr,len-pFN_maxlength); }
                  }
                  else
                  { // Kein Carry.
                    if (FN_V_minusp(x,x_)) // kürzerer Summand war negativ
                      // Dann ist ein negativer Übertrag weiterzutragen
                      // (Beispiel: 00020003 + FFF5 = 0001FFF8)
                      { DS_minus1_plus(midptr,len-pFN_maxlength); }
              }   }
              return DS_to_I(MSDptr,len); // DS wieder zum Integer machen
            }
        }
        else
        { // x ist Bignum
          if (fixnump(y))
            { // x ist Bignum, y ist Fixnum, also x länger
              #if (intDsize==64)
              var sint64 y_ = FN_to_V(y); // Wert von y
              #else
              var sintV y_ = FN_to_V(y); // Wert von y
              #endif
              if (FN_V_zerop(y,y_)) { return x; } // bei y=0 Ergebnis x
              CL_ALLOCA_STACK;
              BN_to_NDS_1(x, MSDptr=,len=,LSDptr=); // NDS zu x bilden.
              // len>=bn_minlength. len>pFN_maxlength erzwingen:
              if ((bn_minlength==pFN_maxlength) && (len==pFN_maxlength))
                { var sintD sign = sign_of_sintD(mspref(MSDptr,0));
                  lsprefnext(MSDptr) = sign; len++;
                }
              // y_ von den oberen pFN_maxlength Digits von x subtrahieren:
              {
                #if (intDsize==64)
                var uint64 x_ = lspref(LSDptr,0);
                var uint64 x_new = x_-(uint64)y_;
                lspref(LSDptr,0) = x_new;
                #else
                var uintV x_ = pFN_maxlength_digits_at(LSDptr);
                var uintV x_new = x_-(uintV)y_;
                set_pFN_maxlength_digits_at(LSDptr,x_new);
                #endif
                var uintD* midptr = LSDptr lspop pFN_maxlength;
                if (x_new > x_)
                  { // Carry.
                    if (!FN_V_minusp(y,y_)) // kürzerer Summand war positiv
                      // Dann ist ein negativer Übertrag weiterzutragen
                      // (Beispiel: 00030003 - 0007 = 0002FFFC)
                      { DS_minus1_plus(midptr,len-pFN_maxlength); }
                  }
                  else
                  { // Kein Carry.
                    if (FN_V_minusp(y,y_)) // kürzerer Summand war negativ
                      // Dann ist ein positiver Übertrag weiterzutragen
                      // (Beispiel: 0002FFF8 - FFF5 = 00030003)
                      { DS_1_plus(midptr,len-pFN_maxlength); }
              }   }
              return DS_to_I(MSDptr,len); // DS wieder zum Integer machen
            }
            else
            { // x und y sind Bignums
              if (TheBignum(x)->length > TheBignum(y)->length)
                { // x das längere von beiden.
                  CL_ALLOCA_STACK;
                  BN_to_NDS_1(x, MSDptr=,len=,LSDptr=); // NDS zu x bilden.
                  var const uintD* yMSDptr;
                  var uintC ylen;
                  var const uintD* yLSDptr;
                  BN_to_NDS_nocopy(y, yMSDptr=,ylen=,yLSDptr=); // NDS zu y bilden.
                  // yMSDptr/ylen/yLSDptr bilden die DS des kürzeren Arguments y.
                  // Es ist len>ylen.
                  // subtrahieren:
                  { var uintD* midptr = LSDptr lspop ylen;
                    var uintD carry = subfrom_loop_lsp(yLSDptr,LSDptr,ylen);
                    if (carry)
                      { // Carry.
                        if ((sintD)mspref(yMSDptr,0) >=0) // kürzerer Summand war positiv
                          // Dann ist ein negativer Übertrag weiterzutragen
                          // (Beispiel: 00030003 - 0007 = 0002FFFC)
                          { DS_minus1_plus(midptr,len-ylen); }
                      }
                      else
                      { // Kein Carry.
                        if ((sintD)mspref(yMSDptr,0) <0) // kürzerer Summand war negativ
                          // Dann ist ein positiver Übertrag weiterzutragen
                          // (Beispiel: 0002FFF8 - FFF5 = 00030003)
                          { DS_1_plus(midptr,len-ylen); }
                  }   }
                  return DS_to_I(MSDptr,len); // DS wieder zum Integer machen
                }
                else
                { // y das längere von beiden.
                  CL_ALLOCA_STACK;
                  BN_to_NDS_1(y, MSDptr=,len=,LSDptr=); // NDS zu y bilden.
                  // vorsorglich 1 Digit mehr belegen:
                  { var sintD sign = sign_of_sintD(mspref(MSDptr,0));
                    lsprefnext(MSDptr) = sign; len++;
                  }
                  // Negierschleife:
                  neg_loop_lsp(LSDptr,len);
                  // MSDigit ist nun = 0x0000 oder = 0xFFFF
                  var const uintD* xMSDptr;
                  var uintC xlen;
                  var const uintD* xLSDptr;
                  BN_to_NDS_nocopy(x, xMSDptr=,xlen=,xLSDptr=); // NDS zu x bilden.
                  // xMSDptr/xlen/xLSDptr bilden die DS des kürzeren Arguments x.
                  // Es ist jetzt len>xlen.
                  // addieren:
                  { var uintD* midptr = LSDptr lspop xlen;
                    var uintD carry = addto_loop_lsp(xLSDptr,LSDptr,xlen);
                    if (carry)
                      { // Carry.
                        if ((sintD)mspref(xMSDptr,0) >=0) // kürzerer Summand war positiv
                          // Dann ist ein positiver Übertrag weiterzutragen
                          // (Beispiel: 0002FFFC + 0007 = 00030003)
                          { DS_1_plus(midptr,len-xlen); }
                      }
                      else
                      { // Kein Carry.
                        if ((sintD)mspref(xMSDptr,0) <0) // kürzerer Summand war negativ
                          // Dann ist ein negativer Übertrag weiterzutragen
                          // (Beispiel: 00020003 + FFF5 = 0001FFF8)
                          { DS_minus1_plus(midptr,len-xlen); }
                  }   }
                  return DS_to_I(MSDptr,len); // DS wieder zum Integer machen
                }
        }   }
}

}  // namespace cln
