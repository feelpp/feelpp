// cl_divide().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/exception.h"

namespace cln {

// Dividiert zwei Integers x,y >=0 und liefert Quotient und Rest
// der Division x/y. Bei y=0 Error.
// cl_divide(x,y)
// > x,y: Integers >=0
// < q,r: Quotient q, Rest r
  const cl_I_div_t cl_divide (const cl_I& x, const cl_I& y)
    { if (fixnump(x))
        // x Fixnum >=0
        { if (fixnump(y))
            // auch y Fixnum >=0
            { var uintV x_ = FN_to_UV(x);
              var uintV y_ = FN_to_UV(y);
              if (y_==0) { throw division_by_0_exception(); }
              elif (x_ < y_)
                // Trivialfall: q=0, r=x
                goto trivial;
              else
                {
                  #if (intVsize>32)
                  if (x_ >= bit(32))
                    { if (y_ < bit(32))
                        // 64-durch-32-Bit-Division
                        { var uint64 q;
                          var uint32 r;
                          divu_6432_6432(x_,y_,q=,r=);
                          return cl_I_div_t(
                            /* result.quotient =  */ UQ_to_I(q),
                            /* result.remainder = */ UL_to_I(r)
                            );
                        }
                      else
                        // volle 64-durch-64-Bit-Division
                        { var uint64 q;
                          var uint64 r;
                          divu_6464_6464(x_,y_,q=,r=);
                          return cl_I_div_t(
                            /* result.quotient =  */ UQ_to_I(q),
                            /* result.remainder = */ UQ_to_I(r)
                            );
                        }
                    }
                  else
                  #endif
                    { if (y_ < bit(16))
                        // 32-durch-16-Bit-Division
                        { var uint32 q;
                          var uint16 r;
                          divu_3216_3216(x_,y_,q=,r=);
                          return cl_I_div_t(
                            /* result.quotient =  */ UL_to_I(q),
                            /* result.remainder = */ L_to_FN((uintL)r)
                            );
                        }
                      else
                        // volle 32-durch-32-Bit-Division
                        { var uint32 q;
                          var uint32 r;
                          divu_3232_3232(x_,y_,q=,r=);
                          return cl_I_div_t(
                            /* result.quotient =  */ UL_to_I(q),
                            /* result.remainder = */ UL_to_I(r)
                            );
                        }
                    }
                }
            }
            else
            // y Bignum >0
            { trivial:
              // Trivialfall: q=0, r=x
              return cl_I_div_t(
                /* result.quotient =  */ 0,
                /* result.remainder = */ x
                );
            }
        }
        else
        // x Bignum -> allgemeine Division:
        { CL_ALLOCA_STACK;
          var const uintD* x_MSDptr;
          var uintC x_len;
          var const uintD* x_LSDptr;
          var const uintD* y_MSDptr;
          var uintC y_len;
          var const uintD* y_LSDptr;
          // x in NDS umwandeln, als UDS auffassen:
          BN_to_NDS_nocopy(x, x_MSDptr=,x_len=,x_LSDptr=);
          // y in NDS umwandeln, als UDS auffassen:
          I_to_NDS_nocopy(y, y_MSDptr=,y_len=,y_LSDptr=,/*true*/false,);
          // dividieren:
         {var DS q;
          var DS r;
          UDS_divide(x_MSDptr,x_len,x_LSDptr,y_MSDptr,y_len,y_LSDptr, &q,&r);
          // q in Integer umwandeln:
          var cl_I quotient = NUDS_to_I(q.MSDptr,q.len);
          // r in Integer umwandeln (jetzt erst, nachdem q verwertet ist!):
          return cl_I_div_t(
            quotient,
            /* result.remainder = */ NUDS_to_I(r.MSDptr,r.len)
            );
        }}
    }
// Bit complexity (N = length(x)): O(M(N)).

}  // namespace cln
