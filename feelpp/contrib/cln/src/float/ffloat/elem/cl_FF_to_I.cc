// cl_FF_to_I().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/ffloat/cl_FF.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

const cl_I cl_FF_to_I (const cl_FF& x)
{
// Methode:
// Falls x=0.0, Ergebnis 0.
// Sonst (ASH Vorzeichen*Mantisse (e-24)).
      // x entpacken:
      var cl_signean sign;
      var sintL exp;
      var uint32 mant;
      FF_decode(x, { return 0; }, sign=,exp=,mant=);
      exp = exp-(FF_mant_len+1);
      return ash(
        // mant >0, <2^(FF_mant_len+1) in ein Fixnum umwandeln:
        #if (FF_mant_len+1 < cl_value_len)
          L_to_FN(sign==0 ? (sintL)mant : -(sintL)mant)
        #else
          L_to_I(sign==0 ? (sintL)mant : -(sintL)mant)
        #endif
        ,exp
        );
}

}  // namespace cln
