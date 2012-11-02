// cl_DF_to_I().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/dfloat/cl_DF.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

const cl_I cl_DF_to_I (const cl_DF& x)
{
// Methode:
// Falls x=0.0, Ergebnis 0.
// Sonst (ASH Vorzeichen*Mantisse (e-53)).
#if (cl_word_size==64)
      // x entpacken:
      var cl_signean sign;
      var sintL exp;
      var sint64 mant;
      DF_decode(x, { return 0; }, sign=,exp=,mant=);
      exp = exp-(DF_mant_len+1);
      // mant mit Vorzeichen versehen:
      if (!(sign==0)) { mant = -mant; }
      // in ein Bignum umwandeln und shiften:
      return ash( Q_to_I(mant), exp );
#else
      // x entpacken:
      var cl_signean sign;
      var sintL exp;
      var sint32 manthi;
      var uint32 mantlo;
      DF_decode2(x, { return 0; }, sign=,exp=,manthi=,mantlo=);
      exp = exp-(DF_mant_len+1);
      // mant mit Vorzeichen versehen:
      if (!(sign==0))
        { manthi = -manthi; mantlo = -mantlo; if (!(mantlo==0)) { manthi -= 1; } }
      // in ein Bignum umwandeln und shiften:
      return ash( L2_to_I(manthi,mantlo), exp );
#endif
}

}  // namespace cln
