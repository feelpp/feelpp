// cl_LF_to_I().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

#include "base/cl_inline.h"
#include "float/lfloat/elem/cl_LF_minusp.cc"

namespace cln {

const cl_I cl_LF_to_I (const cl_LF& x)
{
// Methode:
// Falls x=0.0, Ergebnis 0.
// Sonst (ASH Vorzeichen*Mantisse (e-16n)).
      var uintE uexp = TheLfloat(x)->expo;
      if (uexp==0) { return 0; } // x=0.0 -> Ergebnis 0
      // Mantisse zu einem Integer machen:
      CL_ALLOCA_STACK;
      var uintD* MSDptr;
      var uintD* LSDptr;
      var uintC len = TheLfloat(x)->len;
      var uintC len1 = len+1; // brauche 1 Digit mehr
      num_stack_alloc(len1, MSDptr=,LSDptr=);
      copy_loop_msp(arrayMSDptr(TheLfloat(x)->data,len),MSDptr mspop 1,len); // Mantisse kopieren
      mspref(MSDptr,0) = 0; // und zusätzliches Nulldigit
      // Mantisse ist die UDS MSDptr/len1/LSDptr.
      if (minusp_inline(x))
        // x<0 -> Mantisse negieren:
        { neg_loop_lsp(LSDptr,len1); }
      // Vorzeichen*Mantisse ist die DS MSDptr/len1/LSDptr.
      // (ASH Vorzeichen*Mantisse (- e 16n)) durchführen:
      return ash(DS_to_I(MSDptr,len1),
                 minus(uexp, LF_exp_mid + intDsize*len)
                );
}

}  // namespace cln
