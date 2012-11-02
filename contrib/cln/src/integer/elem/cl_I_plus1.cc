// plus1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_I plus1 (const cl_I& x)
{
	if (fixnump(x))
	  { // x ist Fixnum
	    if (x.word != cl_combine(cl_FN_tag,bit(cl_value_len-1)-1))
		// bleibt Fixnum: direkt 1 addieren
		// This assumes cl_value_shift + cl_value_len == cl_pointer_size.
		{ return cl_I_from_word(x.word + cl_combine(0,1)); }
          }
        // die sichere Methode
        { CL_ALLOCA_STACK;
          var uintD* MSDptr;
          var uintC len;
          var uintD* LSDptr;
          I_to_NDS_1(x, MSDptr=,len=,LSDptr=); // NDS zu x bilden.
          DS_1_plus(LSDptr,len); // zur NDS 1 addieren
          return DS_to_I(MSDptr,len); // wieder zum Integer machen
        }
}

}  // namespace cln
