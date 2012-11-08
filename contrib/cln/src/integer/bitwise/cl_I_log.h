// cl_I internals for logical operations

#ifndef _CL_I_LOG_H
#define _CL_I_LOG_H

#include "cln/number.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

// Liefert die Anzahl Digits, die ein Integer als DS bräuchte.
// (Leicht aufgerundet.)
inline uintC I_to_DS_need (const cl_I& x)
{
	if (fixnump(x))
		return FN_maxlength; // das wird reichen
	else
		return TheBignum(x)->length;
}

// Integer to Digit sequence, n Digits
// I_to_DS_n(obj,n,ptr=);
// Integer obj zu einer Digit sequence MSDptr/n/LSDptr machen,
// die genau n Digits hat (sollte n >= Bedarf und >= FN_maxlength sein).
// Die neue Digit-sequence darf modifiziert werden.
// < ptr: MSDptr der neuen DS
// Dabei wird num_stack erniedrigt.
  #define I_to_DS_n(obj,n,ptr_zuweisung)  \
    {var uintD* destptr;						\
     num_stack_alloc(n,,destptr=);					\
     ptr_zuweisung I_to_DS_n_aux(obj,n,destptr);			\
    }
  extern uintD* I_to_DS_n_aux (const cl_I& obj, uintC n, uintD* destptr);

}  // namespace cln

#endif /* _CL_I_LOG_H */
