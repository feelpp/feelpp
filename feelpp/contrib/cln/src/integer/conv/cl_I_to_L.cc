// cl_I_to_L().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "cln/number.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "cln/io.h"
#include "cln/integer_io.h"
#include "cln/exception.h"
#include <sstream>

namespace cln {

sint32 cl_I_to_L (const cl_I& obj)
{
	if (fixnump(obj)) {
		// Fixnum
		var sintV wert = FN_to_V(obj);
		#if (intVsize>32)
		if ((sintV)(sint32)wert != wert)
			goto bad;
		#endif
		return (sint32)wert;
	} else { // Bignum
		var cl_heap_bignum* bn = TheBignum(obj);
		var uintC len = bn->length;
		if ((sintD)mspref(arrayMSDptr(bn->data,len),0) >= 0) {
			// Bignum > 0
			#define IF_LENGTH(i)  \
			  if ((bn_minlength <= i) /* genau i Digits überhaupt möglich? */	\
			    && (len == i) /* genau i Digits? */					\
			    /* 2^((i-1)*intDsize-1) <= obj < 2^(i*intDsize-1) */		\
			    && ( (i*intDsize > 32)						\
			         && ( ((i-1)*intDsize >= 32)					\
			              || (mspref(arrayMSDptr(bn->data,len),0) >= (uintD)bitc(31-(i-1)*intDsize)) \
			     ) )    )								\
			    goto bad;								\
			    else
			IF_LENGTH(1)
				return get_uint1D_Dptr(arrayLSDptr(bn->data,1));
			IF_LENGTH(2)
				return get_uint2D_Dptr(arrayLSDptr(bn->data,2));
			IF_LENGTH(3)
				return get_uint3D_Dptr(arrayLSDptr(bn->data,3));
			IF_LENGTH(4)
				return get_uint4D_Dptr(arrayLSDptr(bn->data,4));
			#undef IF_LENGTH
		} else {
			// Bignum < 0
			#define IF_LENGTH(i)  \
			  if ((bn_minlength <= i) /* genau i Digits überhaupt möglich? */	\
			    && (len == i) /* genau i Digits? */					\
			    /* - 2^(i*intDsize-1) <= obj < - 2^((i-1)*intDsize-1) */		\
			    && ( (i*intDsize > 32)						\
			         && ( ((i-1)*intDsize >= 32)					\
			              || (mspref(arrayMSDptr(bn->data,len),0) < (uintD)(-bitc(31-(i-1)*intDsize))) \
			     ) )    )								\
			    goto bad;								\
			    else
			IF_LENGTH(1)
				return get_sint1D_Dptr(arrayLSDptr(bn->data,1));
			IF_LENGTH(2)
				return get_sint2D_Dptr(arrayLSDptr(bn->data,2));
			IF_LENGTH(3)
				return get_sint3D_Dptr(arrayLSDptr(bn->data,3));
			IF_LENGTH(4)
				return get_sint4D_Dptr(arrayLSDptr(bn->data,4));
			#undef IF_LENGTH
		}
	}
	bad: // unpassendes Objekt
	std::ostringstream buf;
	fprint(buf, "Not a 32-bit integer: ");
	fprint(buf, obj);
	throw runtime_exception(buf.str());
}

}  // namespace cln
