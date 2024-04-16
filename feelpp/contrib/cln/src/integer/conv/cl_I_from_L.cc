// L_to_I() helper.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"

#if (cl_value_len < 32) || (long_bitsize==32)

#include "base/digitseq/cl_DS.h"

namespace cln {

cl_private_thing cl_I_constructor_from_L (sint32 wert)
{
#if (cl_value_len < 32)
	var uint32 test = wert & minus_bit(cl_value_len-1);
	// test enthält die Bits, die nicht in den Fixnum-Wert >= 0 reinpassen.
	if ((test == 0) || (test == (uint32)minus_bit(cl_value_len-1)))
		{ return (cl_private_thing)(cl_combine(cl_FN_tag,wert)); }
	#if (intDsize==64)
	  // trivially generate a Bignum of length one digit
	  var cl_heap_bignum* ptr = allocate_bignum(1);
	  arrayLSref(ptr->data,1,0) = wert;
	  return (cl_private_thing)(ptr);
	#else
	// Bignum erzeugen:
	// (dessen Länge  bn_minlength <= n <= ceiling(32/intDsize)  erfüllt)
	if (bn_minlength == ceiling(32,intDsize)) {
		#if (intDsize==8)
		goto bignum4;
		#endif
		#if (intDsize==16)
		goto bignum2;
		#endif
		#if (intDsize==32)
		goto bignum1;
		#endif
	}
	if (wert >= 0) {
		#define IF_LENGTH(i)  \
		  if ((bn_minlength <= i) && (i*intDsize <= 32)		\
		    && (!((i+1)*intDsize <= 32)				\
		        || ((uint32)wert < (uint32)bitc(i*intDsize-1))	\
		     ) )
		#if (intDsize <= 32)
		IF_LENGTH(1)
			bignum1:
			{ var cl_heap_bignum* ptr = allocate_bignum(1);
			  arrayLSref(ptr->data,1,0) = wert;
			  return (cl_private_thing)(ptr);
			}
		#if (intDsize <= 16)
		IF_LENGTH(2)
			bignum2:
			{ var cl_heap_bignum* ptr = allocate_bignum(2);
			  arrayLSref(ptr->data,2,0) = (uintD)wert;
			  arrayLSref(ptr->data,2,1) = (uintD)(wert>>intDsize);
			  return (cl_private_thing)(ptr);
			}
		#if (intDsize <= 8)
		IF_LENGTH(3)
			bignum3:
			{ var cl_heap_bignum* ptr = allocate_bignum(3);
			  arrayLSref(ptr->data,3,0) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,3,1) = (uintD)wert;
			  arrayLSref(ptr->data,3,2) = (uintD)(wert>>intDsize);
			  return (cl_private_thing)(ptr);
			}
		IF_LENGTH(4)
			bignum4:
			{ var cl_heap_bignum* ptr = allocate_bignum(4);
			  arrayLSref(ptr->data,4,0) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,4,1) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,4,2) = (uintD)wert;
			  arrayLSref(ptr->data,4,3) = (uintD)(wert>>intDsize);
			  return (cl_private_thing)(ptr);
			}
		#endif
		#endif
		#endif
		#undef IF_LENGTH
	} else {
		#define IF_LENGTH(i)  \
		  if ((bn_minlength <= i) && (i*intDsize <= 32)		\
		    && (!((i+1)*intDsize <= 32)				\
		        || ((uint32)wert >= (uint32)(-bitc(i*intDsize-1))) \
		     ) )
		#if (intDsize <= 32)
		IF_LENGTH(1)
			goto bignum1;
		#if (intDsize <= 16)
		IF_LENGTH(2)
			goto bignum2;
		#if (intDsize <= 8)
		IF_LENGTH(3)
			goto bignum3;
		IF_LENGTH(4)
			goto bignum4;
		#endif
		#endif
		#endif
		#undef IF_LENGTH
	}
	#endif
	NOTREACHED
#else // cl_value_len >= 32
	// All bits fit in a fixnum value.
	return (cl_private_thing)(cl_combine(cl_FN_tag,wert));
#endif
}

}  // namespace cln

#endif
