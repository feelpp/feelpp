// UL_to_I() helper.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"

#if (cl_value_len <= 32) || (long_bitsize==32)

#include "base/digitseq/cl_DS.h"

namespace cln {

cl_private_thing cl_I_constructor_from_UL (uint32 wert)
{
#if (cl_value_len <= 32)
	if ((wert & minus_bit(cl_value_len-1)) == 0)
	   // Bits, die nicht in den Fixnum-Wert >= 0 reinpassen.
		return (cl_private_thing)(cl_combine(cl_FN_tag,wert));
	// Bignum erzeugen:
	// (dessen Länge  bn_minlength <= n <= ceiling((32+1)/intDsize)  erfüllt)
	#define UL_maxlength  ceiling(32+1,intDsize)
	#define IF_LENGTH(i)  \
	  if ((bn_minlength <= i) && (i <= UL_maxlength)	\
	    && (!(i+1 <= UL_maxlength)				\
	        || ((uint32)wert < (uint32)bitc(i*intDsize-1))	\
	     ) )
	IF_LENGTH(1)
		{ var cl_heap_bignum* ptr = allocate_bignum(1);
		  arrayLSref(ptr->data,1,0) = wert;
		  return (cl_private_thing)(ptr);
		}
	#if (intDsize <= 32)
	IF_LENGTH(2)
		{ var cl_heap_bignum* ptr = allocate_bignum(2);
		  arrayLSref(ptr->data,2,0) = (uintD)wert;
		  #if (intDsize>=32)
		  arrayLSref(ptr->data,2,1) = 0;
		  #else
		  arrayLSref(ptr->data,2,1) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	#if (intDsize <= 16)
	IF_LENGTH(3)
		{ var cl_heap_bignum* ptr = allocate_bignum(3);
		  arrayLSref(ptr->data,3,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,3,1) = (uintD)wert;
		  #if (2*intDsize>=32)
		  arrayLSref(ptr->data,3,2) = 0;
		  #else
		  arrayLSref(ptr->data,3,2) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	#if (intDsize <= 8)
	IF_LENGTH(4)
		{ var cl_heap_bignum* ptr = allocate_bignum(4);
		  arrayLSref(ptr->data,4,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,4,1) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,4,2) = (uintD)wert;
		  #if (3*intDsize>=32)
		  arrayLSref(ptr->data,4,3) = 0;
		  #else
		  arrayLSref(ptr->data,4,3) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	IF_LENGTH(5)
		{ var cl_heap_bignum* ptr = allocate_bignum(5);
		  arrayLSref(ptr->data,5,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,5,1) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,5,2) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,5,3) = (uintD)wert;
		  #if (4*intDsize>=32)
		  arrayLSref(ptr->data,5,4) = 0;
		  #else
		  arrayLSref(ptr->data,5,4) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	#endif
	#endif
	#endif
	#undef IF_LENGTH
	#undef UL_maxlength
#else // cl_value_len > 32
	// All bits fit in a fixnum value >= 0.
	return (cl_private_thing)(cl_combine(cl_FN_tag,wert));
#endif
}

}  // namespace cln

#endif
