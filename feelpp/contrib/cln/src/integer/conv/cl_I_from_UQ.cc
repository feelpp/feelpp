// UQ_to_I() helper.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"

#ifdef intQsize

#include "base/digitseq/cl_DS.h"

namespace cln {

cl_private_thing cl_I_constructor_from_UQ (uint64 wert)
{
	if ((wert & (sint64)minus_bit(cl_value_len-1)) == 0)
	   // Bits, die nicht in den Fixnum-Wert >= 0 reinpassen.
		return (cl_private_thing)(cl_combine(cl_FN_tag,wert));
	// Bignum erzeugen:
	// (dessen Länge  bn_minlength <= n <= ceiling((32+1)/intDsize)  erfüllt)
	#define UQ_maxlength  ceiling(64+1,intDsize)
	#define IF_LENGTH(i)  \
	  if ((bn_minlength <= i) && (i <= UQ_maxlength)	\
	    && (!(i+1 <= UQ_maxlength)				\
	        || ((uint64)wert < ((uint64)1 << (i*intDsize-1 < 64 ? i*intDsize-1 : 0))) \
	     ) )
	IF_LENGTH(1)
		{ var cl_heap_bignum* ptr = allocate_bignum(1);
		  arrayLSref(ptr->data,1,0) = (uintD)wert;
		  return (cl_private_thing)(ptr);
		}
	IF_LENGTH(2)
		{ var cl_heap_bignum* ptr = allocate_bignum(2);
		  arrayLSref(ptr->data,2,0) = (uintD)wert;
		  #if (intDsize>=64)
		  arrayLSref(ptr->data,2,1) = 0;
		  #else
		  arrayLSref(ptr->data,2,1) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	#if (intDsize <= 32)
	IF_LENGTH(3)
		{ var cl_heap_bignum* ptr = allocate_bignum(3);
		  arrayLSref(ptr->data,3,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,3,1) = (uintD)wert;
		  #if (2*intDsize>=64)
		  arrayLSref(ptr->data,3,2) = 0;
		  #else
		  arrayLSref(ptr->data,3,2) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	#if (intDsize <= 16)
	IF_LENGTH(4)
		{ var cl_heap_bignum* ptr = allocate_bignum(4);
		  arrayLSref(ptr->data,4,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,4,1) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,4,2) = (uintD)wert;
		  #if (3*intDsize>=64)
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
		  #if (4*intDsize>=64)
		  arrayLSref(ptr->data,5,4) = 0;
		  #else
		  arrayLSref(ptr->data,5,4) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	#if (intDsize <= 8)
	IF_LENGTH(6)
		{ var cl_heap_bignum* ptr = allocate_bignum(6);
		  arrayLSref(ptr->data,6,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,6,1) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,6,2) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,6,3) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,6,4) = (uintD)wert;
		  #if (5*intDsize>=64)
		  arrayLSref(ptr->data,6,5) = 0;
		  #else
		  arrayLSref(ptr->data,6,5) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	IF_LENGTH(7)
		{ var cl_heap_bignum* ptr = allocate_bignum(7);
		  arrayLSref(ptr->data,7,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,7,1) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,7,2) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,7,3) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,7,4) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,7,5) = (uintD)wert;
		  #if (6*intDsize>=64)
		  arrayLSref(ptr->data,7,6) = 0;
		  #else
		  arrayLSref(ptr->data,7,6) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	IF_LENGTH(8)
		{ var cl_heap_bignum* ptr = allocate_bignum(8);
		  arrayLSref(ptr->data,8,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,8,1) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,8,2) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,8,3) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,8,4) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,8,5) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,8,6) = (uintD)wert;
		  #if (7*intDsize>=64)
		  arrayLSref(ptr->data,8,7) = 0;
		  #else
		  arrayLSref(ptr->data,8,7) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	IF_LENGTH(9)
		{ var cl_heap_bignum* ptr = allocate_bignum(9);
		  arrayLSref(ptr->data,9,0) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,9,1) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,9,2) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,9,3) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,9,4) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,9,5) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,9,6) = (uintD)wert; wert >>= intDsize;
		  arrayLSref(ptr->data,9,7) = (uintD)wert;
		  #if (8*intDsize>=64)
		  arrayLSref(ptr->data,9,8) = 0;
		  #else
		  arrayLSref(ptr->data,9,8) = (uintD)(wert>>intDsize);
		  #endif
		  return (cl_private_thing)(ptr);
		}
	#endif
	#endif
	#endif
	#undef IF_LENGTH
	#undef UQ_maxlength
}

}  // namespace cln

#endif
