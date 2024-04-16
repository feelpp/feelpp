// Q_to_I() helper.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"

#ifdef intQsize

#include "base/digitseq/cl_DS.h"

namespace cln {

cl_private_thing cl_I_constructor_from_Q (sint64 wert)
{
	var uint64 test = wert & (sint64)minus_bit(cl_value_len-1);
	// test enthält die Bits, die nicht in den Fixnum-Wert >= 0 reinpassen.
	if ((test == 0) || (test == (uint64)(sint64)minus_bit(cl_value_len-1)))
		return (cl_private_thing)(cl_combine(cl_FN_tag,wert));
	// Bignum erzeugen:
	// (dessen Länge  bn_minlength <= n <= ceiling(32/intDsize)  erfüllt)
	if (wert >= 0) {
		#define IF_LENGTH(i)  \
		  if ((bn_minlength <= i) && (i*intDsize <= 64)		\
		    && (!((i+1)*intDsize <= 64)				\
		        || ((uint64)wert < ((uint64)1 << (i*intDsize-1))) \
		     ) )
		IF_LENGTH(1)
			bignum1:
			{ var cl_heap_bignum* ptr = allocate_bignum(1);
			  arrayLSref(ptr->data,1,0) = (uintD)wert;
			  return (cl_private_thing)(ptr);
			}
		#if (intDsize <= 32)
		IF_LENGTH(2)
			bignum2:
			{ var cl_heap_bignum* ptr = allocate_bignum(2);
			  arrayLSref(ptr->data,2,0) = (uintD)wert;
			  arrayLSref(ptr->data,2,1) = (uintD)(wert>>intDsize);
			  return (cl_private_thing)(ptr);
			}
		#if (intDsize <= 16)
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
		#if (intDsize <= 8)
		IF_LENGTH(5)
			bignum5:
			{ var cl_heap_bignum* ptr = allocate_bignum(5);
			  arrayLSref(ptr->data,5,0) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,5,1) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,5,2) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,5,3) = (uintD)wert;
			  arrayLSref(ptr->data,5,4) = (uintD)(wert>>intDsize);
			  return (cl_private_thing)(ptr);
			}
		IF_LENGTH(6)
			bignum6:
			{ var cl_heap_bignum* ptr = allocate_bignum(6);
			  arrayLSref(ptr->data,6,0) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,6,1) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,6,2) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,6,3) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,6,4) = (uintD)wert;
			  arrayLSref(ptr->data,6,5) = (uintD)(wert>>intDsize);
			  return (cl_private_thing)(ptr);
			}
		IF_LENGTH(7)
			bignum7:
			{ var cl_heap_bignum* ptr = allocate_bignum(7);
			  arrayLSref(ptr->data,7,0) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,7,1) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,7,2) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,7,3) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,7,4) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,7,5) = (uintD)wert;
			  arrayLSref(ptr->data,7,6) = (uintD)(wert>>intDsize);
			  return (cl_private_thing)(ptr);
			}
		IF_LENGTH(8)
			bignum8:
			{ var cl_heap_bignum* ptr = allocate_bignum(8);
			  arrayLSref(ptr->data,8,0) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,8,1) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,8,2) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,8,3) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,8,4) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,8,5) = (uintD)wert; wert >>= intDsize;
			  arrayLSref(ptr->data,8,6) = (uintD)wert;
			  arrayLSref(ptr->data,8,7) = (uintD)(wert>>intDsize);
			  return (cl_private_thing)(ptr);
			}
		#endif
		#endif
		#endif
		#undef IF_LENGTH
	} else {
		#define IF_LENGTH(i)  \
		  if ((bn_minlength <= i) && (i*intDsize <= 64)		\
		    && (!((i+1)*intDsize <= 64)				\
		        || ((uint64)wert >= ((uint64)(-1) << (i*intDsize-1))) \
		     ) )
		IF_LENGTH(1)
			goto bignum1;
		#if (intDsize <= 32)
		IF_LENGTH(2)
			goto bignum2;
		#if (intDsize <= 16)
		IF_LENGTH(3)
			goto bignum3;
		IF_LENGTH(4)
			goto bignum4;
		#if (intDsize <= 8)
		IF_LENGTH(5)
			goto bignum5;
		IF_LENGTH(6)
			goto bignum6;
		IF_LENGTH(7)
			goto bignum7;
		IF_LENGTH(8)
			goto bignum8;
		#endif
		#endif
		#endif
		#undef IF_LENGTH
	}
}

}  // namespace cln

#endif
