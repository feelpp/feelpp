// L2_to_I() helper.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"

#if (cl_word_size < 64)

#include "base/digitseq/cl_DS.h"

namespace cln {

cl_private_thing cl_I_constructor_from_L2 (sint32 wert_hi, uint32 wert_lo)
{
	if (wert_hi == 0) {
		if ((wert_lo & minus_bit(cl_value_len-1)) == 0)
			return (cl_private_thing)(cl_combine(cl_FN_tag,wert_lo));
	}
	elif (wert_hi == ~(sint32)0) {
		if ((~wert_lo & minus_bit(cl_value_len-1)) == 0)
			return (cl_private_thing)(cl_combine(cl_FN_tag,(sint32)wert_lo));
	}
	// Bignum erzeugen:
	// (dessen Länge  bn_minlength <= n <= ceiling(64/intDsize)  erfüllt)
	#define FILL_1_DIGIT(l,i,from) \
		arrayLSref(ptr->data,l,i) = (uintD)from;
	#define FILL_2_DIGIT(l,i,from) \
		arrayLSref(ptr->data,l,i) = (uintD)from; \
		arrayLSref(ptr->data,l,i+1) = (uintD)(from>>intDsize);
	#define FILL_3_DIGIT(l,i,from) \
		arrayLSref(ptr->data,l,i) = (uintD)from; from>>=intDsize; \
		arrayLSref(ptr->data,l,i+1) = (uintD)from; \
		arrayLSref(ptr->data,l,i+2) = (uintD)(from>>intDsize);
	#define FILL_4_DIGIT(l,i,from) \
		arrayLSref(ptr->data,l,i) = (uintD)from; from>>=intDsize; \
		arrayLSref(ptr->data,l,i+1) = (uintD)from; from>>=intDsize; \
		arrayLSref(ptr->data,l,i+2) = (uintD)from; \
		arrayLSref(ptr->data,l,i+3) = (uintD)(from>>intDsize);
	#if (intDsize==64)
	#define FILL_1  FILL_1_DIGIT(1,0,highlow64(wert_hi,wert_lo));
	#endif
	#if (32/intDsize==1)
	#define FILL_1  FILL_1_DIGIT(1,0,wert_lo);
	#define FILL_2  FILL_1_DIGIT(2,1,wert_hi); FILL_1_DIGIT(2,0,wert_lo);
	#endif
	#if (32/intDsize==2)
	#define FILL_1  FILL_1_DIGIT(1,0,wert_lo);
	#define FILL_2  FILL_2_DIGIT(2,0,wert_lo);
	#define FILL_3  FILL_1_DIGIT(3,2,wert_hi); FILL_2_DIGIT(3,0,wert_lo);
	#define FILL_4  FILL_2_DIGIT(4,2,wert_hi); FILL_2_DIGIT(4,0,wert_lo);
	#endif
	#if (32/intDsize==4)
	#define FILL_1  FILL_1_DIGIT(1,0,wert_lo);
	#define FILL_2  FILL_2_DIGIT(2,0,wert_lo);
	#define FILL_3  FILL_3_DIGIT(3,0,wert_lo);
	#define FILL_4  FILL_4_DIGIT(4,0,wert_lo);
	#define FILL_5  FILL_1_DIGIT(5,4,wert_hi); FILL_4_DIGIT(5,0,wert_lo);
	#define FILL_6  FILL_2_DIGIT(6,4,wert_hi); FILL_4_DIGIT(6,0,wert_lo);
	#define FILL_7  FILL_3_DIGIT(7,4,wert_hi); FILL_4_DIGIT(7,0,wert_lo);
	#define FILL_8  FILL_4_DIGIT(8,4,wert_hi); FILL_4_DIGIT(8,0,wert_lo);
	#endif
	if (wert_hi >= 0) {
		#define IF_LENGTH(i)  \
		  if ((bn_minlength <= i) && (i*intDsize <= 64)		\
		    && (!((i+1)*intDsize <= 64)				\
		        || (i*intDsize-1 < 32				\
		            ? ((wert_hi == 0) && (wert_lo < (uint32)bitc(i*intDsize-1))) \
		            : ((uint32)wert_hi < (uint32)bitc(i*intDsize-1-32)) \
		     ) )   )
		#define ALLOC(i)  \
		  var cl_heap_bignum* ptr = allocate_bignum(i);
		#define OK  \
		  return (cl_private_thing)(ptr);
		IF_LENGTH(1)
			bignum1: { ALLOC(1); FILL_1; OK; }
		#if (intDsize <= 32)
		IF_LENGTH(2)
			bignum2: { ALLOC(2); FILL_2; OK; }
		#if (intDsize <= 16)
		IF_LENGTH(3)
			bignum3: { ALLOC(3); FILL_3; OK; }
		IF_LENGTH(4)
			bignum4: { ALLOC(4); FILL_4; OK; }
		#if (intDsize <= 8)
		IF_LENGTH(5)
			bignum5: { ALLOC(5); FILL_5; OK; }
		IF_LENGTH(6)
			bignum6: { ALLOC(6); FILL_6; OK; }
		IF_LENGTH(7)
			bignum7: { ALLOC(7); FILL_7; OK; }
		IF_LENGTH(8)
			bignum8: { ALLOC(8); FILL_8; OK; }
		#endif
		#endif
		#endif
		#undef IF_LENGTH
	} else {
		#define IF_LENGTH(i)  \
		  if ((bn_minlength <= i) && (i*intDsize <= 64)		\
		    && (!((i+1)*intDsize <= 64)				\
		        || (i*intDsize-1 < 32				\
		            ? ((wert_hi == ~(sint32)0) && (wert_lo >= (uint32)(-bitc(i*intDsize-1)))) \
		            : ((uint32)wert_hi >= (uint32)(-bitc(i*intDsize-1-32))) \
		     ) )   )
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
