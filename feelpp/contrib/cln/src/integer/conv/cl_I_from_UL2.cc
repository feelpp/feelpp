// UL2_to_I() helper.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"

#if (cl_word_size < 64)

#include "base/digitseq/cl_DS.h"

namespace cln {

cl_private_thing cl_I_constructor_from_UL2 (uint32 wert_hi, uint32 wert_lo)
{
	if ((wert_hi == 0)
	    && ((wert_lo & minus_bit(cl_value_len-1)) == 0)
	   )
		return (cl_private_thing)(cl_combine(cl_FN_tag,wert_lo));
	// Bignum erzeugen:
	// (dessen Länge  bn_minlength <= n <= ceiling((64+1)/intDsize)  erfüllt)
	#define UL2_maxlength  ceiling(64+1,intDsize)
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
	#define FILL_2  FILL_1_DIGIT(2,1,0); FILL_1_DIGIT(2,0,highlow64(wert_hi,wert_lo));
	#endif
	#if (32/intDsize==1)
	#define FILL_1  FILL_1_DIGIT(1,0,wert_lo);
	#define FILL_2  FILL_1_DIGIT(2,1,wert_hi); FILL_1_DIGIT(2,0,wert_lo);
	#define FILL_3	FILL_1_DIGIT(3,2,0); FILL_1_DIGIT(3,1,wert_hi); FILL_1_DIGIT(3,0,wert_lo);
	#endif
	#if (32/intDsize==2)
	#define FILL_1  FILL_1_DIGIT(1,0,wert_lo);
	#define FILL_2  FILL_2_DIGIT(2,0,wert_lo);
	#define FILL_3  FILL_1_DIGIT(3,2,wert_hi); FILL_2_DIGIT(3,0,wert_lo);
	#define FILL_4  FILL_2_DIGIT(4,2,wert_hi); FILL_2_DIGIT(4,0,wert_lo);
	#define FILL_5	FILL_1_DIGIT(5,4,0); FILL_2_DIGIT(5,2,wert_hi); FILL_2_DIGIT(5,0,wert_lo);
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
	#define FILL_9	FILL_1_DIGIT(9,8,0); FILL_4_DIGIT(9,4,wert_hi); FILL_4_DIGIT(9,0,wert_lo);
	#endif
	#define IF_LENGTH(i)  \
	  if ((bn_minlength <= i) && (i <= UL2_maxlength)	\
	    && (!(i+1 <= UL2_maxlength)				\
	        || (i*intDsize-1 < 32				\
	            ? ((wert_hi == 0) && (wert_lo < (uint32)bitc(i*intDsize-1))) \
	            : (wert_hi < (uint32)bitc(i*intDsize-1-32)) \
	     ) )   )
	#define ALLOC(i)  \
	  var cl_heap_bignum* ptr = allocate_bignum(i);
	#define OK  \
	  return (cl_private_thing)(ptr);
	IF_LENGTH(1)
		{ ALLOC(1); FILL_1; OK; }
	IF_LENGTH(2)
		{ ALLOC(2); FILL_2; OK; }
	#if (intDsize <= 32)
	IF_LENGTH(3)
		{ ALLOC(3); FILL_3; OK; }
	#if (intDsize <= 16)
	IF_LENGTH(4)
		{ ALLOC(4); FILL_4; OK; }
	IF_LENGTH(5)
		{ ALLOC(5); FILL_5; OK; }
	#if (intDsize <= 8)
	IF_LENGTH(6)
		{ ALLOC(6); FILL_6; OK; }
	IF_LENGTH(7)
		{ ALLOC(7); FILL_7; OK; }
	IF_LENGTH(8)
		{ ALLOC(8); FILL_8; OK; }
	IF_LENGTH(9)
		{ ALLOC(9); FILL_9; OK; }
	#endif
	#endif
	#endif
	#undef IF_LENGTH
	#undef UL2_maxlength
}

}  // namespace cln

#endif
