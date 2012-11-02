// cl_I equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "base/cl_N.h"
#include "integer/cl_I.h"

namespace cln {

static inline uint32 equal_hashcode (const cl_FN& x)
{
	var cl_signean sign;
	var uintV x_ = FN_to_V(x); // x als intVsize-Bit-Zahl
	if (FN_V_minusp(x,(sintV)x_)) {
		x_ = -x_;
		sign = -1;
	} else {
		sign = 0;
		if (x_ == 0)
			return 0;
	}
	var uintL s;
        #if (intVsize > 32)
        integerlength64(x_, s = 64 - );
        var uint32 msd = (x_ << s) >> 32;
        var sintL exp = 64-s;
        #else
	integerlength32(x_, s = 32 - );
	var uint32 msd = x_ << s;
	var sintL exp = 32-s;
        #endif
	return equal_hashcode_low(msd,exp,sign);
}

static inline uint32 equal_hashcode (const cl_BN& x)
{
	var const uintD* MSDptr;
	var uintC len;
	BN_to_NDS_nocopy(x, MSDptr = , len = ,);
	// Nicht alle führenden intDsize+1 Bits sind gleich.
#if (intDsize==64)
	var uint64 msd = mspref(MSDptr,0);
	var uint64 msd2 = (len >= 2 ? mspref(MSDptr,1) : 0);
	var cl_signean sign;
	if ((sint64)msd < 0) { // falls <0, negieren
		sign = -1;
		// msd|msd2 := - msd|msd2 - (1 falls noch weitere Bits /= 0)
		msd = ~msd; msd2 = ~msd2;
		if ((len <= 2)
		    || !test_loop_msp(MSDptr mspop 2, len - 2)
		   ) {
			msd2++;
			if (msd2 == 0)
				msd++;
		}
	} else {
		sign = 0;
	}
	var sintC exp = len * intDsize;
	// Nicht alle führenden 65 Bits sind =0.
	if (msd==0) {
		msd = msd2;
		exp -= 64;
	} else {
		var uintL s;
		integerlength64(msd, s = 64 - );
		if (s > 0)
			msd = (msd << s) | (msd2 >> (64-s));
		exp -= s;
	}
	return equal_hashcode_low((uint32)(msd>>32),exp,sign);
#else // (intDsize<=32)
	var uint32 msd;
	var uint32 msd2;
	if (len >= 64/intDsize) {
		msd = get_32_Dptr(MSDptr);
		msd2 = get_32_Dptr(MSDptr mspop 32/intDsize);
	} elif (len > 32/intDsize) {
		msd = get_32_Dptr(MSDptr);
		msd2 = get_max32_Dptr(intDsize*len-32, MSDptr mspop 32/intDsize)
		       << (64-intDsize*len);
	} elif ((32/intDsize == 1) || (len == 32/intDsize)) {
		msd = get_32_Dptr(MSDptr);
		msd2 = 0;
	} else { // (len > 0) && (len < 32/intDsize)
		msd = get_max32_Dptr(intDsize*len,MSDptr) << (32-intDsize*len);
		msd2 = 0;
	}
	var cl_signean sign;
	if ((sint32)msd < 0) { // falls <0, negieren
		sign = -1;
		// msd|msd2 := - msd|msd2 - (1 falls noch weitere Bits /= 0)
		msd = ~msd; msd2 = ~msd2;
		if ((len <= 64/intDsize)
		    || !test_loop_msp(MSDptr mspop 64/intDsize, len - 64/intDsize)
		   ) {
			msd2++;
			if (msd2 == 0)
				msd++;
		}
	} else {
		sign = 0;
	}
	var sintC exp = len * intDsize;
	// Nicht alle führenden intDsize+1 Bits sind =0.
	// Wegen intDsize<=32: Nicht alle führenden 33 Bits sind =0.
	if (msd==0) {
		msd = msd2;
		exp -= 32;
	}
	// Nicht alle führenden 32 Bits sind =0.
	// Führendes Bit auf 1 normalisieren:
	else {
		var uintL s;
		integerlength32(msd, s = 32 - );
		if (s > 0)
			msd = (msd << s) | (msd2 >> (32-s));
		exp -= s;
	}
	return equal_hashcode_low(msd,exp,sign);
#endif
}

CL_INLINE uint32 CL_INLINE_DECL(equal_hashcode) (const cl_I& x)
{
	if (fixnump(x)) {
		DeclareType(cl_FN,x);
		return equal_hashcode(x);
	} else {
		DeclareType(cl_BN,x);
		return equal_hashcode(x);
	}
}

}  // namespace cln
