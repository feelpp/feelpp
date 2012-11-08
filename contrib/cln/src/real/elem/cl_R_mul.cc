// binary operator *

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"
#include "rational/cl_RA.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "float/cl_F.h"
#include "cln/sfloat.h"
#include "float/sfloat/cl_SF.h"
#include "cln/ffloat.h"
#include "float/ffloat/cl_FF.h"
#include "cln/dfloat.h"
#include "float/dfloat/cl_DF.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"
#include "base/cl_N.h"

namespace cln {

ALL_cl_LF_OPERATIONS_SAME_PRECISION()

const cl_R operator* (const cl_R& x, const cl_R& y)
{
	if (eq(x,0)) { return 0; } // 0 * y = exakte 0
	elif (eq(y,0)) { return 0; } // x * 0 = exakte 0
	else
#define mul(a,b) a*b
	realcase6(x
	, /* I */	
			realcase6(y
			, /* I */	return mul(x,y);
			, /* RT */	return mul(x,y);
			, /* SF */	return mul(cl_I_to_SF(x),y);
			, /* FF */	return mul(cl_I_to_FF(x),y);
			, /* DF */	return mul(cl_I_to_DF(x),y);
			, /* LF */	return cl_LF_I_mul(y,x);
			);
	, /* RT */	
			realcase6(y
			, /* I */	return mul(x,y);
			, /* RT */	return mul(x,y);
			, /* SF */	return mul(cl_RA_to_SF(x),y);
			, /* FF */	return mul(cl_RA_to_FF(x),y);
			, /* DF */	return mul(cl_RA_to_DF(x),y);
			, /* LF */	return cl_LF_RA_mul(y,x);
			);
	, /* SF */	
			realcase6(y
			, /* I */	return mul(x,cl_I_to_SF(y));
			, /* RT */	return mul(x,cl_RA_to_SF(y));
			, /* SF */	return mul(x,y);
			, /* FF */	return cl_FF_to_SF(mul(cl_SF_to_FF(x),y));
			, /* DF */	return cl_DF_to_SF(mul(cl_SF_to_DF(x),y));
			, /* LF */	return cl_LF_to_SF(mul(cl_SF_to_LF(x,LFlen0(y)),y));
			);
	, /* FF */	
			realcase6(y
			, /* I */	return mul(x,cl_I_to_FF(y));
			, /* RT */	return mul(x,cl_RA_to_FF(y));
			, /* SF */	return cl_FF_to_SF(mul(x,cl_SF_to_FF(y)));
			, /* FF */	return mul(x,y);
			, /* DF */	return cl_DF_to_FF(mul(cl_FF_to_DF(x),y));
			, /* LF */	return cl_LF_to_FF(mul(cl_FF_to_LF(x,LFlen0(y)),y));
			);
	, /* DF */	
			realcase6(y
			, /* I */	return mul(x,cl_I_to_DF(y));
			, /* RT */	return mul(x,cl_RA_to_DF(y));
			, /* SF */	return cl_DF_to_SF(mul(x,cl_SF_to_DF(y)));
			, /* FF */	return cl_DF_to_FF(mul(x,cl_FF_to_DF(y)));
			, /* DF */	return mul(x,y);
			, /* LF */	return cl_LF_to_DF(mul(cl_DF_to_LF(x,LFlen0(y)),y));
			);
	, /* LF */	
			realcase6(y
			, /* I */	return cl_LF_I_mul(x,y);
			, /* RT */	return cl_LF_RA_mul(x,y);
			, /* SF */	return cl_LF_to_SF(mul(x,cl_SF_to_LF(y,LFlen0(x))));
			, /* FF */	return cl_LF_to_FF(mul(x,cl_FF_to_LF(y,LFlen0(x))));
			, /* DF */	return cl_LF_to_DF(mul(x,cl_DF_to_LF(y,LFlen0(x))));
			, /* LF */	return mul(x,y);
			);
	);
}

}  // namespace cln
