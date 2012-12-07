// binary operator /

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

const cl_R operator/ (const cl_R& x, const cl_R& y)
{
	if (eq(x,0))
		// 0 / y = exakte 0, außer wenn y=0
		{ if (zerop(y))
			{ throw division_by_0_exception(); }
		  else
			return 0;
		}
	else
#define div(a,b) a/b
	realcase6(x
	, /* I */	
			realcase6(y
			, /* I */	return div(x,y);
			, /* RT */	return div(x,y);
			, /* SF */	return div(cl_I_to_SF(x),y);
			, /* FF */	return div(cl_I_to_FF(x),y);
			, /* DF */	return div(cl_I_to_DF(x),y);
			, /* LF */	return div(cl_I_to_LF(x,LFlen0(y)),y); // cf. cl_I_LF_div
			);
	, /* RT */	
			realcase6(y
			, /* I */	return div(x,y);
			, /* RT */	return div(x,y);
			, /* SF */	return div(cl_RA_to_SF(x),y);
			, /* FF */	return div(cl_RA_to_FF(x),y);
			, /* DF */	return div(cl_RA_to_DF(x),y);
			, /* LF */	return cl_RA_LF_div(x,y);
			);
	, /* SF */	
			realcase6(y
			, /* I */	return div(x,cl_I_to_SF(y));
			, /* RT */	return div(x,cl_RA_to_SF(y));
			, /* SF */	return div(x,y);
			, /* FF */	return cl_FF_to_SF(div(cl_SF_to_FF(x),y));
			, /* DF */	return cl_DF_to_SF(div(cl_SF_to_DF(x),y));
			, /* LF */	return cl_LF_to_SF(div(cl_SF_to_LF(x,LFlen0(y)),y));
			);
	, /* FF */	
			realcase6(y
			, /* I */	return div(x,cl_I_to_FF(y));
			, /* RT */	return div(x,cl_RA_to_FF(y));
			, /* SF */	return cl_FF_to_SF(div(x,cl_SF_to_FF(y)));
			, /* FF */	return div(x,y);
			, /* DF */	return cl_DF_to_FF(div(cl_FF_to_DF(x),y));
			, /* LF */	return cl_LF_to_FF(div(cl_FF_to_LF(x,LFlen0(y)),y));
			);
	, /* DF */	
			realcase6(y
			, /* I */	return div(x,cl_I_to_DF(y));
			, /* RT */	return div(x,cl_RA_to_DF(y));
			, /* SF */	return cl_DF_to_SF(div(x,cl_SF_to_DF(y)));
			, /* FF */	return cl_DF_to_FF(div(x,cl_FF_to_DF(y)));
			, /* DF */	return div(x,y);
			, /* LF */	return cl_LF_to_DF(div(cl_DF_to_LF(x,LFlen0(y)),y));
			);
	, /* LF */	
			realcase6(y
			, /* I */	return cl_LF_I_div(x,y);
			, /* RT */	return cl_LF_RA_div(x,y);
			, /* SF */	return cl_LF_to_SF(div(x,cl_SF_to_LF(y,LFlen0(x))));
			, /* FF */	return cl_LF_to_FF(div(x,cl_FF_to_LF(y,LFlen0(x))));
			, /* DF */	return cl_LF_to_DF(div(x,cl_DF_to_LF(y,LFlen0(x))));
			, /* LF */	return div(x,y);
			);
	);
}

}  // namespace cln
