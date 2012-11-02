// cl_F_ln2.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "base/digitseq/cl_DS.h"
#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "float/cl_F.h"

namespace cln {


cl_LF& cl_LF_ln2()
{
	// Mantisse von ln(2) :
	static const uintD ln2_mantisse [64/intDsize] =
		#include "cl_F_ln2_var.h"
	static cl_LF val = encode_LF_array(0,0,ln2_mantisse,64/intDsize);
	return val;
}

// Problem: If someone changes free_hook, the destructor of this
// will call the new hook, passing it some pointer obtained by the old
// malloc_hook. ??

const cl_SF& cl_SF_ln2()
{
	static const cl_SF val = cl_LF_to_SF(cl_LF_ln2());
	return val;
}

const cl_FF& cl_FF_ln2()
{
	static const cl_FF val = cl_LF_to_FF(cl_LF_ln2());
	return val;
}

const cl_DF& cl_DF_ln2()
{
	static const cl_DF val = cl_LF_to_DF(cl_LF_ln2());
	return val;
}

}  // namespace cln

