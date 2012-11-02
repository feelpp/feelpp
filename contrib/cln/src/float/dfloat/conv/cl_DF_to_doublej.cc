// cl_DF_to_double().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/dfloat/cl_DF.h"


// Implementation.

namespace cln {

void cl_DF_to_double (const cl_DF& obj, dfloatjanus* val_)
{
	var dfloat val = TheDfloat(obj)->dfloat_value;
	// Der Exponent muß um DF_exp_mid-1022 erniedrigt werden.
	if (DF_exp_mid>1022)
	  #if (cl_word_size==64)
	  { var uintL exp = (val >> DF_mant_len) & (bit(DF_exp_len)-1); // e
	    if (exp < DF_exp_mid-1022+1)
	      { // produziere denormalisiertes Float
	        val = (val & minus_bit(DF_exp_len+DF_mant_len)) // selbes Vorzeichen
	              | ((sint64)0 << DF_mant_len) // Exponent 0
	              | (((val & (bit(DF_mant_len)-1)) | bit(DF_mant_len)) // Mantisse shiften
	                 >> (DF_exp_mid-1022+1 - exp) // shiften
	                );
	      }
	      else
	      { val -= (sint64)(DF_exp_mid - 1022) << DF_mant_len; }
	  }
	  #else
	  { var uintL exp = (val.semhi >> (DF_mant_len-32)) & (bit(DF_exp_len)-1); // e
	    if (exp < DF_exp_mid-1022+1)
	      { // produziere denormalisiertes Float
	        var uintL shiftcount = DF_exp_mid-1022+1 - exp;
	        val.mlo = val.mlo >> shiftcount; // Mantisse shiften
	        val.mlo |= val.semhi << (32-shiftcount);
	        val.semhi = (val.semhi & minus_bit(DF_exp_len+DF_mant_len-32)) // selbes Vorzeichen
	                    | ((sint32)0 << (DF_mant_len-32)) // Exponent 0
	                    | (((val.semhi & (bit(DF_mant_len-32)-1)) | bit(DF_mant_len-32)) // Mantisse shiften
	                       >> shiftcount // shiften
	                      );
	      }
	      else
	      { val.semhi -= (sint32)(DF_exp_mid - 1022) << (DF_mant_len-32); }
	  }
	  #endif
	val_->eksplicit = val;
}

}  // namespace cln
