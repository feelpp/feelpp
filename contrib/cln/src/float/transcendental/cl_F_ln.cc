// ln().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/transcendental/cl_F_tran.h"
#include "float/cl_F.h"
#include "float/sfloat/cl_SF.h"
#include "cln/integer.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F ln (const cl_F& x)
{
// Methode:
// d := (float-digits x),
// Genauigkeit um sqrt(d)+max(integer-length(e)) Bits erhöhen,
// (m,e) := (decode-float x), so daß 1/2 <= m < 1.
// m<2/3 -> m:=2m, e:=e-1, so daß 2/3 <= m <= 4/3.
// ln(m) errechnen, ln(x)=ln(m)+e*ln(2) als Ergebnis.

	// Rechengenauigkeit erhöhen und m,e,s bestimmen:
	if (longfloatp(x) && (TheLfloat(x)->len >= 110)) {
		DeclareType(cl_LF,x);
		var decoded_lfloat m_e_s = decode_float(extend(x,TheLfloat(x)->len+1));
		var cl_LF& m = m_e_s.mantissa;
		var cl_I& e = m_e_s.exponent;
		if (m < make_SF(0,0+SF_exp_mid,floor(bit(SF_mant_len+2),3))) { // Short-Float 2/3
			m = scale_float(m,1); // m verdoppeln
			e = minus1(e); // e decrementieren
		}
		var cl_F res = lnx_ratseries(m);
		if (!zerop(e))
			res = res + cl_float(e,m)*cl_ln2(m); // ln(m)+e*ln(2)
		return cl_float(res,x);
	} else {
		var decoded_float m_e_s = decode_float(cl_F_extendsqrtx(x));
		var cl_F& m = m_e_s.mantissa;
		var cl_I& e = m_e_s.exponent;
		if (m < make_SF(0,0+SF_exp_mid,floor(bit(SF_mant_len+2),3))) { // Short-Float 2/3
			m = scale_float(m,1); // m verdoppeln
			e = minus1(e); // e decrementieren
		}
		var cl_F res = lnx_naive(m);
		if (!zerop(e))
			res = res + cl_float(e,m)*cl_ln2(m); // ln(m)+e*ln(2)
		return cl_float(res,x);
	}
}

}  // namespace cln
