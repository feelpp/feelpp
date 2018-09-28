// atan().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

// Methode:
// Wert und Branch Cuts nach der Formel CLTL2, S. 307/312/313:
//   arctan(z) = (log(1+iz)-log(1-iz)) / 2i
// Sei z=x+iy, errechne u+iv = artanh(-y+ix) wie oben, Ergebnis v-iu.
// Real- und Imaginärteil des Ergebnisses sind Floats, außer wenn z reell oder
// rein imaginär ist.

inline const cl_C_R _atan (const cl_N& z)
{
	if (realp(z)) {
		DeclareType(cl_R,z);
		return atanh(0,z);
	} else {
		DeclareType(cl_C,z);
		return atanh(-imagpart(z),realpart(z));
	}
}

const cl_N atan (const cl_N& z)
{
	var cl_C_R u_v = _atan(z);
	var cl_R& u = u_v.realpart;
	var cl_R& v = u_v.imagpart;
	return complex(v,-u); // v-iu
}

}  // namespace cln
