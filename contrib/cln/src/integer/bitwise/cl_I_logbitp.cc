// logbitp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

bool logbitp (uintC x, const cl_I& y)
{
    // Methode:
    // Falls x>=intDsize*Länge(y), teste Vorzeichen von y.
    // Sonst x=intDsize*k+i, Teste Bit i vom Worte Nr. k+1 (von oben herab).

	var const uintD* yMSDptr;
	var uintC ylen;
	var const uintD* yLSDptr;
	I_to_NDS_nocopy(y, yMSDptr=,ylen=,yLSDptr=,true, { return false; } ); // DS zu y
	if (x < intDsize*ylen)
		// x ist >=0, < intDsize*ylen
		{ if (lspref(yLSDptr,floor(x,intDsize)) & bit(x%intDsize))
		    return true;
		    else
		    return false;
		}
	// Vorzeichen von y testen
	if (/* (ylen > 0) && */ ((sintD)mspref(yMSDptr,0) < 0))
	    return true;
	    else
	    return false;
}

}  // namespace cln
