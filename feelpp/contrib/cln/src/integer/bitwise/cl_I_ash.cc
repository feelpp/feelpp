// ash().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_I ash (const cl_I& x, sintC y)
{
    // Methode:
    // x = 0 -> 0 als Ergebnis
    // y = 0 -> x als Ergebnis
    // y > 0 -> y = intDsize*k + i, j=k+(1 falls i>0, 0 falls i=0).
    //          j Wörter mehr reservieren, k Nullwörter, dann übertragen,
    //          bei i>0: um i Bits links schieben (i=1 geht einfacher).
    // y < 0 -> y <= - intDsize * (Länge(A0) in Digits) -> Ergebnis = 0 oder -1.
    //          Sonst: -y = intDsize*k + i mit k<Länge(A0).
    //                  Übertrage die (Länge(A0)-k) MSDigits,
    //                  falls i>0: schiebe sie um i Bits nach rechts (i=1 geht einfacher).
	if (zerop(x))
		return 0;		// x=0 -> 0 als Ergebnis
	if (y == 0)
		return x;		// y=0 -> x als Ergebnis
	CL_ALLOCA_STACK;
	if (y >= 0) {
	        // y>0
		var uintC y_ = (uintC)y;
		var uintL i = y_%intDsize; // i = y mod intDsize, >=0, <intDsize
		var uintC k = floor(y_,intDsize); // k = y div intDsize, >=0, <2^intCsize
		var uintD* LSDptr;
		var uintC len;
		var const uintD* x_LSDptr;
		I_to_NDS_nocopy(x, ,len=,x_LSDptr=,false,); // DS zu x bilden.
		if (k >= (uintC)(~len)) // kann len+k+1 Überlauf geben?
			{ throw ash_exception(y); } // ja -> Fehler
		num_stack_alloc_1(len+k,,LSDptr=);
		LSDptr = clear_loop_lsp(LSDptr,k); // k Nulldigits
	       {var uintD* MSDptr = copy_loop_lsp(x_LSDptr,LSDptr,len);
		// Nun ist MSDptr/len/LSDptr die DS zu x.
		// Oberhalb von ihr liegen k Nulldigits, unterhalb ist 1 Digit Platz.
		// MSDptr/len+k/.. ist jetzt die Gesamt-DS.
		// Noch um i Bits nach links schieben:
		if (!(i==0)) // Bei i>0
		  { // noch ein weiteres Digit dazunehmen (Vorzeichen)
		    {var uintD sign = sign_of_sintD(mspref(MSDptr,0));
		     lsprefnext(MSDptr) = sign;
		     len++;
		    }
		    // Schiebeschleife: die unteren len Digits um i Bits schieben
		    if (i==1)
		      { shift1left_loop_lsp(LSDptr,len); }
		    else
		      { shiftleft_loop_lsp(LSDptr,len,i,0); }
		  }
		return DS_to_I(MSDptr,len+k);
	       }
	} else {
		// y<0
		var uintC y_ = (uintC)(-y); // Wert von -y, >0
		var uintL i = y_%intDsize; // i = (-y) mod intDsize, >=0, <intDsize
		var uintC k = floor(y_,intDsize); // k = (-y) div intDsize, >=0
		// DS zu x bilden:
		var uintD* MSDptr;
		var uintC len;
		I_to_NDS(x, MSDptr=,len=,); // DS zu x bilden.
		if (k>=len) goto sign; // -y >= intDsize*len -> Vorzeichen von x zurück
		len -= k; // rechte k Digits einfach streichen
		// Noch ist len>0. Um i Bits nach rechts schieben:
		if (!(i==0)) // Bei i>0:
		  { // Schiebe len Digits ab MSDptr um i Bits nach rechts:
		    if (i==1)
		      { shift1right_loop_msp(MSDptr,len,sign_of_sintD(mspref(MSDptr,0))); }
		      else
		      { shiftrightsigned_loop_msp(MSDptr,len,i); }
		  }
		return DS_to_I(MSDptr,len);
	}
sign:	// Ergebnis ist 0, falls x>=0, und -1, falls x<0:
	return (minusp(x) ? cl_I(-1) : cl_I(0));
}

}  // namespace cln
