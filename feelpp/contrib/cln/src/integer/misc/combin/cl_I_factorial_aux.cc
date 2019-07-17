// cl_I_prod_ungerade().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/misc/combin/cl_I_combin.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I cl_I_prod_ungerade (uintL a, uintL b)
{
	var uintL diff = b-a; // Anzahl der Faktoren
	if (diff <= 4) {
		// Produkt iterativ bilden
		var cl_I faktor = L_to_FN(2*b+1); // 2*b+1 als letzter Faktor
		var cl_I produkt = faktor;
		var uintC count;
		dotimesC(count,diff-1,
		  { faktor = faktor-2; // nächster Faktor
		    produkt = faktor*produkt; // mit bisherigem Produkt multiplizieren
		  });
		return produkt;
	} else {
		// Produkt rekursiv bilden
		var uintL c = floor(a+b,2); // c:=floor((a+b)/2)
		return cl_I_prod_ungerade(a,c) * cl_I_prod_ungerade(c,b); // zwei Teilprodukte
	}
}

}  // namespace cln
