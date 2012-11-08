// print_integer().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer_io.h"


// Implementation.

#include "cln/io.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

void print_integer (std::ostream& stream, unsigned int base, const cl_I& z)
{
	var cl_I abs_z;
	if (minusp(z)) {
		// z<0 -> Vorzeichen ausgeben:
		fprintchar(stream,'-');
		abs_z = -z;
	} else
		abs_z = z;
	CL_ALLOCA_STACK;
	var uintC need = cl_digits_need(abs_z,base);
	var uintB* ziffern = cl_alloc_array(uintB,need); // Platz für die Ziffern
	var cl_digits erg; erg.LSBptr = &ziffern[need];
	I_to_digits(abs_z,(uintD)base,&erg); // Umwandlung in Ziffern
	// Ziffern ausgeben:
	{
		var uintB* ptr = erg.MSBptr;
		var uintC count = erg.len;
		do { fprintchar(stream,*ptr++); } until (--count==0);
	}
}

}  // namespace cln
