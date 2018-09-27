// print_complex().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex_io.h"


// Implementation.

#include "cln/output.h"
#include "cln/complex.h"
#include "complex/cl_C.h"
#include "cln/real_io.h"

namespace cln {

void print_complex (std::ostream& stream, const cl_print_number_flags& flags, const cl_N& z)
{
	if (realp(z)) {
		DeclareType(cl_R,z);
		print_real(stream,flags,z);
	} else {
		DeclareType(cl_C,z);
		var cl_R re = realpart(z);
		var cl_R im = imagpart(z);
		if (flags.complex_readably) {
			// Common Lisp #C(re im) syntax
			fprintchar(stream,'#');
			fprintchar(stream,'C');
			fprintchar(stream,'(');
			print_real(stream,flags,re);
			fprintchar(stream,' ');
			print_real(stream,flags,im);
			fprintchar(stream,')');
		} else {
			// Standard mathematical notation: re + im i
			if (!eq(im,0)) {
				if (!eq(re,0)) {
					// Example: 3-7i
					print_real(stream,flags,re);
					if (minusp(im)) {
						fprintchar(stream,'-');
						print_real(stream,flags,-im);
					} else {
						fprintchar(stream,'+');
						print_real(stream,flags,im);
					}
					fprintchar(stream,'i');
				} else {
					// Example: 6i
					print_real(stream,flags,im);
					fprintchar(stream,'i');
				}
			} else {
				// Example: 8
				print_real(stream,flags,re);
			}
		}
	}
}

}  // namespace cln
