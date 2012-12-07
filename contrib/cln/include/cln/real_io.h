// I/O of real numbers.

#ifndef _CL_REAL_IO_H
#define _CL_REAL_IO_H

#include "cln/number_io.h"
#include "cln/real.h"

namespace cln {

// Undocumented input functions

// The following does strictly the same as the general read_complex.
// It is here only so that you don't need the complex number reader
// in order to read an rational number. ("Treeshaking")
extern const cl_R read_real (const cl_read_flags& flags, const char * string, const char * string_limit, const char * * end_of_parse);
extern const cl_R read_real (std::istream& stream, const cl_read_flags& flags);

// Documented input functions

inline std::istream& operator>> (std::istream& stream, cl_R& result)
{
	extern cl_read_flags cl_R_read_flags;
	result = read_real(stream,cl_R_read_flags);
	return stream;
}


// Undocumented output functions


// Documented output functions

// Gibt eine Zahl aus.
// print_real(stream,flags,z);
// > z: Zahl
// > stream: Stream
// > flags: Ausgabe-Parameter
extern void print_real (std::ostream& stream, const cl_print_flags& flags, const cl_R& z);
extern void print_real (std::ostream& stream, const cl_print_number_flags& flags, const cl_R& z);
extern void print_real (std::ostream& stream, const cl_print_real_flags& flags, const cl_R& z);

// The following does strictly the same as the general `fprint' for numbers.
// It is here only so that you don't need the complex number printer
// in order to print an integer. ("Treeshaking")

inline void fprint (std::ostream& stream, const cl_R& x)
{
	extern cl_print_flags default_print_flags;
	print_real(stream,default_print_flags,x);
}

CL_DEFINE_PRINT_OPERATOR(cl_R)

}  // namespace cln

#endif /* _CL_REAL_IO_H */
