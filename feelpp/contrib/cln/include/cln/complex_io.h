// I/O of complex numbers.

#ifndef _CL_COMPLEX_IO_H
#define _CL_COMPLEX_IO_H

#include "cln/number_io.h"
#include "cln/complex.h"

namespace cln {

// Undocumented input functions

extern const cl_N read_complex (const cl_read_flags& flags, const char * string, const char * string_limit, const char * * end_of_parse);
extern const cl_N read_complex (std::istream& stream, const cl_read_flags& flags);

// Documented input functions

inline std::istream& operator>> (std::istream& stream, cl_N& result)
{
	extern cl_read_flags cl_N_read_flags;
	result = read_complex(stream,cl_N_read_flags);
	return stream;
}


// Undocumented output functions


// Documented output functions

// Gibt eine komplexe Zahl aus.
// print_complex(stream,z);
// > z: komplexe Zahl
// > stream: Stream
extern void print_complex (std::ostream& stream, const cl_print_flags& flags, const cl_N& z);
extern void print_complex (std::ostream& stream, const cl_print_number_flags& flags, const cl_N& z);

inline void fprint (std::ostream& stream, const cl_N& x)
{
	extern cl_print_flags default_print_flags;
	print_complex(stream,default_print_flags,x);
}

CL_DEFINE_PRINT_OPERATOR(cl_N)

}  // namespace cln

#endif /* _CL_COMPLEX_IO_H */
