// I/O of ffloats.

#ifndef _CL_FFLOAT_IO_H
#define _CL_FFLOAT_IO_H

#include "cln/number_io.h"
#include "cln/ffloat.h"

namespace cln {

inline std::istream& operator>> (std::istream& stream, cl_FF& result)
{
	extern cl_read_flags cl_FF_read_flags;
	extern const cl_F read_float (std::istream&, const cl_read_flags&);
	result = As(cl_FF)(read_float(stream,cl_FF_read_flags));
	return stream;
}

// The following does strictly the same as the general `fprint' for floats.
// It is here only so that people don't need to include <cln/float_io.h>.
inline void fprint (std::ostream& stream, const cl_FF& x)
{
	extern void print_float (std::ostream& stream, const cl_print_flags& flags, const cl_F& z);
	extern cl_print_flags default_print_flags;
	print_float(stream,default_print_flags,x);
}
CL_DEFINE_PRINT_OPERATOR(cl_FF)

}  // namespace cln

#endif /* _CL_FFLOAT_IO_H */
