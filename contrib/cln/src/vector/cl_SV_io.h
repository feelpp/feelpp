// I/O of vectors.

#ifndef _CL_SV_IO_H
#define _CL_SV_IO_H

#include "cln/number_io.h"
#include "cln/SV.h"
#include "cln/SV_number.h"

namespace cln {

// Gibt einen Vektor aus.
// print_vector(stream,flags,fun,z);
// > stream: Stream
// > flags: Flags
// > fun: Ausgabefunktion für die einzelnen Elemente
// > vector: Vektor
extern void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* fun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_SV_number& vector);

}  // namespace cln

#endif /* _CL_SV_IO_H */
