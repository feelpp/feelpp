// I/O of floats.

#ifndef _CL_FLOAT_IO_H
#define _CL_FLOAT_IO_H

#include "cln/number_io.h"
#include "cln/float.h"

namespace cln {

// Undocumented input functions

// Wandelt eine Zeichenkette mit Float-Syntax in ein Float um.
// read_float(base,sign,string,index1,index4,index2,index3)
// > base: Lesebasis (=10)
// > sign: Vorzeichen (/=0 falls negativ)
// > string: Simple-String (enthält Ziffern und evtl. Punkt und Exponentmarker)
// > index1: Index vom Mantissenanfang (excl. Vorzeichen)
// > index4: Index nach dem Mantissenende
// > index2: Index beim Ende der Characters
// > index3: Index nach dem Dezimalpunkt (=index4 falls keiner da)
//   (also Mantisse mit index4-index1 Characters: Ziffern und max. 1 '.')
//   (also index4-index3 Nachkommaziffern)
//   (also bei index4<index2: index4 = Index des Exponent-Markers,
//    index4+1 = Index des Exponenten-Vorzeichens oder der ersten
//    Exponenten-Ziffer)
// < ergebnis: Float
extern const cl_F read_float (unsigned int base, float_format_t prec,
                  cl_signean sign, const char * string, uintC index1, uintC index4, uintC index2, uintC index3);

// The following does strictly the same as the general read_complex.
// It is here only so that you don't need the complex and rational number
// readers in order to read a float number. ("Treeshaking")
extern const cl_F read_float (const cl_read_flags& flags, const char * string, const char * string_limit, const char * * end_of_parse);
extern const cl_F read_float (std::istream& stream, const cl_read_flags& flags);

// Documented input functions

inline std::istream& operator>> (std::istream& stream, cl_F& result)
{
	extern cl_read_flags cl_F_read_flags;
	result = read_float(stream,cl_F_read_flags);
	return stream;
}


// Undocumented output functions


// Documented output functions

// Gibt ein Float aus.
// print_float(stream,z);
// > z: Float
// > stream: Stream
extern void print_float (std::ostream& stream, const cl_print_flags& flags, const cl_F& z);
extern void print_float (std::ostream& stream, const cl_print_number_flags& flags, const cl_F& z);
extern void print_float (std::ostream& stream, const cl_print_real_flags& flags, const cl_F& z);
extern void print_float (std::ostream& stream, const cl_print_float_flags& flags, const cl_F& z);

// Gibt ein Float binär (sehr primitiv) aus.
// print_float_binary(stream,z);
// > z: Float
// > stream: Stream
extern void print_float_binary (std::ostream& stream, const cl_F& z);

// The following does strictly the same as the general `fprint' for numbers.
// It is here only so that you don't need the complex printer
// in order to print a float. ("Treeshaking")

inline void fprint (std::ostream& stream, const cl_F& x)
{
	extern cl_print_flags default_print_flags;
	print_float(stream,default_print_flags,x);
}

CL_DEFINE_PRINT_OPERATOR(cl_F)

}  // namespace cln

#endif /* _CL_FLOAT_IO_H */
