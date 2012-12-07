// I/O of rational numbers.

#ifndef _CL_RATIONAL_IO_H
#define _CL_RATIONAL_IO_H

#include "cln/number_io.h"
#include "cln/rational.h"

namespace cln {

// Undocumented input functions

// Wandelt eine Zeichenkette mit Rational-Syntax in eine rationale Zahl um.
// read_rational(base,sign,string,index1,index3,index2)
// > base: Lesebasis (>=2, <=36)
// > sign: Vorzeichen (/=0 falls negativ)
// > string: Simple-String (enthält Ziffern mit Wert <base und Bruchstrich)
// > index1: Index der ersten Ziffer
// > index3: Index von '/'
// > index2: Index nach der letzten Ziffer
//   (also index3-index1 Zähler-Ziffern, index2-index3-1 Nenner-Ziffern)
// < ergebnis: rationale Zahl
extern const cl_RA read_rational (unsigned int base,
                   cl_signean sign, const char * string, uintC index1, uintC index3, uintC index2);

// The following does strictly the same as the general read_complex.
// It is here only so that you don't need the complex and float number
// readers in order to read an rational number. ("Treeshaking")
extern const cl_RA read_rational (const cl_read_flags& flags, const char * string, const char * string_limit, const char * * end_of_parse);
extern const cl_RA read_rational (std::istream& stream, const cl_read_flags& flags);

// Documented input functions

inline std::istream& operator>> (std::istream& stream, cl_RA& result)
{
	extern cl_read_flags cl_RA_read_flags;
	result = read_rational(stream,cl_RA_read_flags);
	return stream;
}


// Undocumented output functions

// Gibt eine rationale Zahl aus.
// print_rational(stream,base,z);
// > z: rationale Zahl
// > base: Basis (>=2, <=36)
// > stream: Stream
extern void print_rational (std::ostream& stream, unsigned int base, const cl_RA& z);


// Documented output functions

// Gibt eine Zahl aus.
// print_rational(stream,flags,z);
// > z: Zahl
// > stream: Stream
// > flags: Ausgabe-Parameter
extern void print_rational (std::ostream& stream, const cl_print_flags& flags, const cl_RA& z);
extern void print_rational (std::ostream& stream, const cl_print_number_flags& flags, const cl_RA& z);
extern void print_rational (std::ostream& stream, const cl_print_real_flags& flags, const cl_RA& z);
extern void print_rational (std::ostream& stream, const cl_print_rational_flags& flags, const cl_RA& z);

// The following does strictly the same as the general `fprint' for numbers.
// It is here only so that you don't need the complex and long-float number
// printers in order to print an integer. ("Treeshaking")

inline void fprint (std::ostream& stream, const cl_RA& x)
{
	extern cl_print_flags default_print_flags;
	print_rational(stream,default_print_flags,x);
}

CL_DEFINE_PRINT_OPERATOR(cl_RA)

}  // namespace cln

#endif /* _CL_RATIONAL_IO_H */
