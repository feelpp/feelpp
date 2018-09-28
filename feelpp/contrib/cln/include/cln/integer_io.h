// I/O of integers.

#ifndef _CL_INTEGER_IO_H
#define _CL_INTEGER_IO_H

#include "cln/number_io.h"
#include "cln/integer_class.h"

namespace cln {

// Undocumented input functions

// Wandelt eine Zeichenkette mit Integer-Syntax in ein Integer um.
// Punkte werden überlesen.
// read_integer(base,sign,string,index1,index2)
// > base: Lesebasis (>=2, <=36)
// > sign: Vorzeichen (/=0 falls negativ)
// > string: Simple-String (enthält Ziffern mit Wert <base und evtl. Punkt)
// > index1: Index der ersten Ziffer
// > index2: Index nach der letzten Ziffer
//   (also index2-index1 Ziffern, incl. evtl. Dezimalpunkt am Schluß)
// < ergebnis: Integer
extern const cl_I read_integer (unsigned int base,
                  cl_signean sign, const char * string, uintC index1, uintC index2);

// The following does strictly the same as the general read_complex.
// It is here only so that you don't need the rational, complex and float number
// readers in order to read an integer. ("Treeshaking")
extern const cl_I read_integer (const cl_read_flags& flags, const char * string, const char * string_limit, const char * * end_of_parse);
extern const cl_I read_integer (std::istream& stream, const cl_read_flags& flags);

// Documented input functions

inline std::istream& operator>> (std::istream& stream, cl_I& result)
{
	extern cl_read_flags cl_I_read_flags;
	result = read_integer(stream,cl_I_read_flags);
	return stream;
}


// Undocumented output functions

// Liefert zu einem Integer >=0  (write-to-string integer :base 10 :radix nil),
// also die Ziffernfolge als String.
// Mit malloc_hook() alloziert, mit free_hook() freizugeben.
extern char * cl_decimal_string (const cl_I& x);

// Gibt ein Integer aus.
// print_integer(stream,base,z);
// > z: Integer
// > base: Basis (>=2, <=36)
// > stream: Stream
extern void print_integer (std::ostream& stream, unsigned int base, const cl_I& z);
// Dasselbe als String. Mit malloc_hook() alloziert, mit free_hook() freizugeben.
extern char * print_integer_to_string (unsigned int base, const cl_I& z);


// Documented output functions

inline void fprintdecimal (std::ostream& stream, const cl_I& x)
{
	print_integer(stream,10,x);
}

inline void fprintbinary (std::ostream& stream, const cl_I& x)
{
	print_integer(stream,2,x);
}

inline void fprintoctal (std::ostream& stream, const cl_I& x)
{
	print_integer(stream,8,x);
}

inline void fprinthexadecimal (std::ostream& stream, const cl_I& x)
{
	print_integer(stream,16,x);
}

// Gibt eine Zahl aus.
// print_integer(stream,flags,z);
// > z: Zahl
// > stream: Stream
// > flags: Ausgabe-Parameter
extern void print_integer (std::ostream& stream, const cl_print_flags& flags, const cl_I& z);
extern void print_integer (std::ostream& stream, const cl_print_number_flags& flags, const cl_I& z);
extern void print_integer (std::ostream& stream, const cl_print_real_flags& flags, const cl_I& z);
extern void print_integer (std::ostream& stream, const cl_print_rational_flags& flags, const cl_I& z);

// The following does strictly the same as the general `fprint' for numbers.
// It is here only so that you don't need the rational number printer
// in order to print an integer. ("Treeshaking")

inline void fprint (std::ostream& stream, const cl_I& x)
{
	extern cl_print_flags default_print_flags;
	print_integer(stream,default_print_flags,x);
}

CL_DEFINE_PRINT_OPERATOR(cl_I)

}  // namespace cln

#endif /* _CL_INTEGER_IO_H */
