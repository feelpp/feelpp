// Formatted output functions à la Common Lisp.

#ifndef _CL_FORMAT_H
#define _CL_FORMAT_H

#include "cln/number.h"
#include "cln/io.h"
#include "cln/float.h"

namespace cln {

// gibt arg als römische Zahl auf stream aus, z.B. 4 als IIII.
extern void format_old_roman (std::ostream& stream, const cl_I& arg);

// gibt arg als römische Zahl auf stream aus, z.B. 4 als IV.
extern void format_new_roman (std::ostream& stream, const cl_I& arg);

extern const char * const cl_format_tens [10];

// gibt die ganze Zahl arg im Klartext auf englisch auf den Stream aus.
extern void format_cardinal (std::ostream& stream, const cl_I& arg);

// gibt eine ganze Zahl arg als Abzählnummer im Klartext auf englisch
// auf den stream aus.
extern void format_ordinal (std::ostream& stream, const cl_I& arg);

// gibt count (>=0) Zeichen ch auf stream aus.
inline void format_padding (std::ostream& stream, sintL count, char ch)
{
	for (; count >= 0; count--)
		fprintchar(stream,ch);
}

// gibt auf den Stream stream aus:
// den String str, eventuell aufgefüllt mit Padding characters padchar.
// Und zwar so, daß die Breite mindestens mincol ist. Um das zu erreichen,
// werden mindestens minpad Zeichen eingefügt, eventuelle weitere dann in
// Blöcken à colinc Zeichen. Falls padleftflag, werden sie links eingefügt,
// sonst rechts vom String.
extern void format_padded_string (std::ostream& stream, sintL mincol, sintL colinc, sintL minpad, char padchar, bool padleftflag, const char * str);

// gibt den Integer arg auf den Stream aus:
// in Zahlenbasis base, mit Vorzeichen (+ nur falls >0 und positive-sign-flag),
// bei commaflag alle drei Stellen unterbrochen durch ein Zeichen commachar.
// Das Ganze links aufgefüllt mit padchar's, so daß die Gesamtbreite mindestens
// mincol ist.
extern void format_integer (std::ostream& stream, const cl_I& arg, unsigned int base, sintL mincol, char padchar, char commachar, uintL commainterval, bool commaflag, bool positive_sign_flag);

// format_scale_exponent(arg) liefert zur Floating-Point-Zahl arg
// drei Werte: mantissa und n, mit
// ganzem n und mantissa floating-point, 0.1 <= mantissa < 1,
// arg = mantissa * 10^n * sign (also 10^(n-1) <= abs(arg) < 10^n ).
// (Bei arg=0.0: 0.0 und n=0.)
extern const decoded_float format_scale_exponent (const cl_F& arg);

// format_float_to_string(arg,width,d,k,dmin)
// ergibt einen String zum Floating-point arg:
// er hat den Wert von abs(arg)*expt(10,k), dabei mind. d Nachkommastellen
// und höchstens die Länge width (width<=0 -> keine Einschränkung).
// Trotzdem wird nicht auf weniger als dmin Stellen gerundet.
struct digits_with_dot {
	char * string; // Mit malloc_hook() alloziert, mit free_hook() freizugeben.
	uintL length; // strlen(string)
	bool dot_comes_first; // string[0] == '.' ?
	bool dot_comes_last; // string[strlen(string)-1] == '.' ?
	uintL dot_position; // string[dot_position] is '.'
// Constructor.
	digits_with_dot (char* s, uintL l, bool df, bool dl, uintL dp)
		: string(s), length(l), dot_comes_first(df), dot_comes_last(dl), dot_position(dp) {}
};
extern const digits_with_dot format_float_to_string (const cl_F& arg, const sintL width, const sintL d, const sintL k, const sintL dmin);

}  // namespace cln

#endif /* _CL_FORMAT_H */
