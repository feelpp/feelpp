// format_float_to_string().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "real/format-output/cl_format.h"


// Implementation.

// BUGS:
// - This is slow.

#include "cln/output.h"
#include "cln/malloc.h"
#include "cln/float.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "base/string/cl_spushstring.h"

namespace cln {

// format_float_to_string(arg,width,d,k,dmin)
// ergibt einen String zum Floating-point arg:
// er hat den Wert von abs(arg)*expt(10,k), dabei mind. d Nachkommastellen
// und höchstens die Länge width (width<=0 -> keine Einschränkung).
// Trotzdem wird nicht auf weniger als dmin Stellen gerundet.

const digits_with_dot format_float_to_string (const cl_F& arg, const sintL width, const sintL d, const sintL k, const sintL dmin)
{
	// One pre-allocated buffer. This reduces the allocation/free cost.
	static cl_spushstring digitstring;

	if (zerop(arg)) {
		var sintL places = (d < dmin ? dmin : d);
		if (width > 0)
			// width angegeben -> places := min(places,width-1)
			if (places >= width)
				places = width-1;
		// ein Punkt und places Nullen
		var char* string = (char *) malloc_hook(1+places+1);
		string[0] = '.';
		for (sintL i = 1; i <= places; i++) string[i] = '0';
		string[1+places] = '\0';
		return digits_with_dot(string, 1+places,
				true, (places==0), 0
			);
	}
	// significand : Integer >0
	// expon : Integer
	// mantprec : Anzahl der echten Mantissenbits von significand
	// (also 2^mantprec <= significand < 2^(mantprec+1))
	// width : Anzahl Stellen, die die Zahl (inklusive Punkt) nicht
	//         überschreiten soll, oder 0
	// d : Mindestanzahl Nachkommastellen oder 0
	// k : Skalierungsfaktor (siehe CLTL S.394)
	// dmin : Mindestanzahl von Dezimaltellen, die (trotz Angabe von width
	//        oder d) nicht gerundet werden dürfen.
	//        (Nur interessant, falls d <= dmin <= (precision der Zahl).)
	// wandelt die Zahl significand*2^expon um in einen Dezimalstring um.
	// Es ist kein Exponent dabei.
	var cl_idecoded_float decoded = integer_decode_float(arg);
	var const cl_I& significand = decoded.mantissa;
	var const cl_I& expon = decoded.exponent;
	var uintC mantprec = float_digits(arg)-1;
	var cl_I numerator = significand;
	var cl_I denominator = 1;
	var cl_I abrund_einh = 1; // Abrundungseinheit:
	       // Abrunden um 1 in der letzten abrundbaren Stelle entspricht
	       // einer Erniedrigung von numerator um abrund_einh.
	var cl_I aufrund_einh = 1; // Aufrundungseinheit:
	       // Aufrunden um 1 in der letzten aufrundbaren Stelle entspricht
	       // einer Erhöhung von numerator um aufrund_einh.
	digitstring.reset();
	if (expon > 0) {
		numerator = numerator << expon;
		aufrund_einh = abrund_einh = 1 << expon;
	}
	elif (expon < 0) {
		denominator = denominator << -expon;
		// aufrund_einh = abrund_einh = 1;
	}
	// Zahl = numerator/denominator
	if (significand == ash(1,mantprec)) {
		// Ist der Significand=2^mantprec, so ist abrund-einh zu halbieren.
		// Man kann stattdessen auch alle 3 anderen Grössen verdoppeln:
		aufrund_einh = aufrund_einh << 1;
		numerator = numerator << 1;
		denominator = denominator << 1;
	}
	// Defaultmäßig: Auf-/Abrunde-Einheit = eine Einheit in der letzten
	// BINÄRstelle.
	// Zahl = numerator/denominator
	// Skalierungsfaktor k in die Zahl mit einbeziehen (vgl. CLTL S.394)
	// k<0 -> Mantisse durch 10^|k| dividieren
	// k>0 -> Mantisse mit 10^k multiplizieren
	// Dabei aufrund-einh, abrund-einh im Verhältnis zu numerator beibehalten.
	if (k != 0) {
		if (k < 0) {
			var cl_I skal_faktor = expt_pos(10,-k);
			denominator = denominator * skal_faktor;
		}
		elif (k > 0) {
			var cl_I skal_faktor = expt_pos(10,k);
			numerator = numerator * skal_faktor;
			aufrund_einh = aufrund_einh * skal_faktor;
			abrund_einh = abrund_einh * skal_faktor;
		}
	}
	// Stellen: 0 = 1. Stelle vor dem Punkt, -1 = 1. Stelle nach dem Punkt.
	var sintL stelle = 0; // Stelle der als nächstes auszugebenden Ziffer
	// auf >= 1/10 adjustieren:
	// (jeweils numerator mit 10 multiplizieren, eine führende 0 mehr vorsehen)
	until (10*numerator >= denominator) {
		stelle = stelle-1;
		numerator = numerator * 10;
		aufrund_einh = aufrund_einh * 10;
		abrund_einh = abrund_einh * 10;
	}
	// stelle = Stelle der letzten führenden 0
	//        = 1 + Stelle der 1. signifikanten Ziffer
	//        oder =0, falls k>=0
	// Ausführung der Rundung:
	var bool letzte_stelle_p = false; // d oder width angegeben?
	var sintL letzte_stelle = 0; // falls d oder width angegeben waren:
				     // Stelle der letzten signifikanten Ziffer
	var bool halbzahlig = false; // zeigt an, ob hinten genau ein 0.500000 wegfällt
	do {
		// Solange das Ergebnis auch nach Aufrundung >= 1 bliebe,
		// eine Vorkommastelle mehr einplanen:
		until (((numerator << 1) + aufrund_einh) < (denominator << 1)) {
			denominator = denominator * 10;
			stelle = stelle+1;
		}
		// Falls d oder width angegeben:
		// letzte_stelle ausrechnen
		if (d != 0) {
			// Falls dmin angegeben: min(-d,-dmin) = -max(d,dmin).
			// Sonst -d.
			letzte_stelle = -d;
			if (dmin > 0)
				if (letzte_stelle > -dmin)
					letzte_stelle = -dmin;
			letzte_stelle_p = true;
		}
		elif (width > 0) {
			// Falls nicht d, nur width angegeben:
			if (stelle < 0)
				// Es kommen führende Nullen nach dem Punkt -> d:=width-1
				letzte_stelle = 1-width;
			else
				// Es kommen keine führenden Nullen nach dem Punkt ->
				// Es wird stelle Vorkommaziffern geben, d:=width-1-stelle
				letzte_stelle = 1+stelle-width;
			// also letzte_stelle = -(width-1 - max(stelle,0))
			// wieder dmin berücksichtigen:
			if (dmin > 0)
				if (letzte_stelle > -dmin)
					letzte_stelle = -dmin;
			letzte_stelle_p = true;
		}
		if (letzte_stelle_p) {
			var sintL ziffernzahl = letzte_stelle - stelle;
			// ziffernzahl = - Zahl signifikanter Stellen oder >=0.
			var cl_I dezimal_einh = denominator;
			// dezimal-einh := ceiling(dezimal_einh*expt(10,ziffernzahl))
			if (ziffernzahl > 0)
				dezimal_einh = dezimal_einh*expt_pos(10,ziffernzahl);
			elif (ziffernzahl < 0)
				dezimal_einh = ceiling1(dezimal_einh,expt_pos(10,-ziffernzahl));
			// dezimal-einh = Um wieviel numerator erhöht bzw. erniedigt werden
			// müßte, damit sich die Dezimaldarstellung um genau 1 an der
			// Position letzte_stelle verändert.
			if (abrund_einh < dezimal_einh)
				abrund_einh = dezimal_einh;
			if (aufrund_einh < dezimal_einh)
				aufrund_einh = dezimal_einh;
			// Jetzt darf auch um eine (halbe) DEZIMAL-Einheit gerundet werden.
			if (aufrund_einh == dezimal_einh)
				halbzahlig = true;
		}
	} until (((numerator << 1) + aufrund_einh) < (denominator << 1));
	// stelle = Position der ersten signifikanten Stelle + 1
	var uintL digit_count = 0; // Zahl der bisher in digit-string
	       // ausgegebenen Ziffern (exklusive den Punkt)
	var uintL point_pos = 0; // Punkt-Position = Zahl führender Stellen
			         // = Zahl der Ziffern vor dem Punkt
	// Führenden Punkt und nachfolgende Nullen ausgeben:
	if (stelle < 0) {
		digitstring.push('.');
		point_pos = digit_count;
		for (int i = -stelle; i >= 0; i--) {
			digitstring.push('0');
			digit_count++;
		}
	}
	// Ziffern der Mantisse ausgeben:
	var uintL digit; // die laufende Ziffer, >=0, <10
	var bool abrunden; // letzte Ziffer abzurunden?
	var bool aufrunden; // letzte Ziffer aufzurunden?
	for (;;) {
		if (stelle == 0) {
			digitstring.push('.');
			point_pos = digit_count;
		}
		stelle = stelle-1;
		var cl_I_div_t div = cl_divide(numerator*10,denominator);
		digit = cl_I_to_UL(div.quotient);
		numerator = div.remainder;
		abrund_einh = abrund_einh*10;
		aufrund_einh = aufrund_einh*10;
		abrunden = ((numerator<<1) < abrund_einh);
		aufrunden = (halbzahlig
			     ? (numerator<<1) >= (denominator<<1) - aufrund_einh
			     : (numerator<<1) > (denominator<<1) - aufrund_einh
			    );
		if (abrunden || aufrunden
		    || (letzte_stelle_p && (stelle <= letzte_stelle))
		   )
		   break;
		digitstring.push("0123456789"[digit]);
		digit_count++;
	}
	// letzte signifikante Ziffer ausgeben:
	if (letzte_stelle_p ? (stelle >= letzte_stelle) : true) {
		digit = (abrunden && !aufrunden ? digit :
			 aufrunden && !abrunden ? digit+1 :
			 (numerator<<1) <= denominator ? digit : digit+1);
		digitstring.push("0123456789"[digit]);
		digit_count++;
	}
	// Nachfolgende Nullen und Punkt ausgeben
	if (stelle >= 0) {
		for (int i = stelle; i >= 0; i--) {
			digitstring.push('0');
			digit_count++;
		}
		digitstring.push('.');
		point_pos = digit_count;
	}
	if (d != 0)
		for (int i = d - (digit_count - point_pos); i >= 0; i--) {
			digitstring.push('0');
			digit_count++;
		}
	return digits_with_dot(digitstring.contents(), digit_count+1,
			point_pos==0,
			point_pos==digit_count,
			point_pos
		);
}

}  // namespace cln

