// logp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "integer/cl_I.h"
#include "rational/cl_RA.h"
#include "base/cl_xmacros.h"

namespace cln {

bool logp (const cl_I& a, const cl_I& b, cl_RA* l)
{
// Methode:
//   log(a,b) soll Bruch c/d mit teilerfremdem c>=0,d>0 ergeben.
//   a=1 -> c=0, d=1.
//   a>=b -> Dividiere a durch b. Rest da -> geht nicht.
//           Sonst log(a,b) = 1+log(a/b,b).
//           Berechne also c/d := log(a/b,b) und setze c:=c+d.
//   1<a<b -> log(a,b) = 1/log(b,a).
//           Berechne also c/d := log(b,a) und vertausche c und d.
// Man konstruiert hierbei eigentlich die Kettenbruchentwicklung von c/d.
// Wegen a>=2^c, b>=2^d sind c,d < (integer-length a,b) < intDsize*2^intCsize.
// In Matrizenschreibweise:
//   Wenn eine Folge von Divisionsschritten D und Vertauschungsschritten V
//   ausgeführt werden muß, z.B. (a,b) V D D = (1,*), so ist
//     ( c )           ( 0 )             ( 1 1 )           ( 0 1 )
//     ( d )  =  V D D ( 1 )  wobei  D = ( 0 1 )  und  V = ( 1 0 ).
//   Man baut diese Matrizen nun von links nach rechts auf, zum Schluß von
//              ( 0 )
//   rechts mit ( 1 ) multiplizieren.
// Entrekursiviert:
//   Wir werden (a,b) und damit auch c/d = log(a/b) verändern.
//   Invariante: Statt (c,d) wollen wir (uc*c+ud*d,vc*c+vd*d) zurückliefern.
//                                           ( uc ud )
//   D.h. die bisherige Matrix von links ist ( vc vd ).
//   uc:=1, ud:=0, vc:=0, vd:=1.
//   Solange a>1,
//     a>=b -> Dividiere a durch b. Rest da -> geht nicht.
//             Sonst a:=a/b, und (für später c:=c+d) ud:=uc+ud, vd:=vc+vd.
//     1<a<b -> vertausche a und b, uc und ud, vc und vd.
//   Liefere (ud,vd), der Bruch ud/vd ist gekürzt.
 {	Mutable(cl_I,a);
	Mutable(cl_I,b);
	var uintL uc = 1;
	var uintL ud = 0;
	var uintL vc = 0;
	var uintL vd = 1;
	loop {
		if (eq(a,1)) // a=1 -> Rekursion zu Ende
			break;
		if (a >= b) {
			var cl_I_div_t div = cl_divide(a,b); // a durch b dividieren
			if (!eq(div.remainder,0)) // Rest /=0 ?
				return false; // -> fertig
			a = div.quotient;  // a := a/b
			ud = uc + ud; vd = vc + vd;
		} else {
			// 1<a<b -> a und b vertauschen
			swap(cl_I, a, b);
			swap(uintL, uc, ud); swap(uintL, vc, vd);
		}
	}
	// a=1 -> c=0,d=1 -> Ergebnis ud/vd
	*l = I_I_to_RA(UL_to_I(ud),UL_to_I(vd)); return true;
}}

}  // namespace cln
