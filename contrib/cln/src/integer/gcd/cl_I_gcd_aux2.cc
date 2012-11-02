// partial_gcd().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/integer.h"
#include "base/digit/cl_D.h"

namespace cln {

// Dasselbe wie partial_gcd(z1,z2,erg), nur daß z1 und z2 Doppelworte sind.
// Bevor im Ergebnis erg ein Überlauf eintritt, wird abgebrochen.

#if HAVE_DD

static inline uintDD muluDD_unchecked(uintD q, uintDD a)
{
	return muluD(q,lowD(a)) + highlowDD_0(muluD_unchecked(q,highD(a)));
}

// Division: liefert min(floor(x / y), 2^intDsize-1).
// Vorausgesetzt wird, daß  x >= 2 * y > 0.
static uintD floorDD (uintDD x, uintDD y)
{
	// vgl. Algorithmus für divu_3232_3232().
	var uintD q;
	if (y < ((uintDD)1 << intDsize)) {
		if (highD(x) >= y)
			q = ~(uintD)0; // instead of overflow
		else
			divuD(x, (uintD)y, q =, );
		return q;
	}
	{
		var uintC shift;
		integerlengthD(highD(y), shift=);
		// NB: 0 < shift < intDsize since 2^intDsize <= y < 2^(2*intDsize-1).
		// Determine q := floor((x>>shift) / ((y>>shift)+1)).
		var uintD y_shifted = y >> shift;
		y_shifted += 1;
		if (y_shifted == 0)
			q = highD(x) >> shift;
		else
			divuD(highD(x) >> shift, y_shifted, q =, );
	}
	// May need to increment q at most twice.
	{
		var uintDD p = muluDD_unchecked(q,y);
		#ifdef DEBUG_GCD
		if (x < p)
			throw runtime_exception();
		#endif
		x -= p;
	}
	if (x >= y) {
		q += 1;
		x -= y;
		if (x >= y) {
			q += 1;
			#ifdef DEBUG_GCD
			x -= y;
			if (x >= y)
				throw runtime_exception();
			#endif
		}
	}
	return q;
}

void partial_gcd (uintDD z1, uintDD z2, partial_gcd_result* erg)
{
    var uintD x1 = 1;
    var uintD y1 = 0;
    var uintD x2 = 0;
    var uintD y2 = 1;
    for (;;) {
	// Hier ist z1-y1>=z2+y2.
	// Bestimme q := floor((z1-y1)/(z2+y2)) >= 1 :
	{
		var uintDD zaehler = z1 - (uintDD)y1;
		var uintDD nenner = z2 + (uintDD)y2;
		// z2+y2 <= z1-y1 < beta^2 !
		if (x2 > (~x1) >> 3) // x1 + 8*x2 >= beta ?
			goto do_subtract_1;
		if (y2 > (~y1) >> 3) // y1 + 8*y2 >= beta ?
			goto do_subtract_1;
		// Ist floor(zaehler,8) >= nenner ? Wenn nein -> do_subtract_1.
		if ((zaehler >> 3) < nenner)
			goto do_subtract_1;
		if (1) {
			var uintD q = floorDD(zaehler,nenner);
		    repeat_1:
			var uintDD qx = muluD(q,x2);
			if (qx > (uintDD)(~x1)) {
				// Choose a smaller value for q, to avoid overflow of x1.
				q = floorD(~x1,x2);
				goto repeat_1;
			}
			var uintDD qy = muluD(q,y2);
			if (qy > (uintDD)(~y1)) {
				// Choose a smaller value for q, to avoid overflow of y1.
				q = floorD(~y1,y2);
				goto repeat_1;
			}
			x1 += (uintD)qx;
			y1 += (uintD)qy;
			z1 -= muluDD_unchecked(q,z2);
		} else {
		    do_subtract_1:
			do {
				// Checks to avoid overflow.
				if (x2 > ~x1) goto done;
				if (y2 > ~y1) goto done;
				#ifdef DEBUG_GCD
				if (z1 < z2) throw runtime_exception();
				#endif
				// Now really subtract.
				x1 += x2;
				y1 += y2;
				z1 -= z2;
			} while (z1 - (uintDD)y1 >= nenner);
		}
	}
	if (z2 - (uintDD)x2 <= z1 + (uintDD)(x1-1)) goto done;
	// Hier ist z2-x2>=z1+x1.
	// Bestimme q := floor((z2-x2)/(z1+x1)) >= 1 :
	{
		var uintDD zaehler = z2 - (uintDD)x2;
		var uintDD nenner = z1 + (uintDD)x1;
		// z1+x1 <= z2-x2 < beta^2 !
		if (x1 > (~x2) >> 3) // x2 + 8*x1 >= beta ?
			goto do_subtract_2;
		if (y1 > (~y2) >> 3) // y2 + 8*y1 >= beta ?
			goto do_subtract_2;
		// Ist floor(zaehler,8) >= nenner ? Wenn nein -> do_subtract_2.
		if ((zaehler >> 3) < nenner)
			goto do_subtract_2;
		if (1) {
			var uintD q = floorDD(zaehler,nenner);
		    repeat_2:
			var uintDD qx = muluD(q,x1);
			if (qx > (uintDD)(~x2)) {
				// Choose a smaller value for q, to avoid overflow of x2.
				q = floorD(~x2,x1);
				goto repeat_2;
			}
			var uintDD qy = muluD(q,y1);
			if (qy > (uintDD)(~y2)) {
				// Choose a smaller value for q, to avoid overflow of y2.
				q = floorD(~y2,y1);
				goto repeat_2;
			}
			x2 += (uintD)qx;
			y2 += (uintD)qy;
			z2 -= muluDD_unchecked(q,z1);
		} else {
		    do_subtract_2:
			do {
				// Checks to avoid overflow.
				if (x1 > ~x2) goto done;
				if (y1 > ~y2) goto done;
				#ifdef DEBUG_GCD
				if (z2 < z1) throw runtime_exception();
				#endif
				// Now really subtract.
				x2 += x1;
				y2 += y1;
				z2 -= z1;
			} while (z2 - (uintDD)x2 >= nenner);
		}
	}
	if (z1 - (uintDD)y1 <= z2 + (uintDD)(y2-1)) goto done;
    }
    done:
	// Keine Subtraktion (ohne Überlauf) mehr möglich.
	erg->x1 = x1; erg->y1 = y1; erg->x2 = x2; erg->y2 = y2; // Ergebnis
}

#else

// Division: liefert min(floor(xhi|xlo / yhi|ylo), 2^intDsize-1).
// Vorausgesetzt wird, daß  xhi|xlo >= 2 * yhi|ylo > 0.
static uintD floorDD (uintD xhi, uintD xlo, uintD yhi, uintD ylo)
{
	// vgl. Algorithmus für divu_3232_3232().
	var uintD q;
	if (yhi == 0) {
		if (xhi >= ylo)
			q = ~(uintD)0; // instead of overflow
		else
			divuD(xhi,xlo, ylo, q =, );
		return q;
	}
	{
		var uintC shift;
		integerlengthD(yhi, shift=);
		// NB: 0 < shift < intDsize since 2^intDsize <= y < 2^(2*intDsize-1).
		// Determine q := floor((x>>shift) / ((y>>shift)+1)).
		var uintD y_shifted = (uintD)(yhi << (intDsize-shift)) | (ylo >> shift);
		y_shifted += 1;
		if (y_shifted == 0)
			q = xhi >> shift;
		else
			divuD(xhi >> shift, (uintD)(xhi << (intDsize-shift)) | (xlo >> shift),
			      y_shifted,
			      q =, );
	}
	// May need to increment q at most twice.
	{
		var uintD phi;
		var uintD plo;
		muluD(q,ylo, phi =, plo =);
		muluD(q,yhi,      , phi +=);
		#ifdef DEBUG_GCD
		if ((xhi < phi) || ((xhi == phi) && (xlo < plo)))
			throw runtime_exception();
		#endif
		xhi = xhi - phi;
		if (xlo < plo)
			xhi -= 1;
		xlo = xlo - plo;
	}
	if ((xhi > yhi) || ((xhi == yhi) && (xlo >= ylo))) {
		q += 1;
		xhi = xhi - yhi;
		if (xlo < ylo)
			xhi -= 1;
		xlo = xlo - ylo;
		if ((xhi > yhi) || ((xhi == yhi) && (xlo >= ylo))) {
			q += 1;
			#ifdef DEBUG_GCD
			xhi = xhi - yhi;
			if (xlo < ylo)
				xhi -= 1;
			xlo = xlo - ylo;
			if ((xhi > yhi) || ((xhi == yhi) && (xlo >= ylo)))
				throw runtime_exception();
			#endif
		}
	}
	return q;
}

void partial_gcd (uintD z1hi, uintD z1lo, uintD z2hi, uintD z2lo, partial_gcd_result* erg)
{
    var uintD x1 = 1;
    var uintD y1 = 0;
    var uintD x2 = 0;
    var uintD y2 = 1;
    for (;;) {
	// Hier ist z1-y1>=z2+y2.
	// Bestimme q := floor((z1-y1)/(z2+y2)) >= 1 :
	{
		var uintD zaehlerhi = z1hi;
		var uintD zaehlerlo = z1lo - y1;
		if (zaehlerlo > z1lo) { zaehlerhi--; }
		var uintD nennerhi = z2hi;
		var uintD nennerlo = z2lo + y2;
		if (nennerlo < z2lo) { nennerhi++; }
		// z2+y2 <= z1-y1 < beta^2 !
		if (x2 > (~x1) >> 3) // x1 + 8*x2 >= beta ?
			goto do_subtract_1;
		if (y2 > (~y1) >> 3) // y1 + 8*y2 >= beta ?
			goto do_subtract_1;
		// Ist floor(zaehler,8) >= nenner ? Wenn nein -> do_subtract_1.
		if ((zaehlerhi >> 3) < nennerhi)
			goto do_subtract_1;
		if ((zaehlerhi >> 3) == nennerhi)
			if (((uintD)(zaehlerhi << (intDsize-3)) | (zaehlerlo >> 3)) < nennerlo)
				goto do_subtract_1;
		if (1) {
			var uintD q = floorDD(zaehlerhi,zaehlerlo,nennerhi,nennerlo);
		    repeat_1:
			var uintD qx;
			var uintD qy;
			{
				var uintD qxhi;
				muluD(q,x2, qxhi =, qx =);
				if ((qxhi > 0) || (qx > ~x1)) {
					// Choose a smaller value for q, to avoid overflow of x1.
					q = floorD(~x1,x2);
					goto repeat_1;
				}
			}
			{
				var uintD qyhi;
				muluD(q,y2, qyhi =, qy =);
				if ((qyhi > 0) || (qy > ~y1)) {
					// Choose a smaller value for q, to avoid overflow of y1.
					q = floorD(~y1,y2);
					goto repeat_1;
				}
			}
			x1 += qx;
			y1 += qy;
			{
				var uintD qzhi;
				var uintD qzlo;
				muluD(q,z2lo, qzhi =, qzlo =);
				muluD(q,z2hi,       , qzhi +=);
				z1hi -= qzhi;
				if (z1lo < qzlo)
					z1hi -= 1;
				z1lo -= qzlo;
			}
		} else {
		    do_subtract_1:
			for (;;) {
				// Checks to avoid overflow.
				if (x2 > ~x1) goto done;
				if (y2 > ~y1) goto done;
				#ifdef DEBUG_GCD
				if (z1hi < z2hi) throw runtime_exception();
				if (z1hi == z2hi) if (z1lo < z2lo) throw runtime_exception();
				#endif
				// Now really subtract.
				x1 += x2;
				y1 += y2;
				z1hi -= z2hi;
				if (z1lo < z2lo)
					z1hi -= 1;
				z1lo -= z2lo;
				var uintD z1dec_hi = z1hi;
				var uintD z1dec_lo = z1lo - y1;
				if (z1lo < y1)
					z1dec_hi -= 1;
				if (z1dec_hi < nennerhi)
					break;
				if (z1dec_hi == nennerhi)
					if (z1dec_lo < nennerlo)
						break;
			}
		}
	}
	{
		var uintD z1inc_hi = z1hi;
		var uintD z1inc_lo = z1lo + x1-1;
		if (z1inc_lo < z1lo)
			z1inc_hi += 1;
		var uintD z2dec_hi = z2hi;
		var uintD z2dec_lo = z2lo - x2;
		if (z2dec_lo > z2lo)
			z2dec_hi -= 1;
		if (z2dec_hi < z1inc_hi) goto done;
		if (z2dec_hi == z1inc_hi) if (z2dec_lo <= z1inc_lo) goto done;
	}
	// Hier ist z2-x2>=z1+x1.
	// Bestimme q := floor((z2-x2)/(z1+x1)) >= 1 :
	{
		var uintD zaehlerhi = z2hi;
		var uintD zaehlerlo = z2lo - x2;
		if (zaehlerlo > z2lo) { zaehlerhi--; }
		var uintD nennerhi = z1hi;
		var uintD nennerlo = z1lo + x1;
		if (nennerlo < z1lo) { nennerhi++; }
		// z1+x1 <= z2-x2 < beta^2 !
		if (x1 > (~x2) >> 3) // x2 + 8*x1 >= beta ?
			goto do_subtract_2;
		if (y1 > (~y2) >> 3) // y2 + 8*y1 >= beta ?
			goto do_subtract_2;
		// Ist floor(zaehler,8) >= nenner ? Wenn nein -> do_subtract_2.
		if ((zaehlerhi >> 3) < nennerhi)
			goto do_subtract_2;
		if ((zaehlerhi >> 3) == nennerhi)
			if (((uintD)(zaehlerhi << (intDsize-3)) | (zaehlerlo >> 3)) < nennerlo)
				goto do_subtract_2;
		if (1) {
			var uintD q = floorDD(zaehlerhi,zaehlerlo,nennerhi,nennerlo);
		    repeat_2:
			var uintD qx;
			var uintD qy;
			{
				var uintD qxhi;
				muluD(q,x1, qxhi =, qx =);
				if ((qxhi > 0) || (qx > ~x2)) {
					// Choose a smaller value for q, to avoid overflow of x2.
					q = floorD(~x2,x1);
					goto repeat_2;
				}
			}
			{
				var uintD qyhi;
				muluD(q,y1, qyhi =, qy =);
				if ((qyhi > 0) || (qy > ~y2)) {
					// Choose a smaller value for q, to avoid overflow of y2.
					q = floorD(~y2,y1);
					goto repeat_2;
				}
			}
			x2 += qx;
			y2 += qy;
			{
				var uintD qzhi;
				var uintD qzlo;
				muluD(q,z1lo, qzhi =, qzlo =);
				muluD(q,z1hi,       , qzhi +=);
				z2hi -= qzhi;
				if (z2lo < qzlo)
					z2hi -= 1;
				z2lo -= qzlo;
			}
		} else {
		    do_subtract_2:
			for (;;) {
				// Checks to avoid overflow.
				if (x1 > ~x2) goto done;
				if (y1 > ~y2) goto done;
				#ifdef DEBUG_GCD
				if (z2hi < z1hi) throw runtime_exception();
				if (z2hi == z1hi) if (z2lo < z1lo) throw runtime_exception();
				#endif
				// Now really subtract.
				x2 += x1;
				y2 += y1;
				z2hi -= z1hi;
				if (z2lo < z1lo)
					z2hi -= 1;
				z2lo -= z1lo;
				var uintD z2dec_hi = z2hi;
				var uintD z2dec_lo = z2lo - x2;
				if (z2lo < x2)
					z2dec_hi -= 1;
				if (z2dec_hi < nennerhi)
					break;
				if (z2dec_hi == nennerhi)
					if (z2dec_lo < nennerlo)
						break;
			}
		}
	}
	{
		var uintD z2inc_hi = z2hi;
		var uintD z2inc_lo = z2lo + y2-1;
		if (z2inc_lo < z2lo)
			z2inc_hi += 1;
		var uintD z1dec_hi = z1hi;
		var uintD z1dec_lo = z1lo - y1;
		if (z1dec_lo > z1lo)
			z1dec_hi -= 1;
		if (z1dec_hi < z2inc_hi) goto done;
		if (z1dec_hi == z2inc_hi) if (z1dec_lo <= z2inc_lo) goto done;
	}
    }
    done:
	// Keine Subtraktion (ohne Überlauf) mehr möglich.
	erg->x1 = x1; erg->y1 = y1; erg->x2 = x2; erg->y2 = y2; // Ergebnis
}

#endif

}  // namespace cln
