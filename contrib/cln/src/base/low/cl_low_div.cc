// Low level: division.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/cl_low.h"


// Implementation.

#ifdef NEED_VAR_divu_16_rest
uint16 divu_16_rest;
#endif

#ifdef NEED_FUNCTION_divu_3216_1616_
uint16 divu_16_rest;
namespace cln {
#if 1
// Most processors have a good 32 by 32 bit division, use that.
uint16 divu_3216_1616_ (uint32 x, uint16 y)
{
	var uint16 q = floor(x,(uint32)y);
	divu_16_rest = x - (uint32)q * (uint32)y;
	return q;
}
#else
// On processors without a hardware division, we have to produce the quotient's
// bits individually. Basically the algorithm is like this:
//     q := 0; r := x;
//     if (r >= 2^15*y) { r -= 2^15*y; q |= 2^15; }
//     if (r >= 2^14*y) { r -= 2^14*y; q |= 2^14; }
//     ...
//     if (r >= 2^0*y) { r -= 2^0*y; q |= 2^0; }
// We don't want to shift 2^k*y to the right. Instead, we shift r to the left:
//     q := 0; r := x;
//     if (r >= 2^15*y) { r -= 2^15*y; q |= 2^15; }
//     if (2*r >= 2^15*y) { 2*r -= 2^15*y; q |= 2^14; }
//     ...
//     if (2^15*r >= 2^15*y) { 2^15*r -= 2^15*y; q |= 2^0; }
// In other terms:
//     q := 0; r := x; s := 2^15*y;
//     if (r >= s) { r -= s; q |= 2^15; }
//     r := 2*r;
//     if (r >= s) { r -= s; q |= 2^14; }
//     r := 2*r;
//     ...
//     r := 2*r;
//     if (r >= s) { r -= s; q |= 2^0; }
//     r := r >> 15;
// Now, we combine r and q into a single register. The bits of q won't disturb
// the "r >= s" comparisons because up to the last comparisons only max 15
// bits have accumulated, and s is a multiple of 2^15.
// We can thus shift r and q with a single instruction.
//     q := 0; r := x; s := 2^15*y;
//     if (r >= s) { r -= s; r := 2*r+1; } else r := 2*r;
//     if (r >= s) { r -= s; r := 2*r+1; } else r := 2*r;
//     ...
//     if (r >= s) { r -= s; r := 2*r+1; } else r := 2*r;
//     q := r & (2^16-1); r := r >> 16;
// Up to now, this is pretty standard.
// On most hardware the comparison "r >= s" already implies doing a subtraction
// r - s. It is wasteful to do the subtraction twice. Better do it once and
// test the carry afterwards:
//     q := 0; r := x; s := 2^15*y;
//     r -= s; if (no subcarry) { r := 2*r+1; } else { r := 2*r+2*s; }
//     r -= s; if (no subcarry) { r := 2*r+1; } else { r := 2*r+2*s; }
//     ...
//     r -= s; if (no subcarry) { r := 2*r+1; } else { r := 2*r+2*s; }
//     q := r & (2^16-1); r := r >> 16;
// In the case of carry we can combine the "+2*s" with the next "-s" operation.
// So this becomes "r := 2*r+s". But note about the carries: the case
// "(2*r+2*s)-s gives no subcarry" is equivalent to "2*r+s gives an addcarry".
// On most processors, subcarry and addcarry are the same bit. So we turn
// the subcarry into an addcarry by writing "r += -s" instead of "r -= s"
// (or vice versa: writing "2*r-(-s)" instead of "2*r+s").
//     q := 0; r := x; s := 2^15*y;
//     r += -s;
//     if (addcarry) { r := 2*r+1; r += -s; } else { r := 2*r+s; }
//     if (addcarry) { r := 2*r+1; r += -s; } else { r := 2*r+s; }
//     ...
//     if (addcarry) { r := 2*r+1; } else { r := 2*r+2*s; }
//     q := r & (2^16-1); r := r >> 16;
// This algorithm is implemented in cl_asm_arm.cc and (in slightly modified
// form) in cl_asm_sparc.cc.
#endif
}  // namespace cln
#endif

#ifdef NEED_FUNCTION_divu_3232_3232_
namespace cln {
// Dies dient nur noch als Hilfsfunktion für floorD().
// Die Rückgabe des Restes in divu_32_rest ist also hier nicht nötig.
uint32 divu_3232_3232_(uint32 x, uint32 y)
{
	var uint32 q;
	divu_3232_3232(x,y,q=,);
	return q;
}
}  // namespace cln
#endif

#ifdef NEED_VAR_divu_32_rest
uint32 divu_32_rest;
#endif

#ifdef NEED_FUNCTION_divu_6432_3232_
uint32 divu_32_rest;
namespace cln {
uint32 divu_6432_3232_(uint32 xhi, uint32 xlo, uint32 y)
// Methode:
// Wie UDS_divide mit intDsize=16, a_len=4, b_len=2.
{
    if (y <= (uint32)(bit(16)-1))
        // 48-durch-16-Bit-Division,
        // aufgebaut aus zwei 32-durch-16-Bit-Divisionen:
        { var uint16 q1;
          var uint16 q0;
          var uint16 r1;
          divu_3216_1616(highlow32(low16(xhi),high16(xlo)),y, q1=,r1=);
          divu_3216_1616(highlow32(r1,low16(xlo)),y, q0=,divu_32_rest=);
          return highlow32(q1,q0);
        }
    // y>=2^16
    {// y shiften:
      var uintL s = 0;
      while ((sint32)y >= 0) { y = y<<1; s++; }
      // x entsprechend shiften:
      if (!(s==0))
        { xhi = (xhi << s) | (xlo >> (32-s)); xlo = xlo << s; }
      // 64-durch-32-Bit-Division,
      // aufgebaut aus zwei 48-durch-32-Bit-Divisionen.
      // Methode für eine 48-durch-32-Bit-Division x/y mit 0 <= x < 2^16*y :
      // (beta = 2^n = 2^16, n = 16)
      // Wir wissen beta^2/2 <= y < beta^2, Quotient  q = floor(x/y) < beta.
      // Schreibe  x = beta*x1 + x0  mit  x1 := floor(x/beta)
      // und       y = beta*y1 + y0  mit  y1 := floor(y/beta)
      // und bilde den Näherungs-Quotienten floor(x1/y1)
      // oder (noch besser) floor(x1/(y1+1)).
      // Wegen 0 <= x1 < 2^(2n) und 0 < 2^(n-1) <= y1 < 2^n
      // und  x1/(y1+1) <= x/y < x1/(y1+1) + 2
      // (denn x1/(y1+1) = (x1*beta)/((y1+1)*beta) <= (x1*beta)/y <= x/y
      // und x/y - x1/(y1+1) = (x+x*y1-x1*y)/(y*(y1+1))
      // = (x+x0*y1-x1*y0)/(y*(y1+1)) <= (x+x0*y1)/(y*(y1+1))
      // <= x/(y*(y1+1)) + x0/y = (x/y)/(y1+1) + x0/y
      // <= 2^n/(2^(n-1)+1) + 2^n/2^(2n-1) = 2^n/(2^(n-1)+1) + 2^(1-n) < 2 )
      // gilt  floor(x1/(y1+1)) <= floor(x/y) <= floor(x1/(y1+1)) + 2  .
      // Man bildet also  q:=floor(x1/(y1+1))  (ein Shift um n Bit oder
      // eine (2n)-durch-n-Bit-Division, mit Ergebnis q <= floor(x/y) < beta)
      // und x-q*y und muß hiervon noch höchstens 2 mal y abziehen und q
      // incrementieren, um den Quotienten  q = floor(x/y)  und den Rest
      // x-floor(x/y)*y  der Division zu bekommen.
      { var uint16 y1_1 = high16(y)+1; // y1+1
        var uint16 q1;
        var uint16 q0;
        var uint32 r;
        // 2^16*xhi+high16(xlo) durch y dividieren:
       {var uint16 r16;
        var uint32 r2;
        if (y1_1==0)
          { q1 = high16(xhi); r16 = low16(xhi); }
          else
          { divu_3216_1616(xhi,y1_1, q1=,r16=); }
        // q1 = floor(xhi/(y1+1)), r16 = xhi - (y1+1)*q1 (>=0, <=y1)
        // Bilde r := (2^16*xhi+high16(xlo)) - y*q1
        //          = 2^16*(xhi-y1*q1) + high16(xlo) - y0*q1
        //          = 2^16*r16 + 2^16*q1 + high16(xlo) - y0*q1 (>=0)
        // Dies ist < 2^16*y1 + 2^32 <= y + 2^32 <= 3*y, kann überlaufen!
        r = highlow32(r16,high16(xlo)); // 2^16*r16 + high16(xlo) < 2^32
        r2 = highlow32_0(q1) - mulu16(low16(y),q1); // 2^16*q1 - y0*q1 < 2^32
        // 0 <= r+r2 < 3*y. Bei der Addition auf Carry testen!
        // Carry -> jedenfalls y <= r+r2 < y + 2^32 <= 3*y.
        // kein Carry -> jedenfalls 0 <= r+r2 < 2^32 <= 2*y.
        if ((r += r2) < r2) // addieren, r >= 2^32 ?
          { q1 += 1; r -= y; }
        // jetzt noch 0 <= r < 2^32 <= 2*y
        if (r >= y)
          { q1 += 1; r -= y; }
       }// Quotient q1, Rest r fertig.
        // 2^16*r+low16(xlo) durch y dividieren:
       {var uint16 r16;
        var uint32 r2;
        if (y1_1==0)
          { q0 = high16(r); r16 = low16(r); }
          else
          { divu_3216_1616(r,y1_1, q0=,r16=); }
        // q0 = floor(r/(y1+1)), r16 = r - (y1+1)*q0 (>=0, <=y1)
        // Bilde r := (2^16*r+low16(xlo)) - y*q0
        //          = 2^16*(r-y1*q0) + low16(xlo) - y0*q0
        //          = 2^16*r16 + 2^16*q0 + low16(xlo) - y0*q0 (>=0)
        // Dies ist < 2^16*y1 + 2^32 <= y + 2^32 <= 3*y, kann überlaufen!
        r = highlow32(r16,low16(xlo)); // 2^16*r16 + low16(xlo) < 2^32
        r2 = highlow32_0(q0) - mulu16(low16(y),q0); // 2^16*q0 - y0*q0 < 2^32
        // 0 <= r+r2 < 3*y. Bei der Addition auf Carry testen!
        // Carry -> jedenfalls y <= r+r2 < y + 2^32 <= 3*y.
        // kein Carry -> jedenfalls 0 <= r+r2 < 2^32 <= 2*y.
        if ((r += r2) < r2) // addieren, r >= 2^32 ?
          { q0 += 1; r -= y; }
        // jetzt noch 0 <= r < 2^32 <= 2*y
        if (r >= y)
          { q0 += 1; r -= y; }
       }// Quotient q0, Rest r fertig.
        divu_32_rest = r >> s; // Rest
        return highlow32(q1,q0); // Quotient
}   } }
}  // namespace cln
#endif

#ifdef NEED_VAR_divu_64_rest
uint64 divu_64_rest;
#endif

#ifdef NEED_FUNCTION_divu_6464_6464_
namespace cln {
uint64 divu_6464_6464_(uint64 x, uint64 y)
// Methode: (beta = 2^n = 2^32, n = 32)
// Falls y < beta, handelt es sich um eine 64-durch-32-Bit-Division.
// Falls y >= beta:
// Quotient  q = floor(x/y) < beta  (da 0 <= x < beta^2, y >= beta).
// y habe genau n+k Bits (1 <= k <= n), d.h. 2^(n+k-1) <= y < 2^(n+k).
// Schreibe  x = 2^k*x1 + x0  mit  x1 := floor(x/2^k)
// und       y = 2^k*y1 + y0  mit  y1 := floor(y/2^k)
// und bilde den Näherungs-Quotienten floor(x1/y1)
// oder (noch besser) floor(x1/(y1+1)).
// Wegen 0 <= x1 < 2^(2n) und 0 < 2^(n-1) <= y1 < 2^n
// und  x1/(y1+1) <= x/y < x1/(y1+1) + 2
// (denn x1/(y1+1) = (x1*2^k)/((y1+1)*2^k) <= (x1*2^k)/y <= x/y
// und x/y - x1/(y1+1) = (x+x*y1-x1*y)/(y*(y1+1))
// = (x+x0*y1-x1*y0)/(y*(y1+1)) <= (x+x0*y1)/(y*(y1+1))
// <= x/(y*(y1+1)) + x0/y
// <= 2^(2n)/(2^(n+k-1)*(2^(n-1)+1)) + 2^k/2^(n+k-1)
// = 2^(n-k+1)/(2^(n-1)+1) + 2^(1-n) <= 2^n/(2^(n-1)+1) + 2^(1-n) < 2 )
// gilt  floor(x1/(y1+1)) <= floor(x/y) <= floor(x1/(y1+1)) + 2  .
// Man bildet also  q:=floor(x1/(y1+1))  (ein Shift um n Bit oder
// eine (2n)-durch-n-Bit-Division, mit Ergebnis q <= floor(x/y) < beta)
// und x-q*y und muss hiervon noch höchstens 2 mal y abziehen und q
// incrementieren, um den Quotienten  q = floor(x/y)  und den Rest
// x-floor(x/y)*y  der Division zu bekommen.
{
  if (y <= (uint64)(((uint64)1<<32)-1))
    { var uint32 q1;
      var uint32 q0;
      var uint32 r1;
      divu_6432_3232(0,high32(x),y, q1 = , r1 = );
      divu_6432_3232(r1,low32(x),y, q0 = , divu_64_rest = );
      return highlow64(q1,q0);
    }
    else
    { var uint64 x1 = x; // x1 := x
      var uint64 y1 = y; // y1 := y
      var uint32 q;
      do { x1 = floor(x1,2); y1 = floor(y1,2); } // k erhöhen
         while (!(y1 <= (uint64)(((uint64)1<<32)-1))); // bis y1 < beta
      { var uint32 y2 = low32(y1)+1; // y1+1 bilden
        if (y2==0)
          { q = high32(x1); } // y1+1=beta -> ein Shift
          else
          { divu_6432_3232(high32(x1),low32(x1),y2,q=,); } // Division von x1 durch y1+1
      }
      // q = floor(x1/(y1+1))
      // x-q*y bilden (eine 32-mal-64-Bit-Multiplikation ohne Überlauf):
      x -= highlow64_0(mulu32_w(q,high32(y))); // q * high32(y) * beta
      // gefahrlos, da q*high32(y) <= q*y/beta <= x/beta < beta
      x -= mulu32_w(q,low32(y)); // q * low32(y)
      // gefahrlos, da q*high32(y)*beta + q*low32(y) = q*y <= x
      // Noch höchstens 2 mal y abziehen:
      if (x >= y)
        { q += 1; x -= y;
          if (x >= y)
            { q += 1; x -= y;
        }   }
      divu_64_rest = x;
      return (uint64)q;
   }
}
}  // namespace cln
#endif

#ifdef NEED_FUNCTION_divu_12864_6464_
namespace cln {
uint64 divu_12864_6464_(uint64 xhi, uint64 xlo, uint64 y)
// Methode:
// Wie UDS_divide mit intDsize=32, a_len=4, b_len=2.
{
    if (y <= (uint64)(bit(32)-1))
        // 96-durch-32-Bit-Division,
        // aufgebaut aus zwei 64-durch-32-Bit-Divisionen:
        { var uint32 q1;
          var uint32 q0;
          var uint32 r1;
          divu_6432_3232(low32(xhi),high32(xlo),y, q1=,r1=);
          divu_6432_3232(r1,low32(xlo),y, q0=,divu_64_rest=);
          return highlow64(q1,q0);
        }
    // y>=2^32
    {// y shiften:
      var uintL s = 0;
      while ((sint64)y >= 0) { y = y<<1; s++; }
      // x entsprechend shiften:
      if (!(s==0))
        { xhi = (xhi << s) | (xlo >> (64-s)); xlo = xlo << s; }
      // 128-durch-64-Bit-Division,
      // aufgebaut aus zwei 96-durch-64-Bit-Divisionen.
      // Methode für eine 96-durch-64-Bit-Division x/y mit 0 <= x < 2^32*y :
      // (beta = 2^n = 2^32, n = 32)
      // Wir wissen beta^2/2 <= y < beta^2, Quotient  q = floor(x/y) < beta.
      // Schreibe  x = beta*x1 + x0  mit  x1 := floor(x/beta)
      // und       y = beta*y1 + y0  mit  y1 := floor(y/beta)
      // und bilde den Näherungs-Quotienten floor(x1/y1)
      // oder (noch besser) floor(x1/(y1+1)).
      // Wegen 0 <= x1 < 2^(2n) und 0 < 2^(n-1) <= y1 < 2^n
      // und  x1/(y1+1) <= x/y < x1/(y1+1) + 2
      // (denn x1/(y1+1) = (x1*beta)/((y1+1)*beta) <= (x1*beta)/y <= x/y
      // und x/y - x1/(y1+1) = (x+x*y1-x1*y)/(y*(y1+1))
      // = (x+x0*y1-x1*y0)/(y*(y1+1)) <= (x+x0*y1)/(y*(y1+1))
      // <= x/(y*(y1+1)) + x0/y = (x/y)/(y1+1) + x0/y
      // <= 2^n/(2^(n-1)+1) + 2^n/2^(2n-1) = 2^n/(2^(n-1)+1) + 2^(1-n) < 2 )
      // gilt  floor(x1/(y1+1)) <= floor(x/y) <= floor(x1/(y1+1)) + 2  .
      // Man bildet also  q:=floor(x1/(y1+1))  (ein Shift um n Bit oder
      // eine (2n)-durch-n-Bit-Division, mit Ergebnis q <= floor(x/y) < beta)
      // und x-q*y und muß hiervon noch höchstens 2 mal y abziehen und q
      // incrementieren, um den Quotienten  q = floor(x/y)  und den Rest
      // x-floor(x/y)*y  der Division zu bekommen.
      { var uint32 y1_1 = high32(y)+1; // y1+1
        var uint32 q1;
        var uint32 q0;
        var uint64 r;
        // 2^32*xhi+high32(xlo) durch y dividieren:
       {var uint32 r32;
        var uint64 r2;
        if (y1_1==0)
          { q1 = high32(xhi); r32 = low32(xhi); }
          else
          { divu_6432_3232_w(xhi,y1_1, q1=,r32=); }
        // q1 = floor(xhi/(y1+1)), r32 = xhi - (y1+1)*q1 (>=0, <=y1)
        // Bilde r := (2^32*xhi+high32(xlo)) - y*q1
        //          = 2^32*(xhi-y1*q1) + high32(xlo) - y0*q1
        //          = 2^32*r32 + 2^32*q1 + high32(xlo) - y0*q1 (>=0)
        // Dies ist < 2^32*y1 + 2^64 <= y + 2^64 <= 3*y, kann überlaufen!
        r = highlow64(r32,high32(xlo)); // 2^32*r32 + high32(xlo) < 2^64
        r2 = highlow64_0(q1) - mulu32_w(low32(y),q1); // 2^32*q1 - y0*q1 < 2^64
        // 0 <= r+r2 < 3*y. Bei der Addition auf Carry testen!
        // Carry -> jedenfalls y <= r+r2 < y + 2^64 <= 3*y.
        // kein Carry -> jedenfalls 0 <= r+r2 < 2^64 <= 2*y.
        if ((r += r2) < r2) // addieren, r >= 2^64 ?
          { q1 += 1; r -= y; }
        // jetzt noch 0 <= r < 2^64 <= 2*y
        if (r >= y)
          { q1 += 1; r -= y; }
       }// Quotient q1, Rest r fertig.
        // 2^32*r+low32(xlo) durch y dividieren:
       {var uint32 r32;
        var uint64 r2;
        if (y1_1==0)
          { q0 = high32(r); r32 = low32(r); }
          else
          { divu_6432_3232_w(r,y1_1, q0=,r32=); }
        // q0 = floor(r/(y1+1)), r32 = r - (y1+1)*q0 (>=0, <=y1)
        // Bilde r := (2^32*r+low32(xlo)) - y*q0
        //          = 2^32*(r-y1*q0) + low32(xlo) - y0*q0
        //          = 2^32*r32 + 2^32*q0 + low32(xlo) - y0*q0 (>=0)
        // Dies ist < 2^32*y1 + 2^64 <= y + 2^64 <= 3*y, kann überlaufen!
        r = highlow64(r32,low32(xlo)); // 2^32*r32 + low32(xlo) < 2^64
        r2 = highlow64_0(q0) - mulu32_w(low32(y),q0); // 2^32*q0 - y0*q0 < 2^64
        // 0 <= r+r2 < 3*y. Bei der Addition auf Carry testen!
        // Carry -> jedenfalls y <= r+r2 < y + 2^64 <= 3*y.
        // kein Carry -> jedenfalls 0 <= r+r2 < 2^64 <= 2*y.
        if ((r += r2) < r2) // addieren, r >= 2^64 ?
          { q0 += 1; r -= y; }
        // jetzt noch 0 <= r < 2^64 <= 2*y
        if (r >= y)
          { q0 += 1; r -= y; }
       }// Quotient q0, Rest r fertig.
        divu_64_rest = r >> s; // Rest
        return highlow64(q1,q0); // Quotient
}   } }
}  // namespace cln
#endif

