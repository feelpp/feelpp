// Internals for transcendental functions on floating-point numbers

#ifndef _CL_F_TRAN_H
#define _CL_F_TRAN_H

#include "cln/number.h"
#include "cln/float.h"

namespace cln {

// pi.
extern const cl_SF& cl_SF_pi();
extern const cl_FF& cl_FF_pi();
extern const cl_DF& cl_DF_pi();
extern cl_LF& cl_LF_pi(); // as long as it has ever been computed
extern const cl_LF pi (uintC len); // computes it even further

// cl_exp_aux(p,lq,len) liefert die Zahl exp(p/2^lq) mit len Digits.
// 0 < |p| < 2^lq.
// Es sollte |p|^2 < 2^lq sein, sonst ist das nicht effizient.
extern const cl_LF cl_exp_aux (const cl_I& p, uintE lq, uintC len);

// cl_cossin_aux(p,lq,len) liefert cos(p/2^lq) und sin(p/2^lq) mit len Digits.
// 0 < |p| < 2^lq.
// Es sollte |p|^2 < 2^lq sein, sonst ist das nicht effizient.
struct cl_LF_cos_sin_t {
	cl_LF cos;
	cl_LF sin;
// Constructor:
	cl_LF_cos_sin_t (const cl_LF& u, const cl_LF& v) : cos (u), sin (v) {}
	cl_LF_cos_sin_t () {}
};
extern const cl_LF_cos_sin_t cl_cossin_aux (const cl_I& p, uintE lq, uintC len);

// cl_coshsinh_aux(p,lq,len) liefert cosh(p/2^lq) und sinh(p/2^lq) mit len Digits.
// 0 < |p| < 2^lq.
// Es sollte |p|^2 < 2^lq sein, sonst ist das nicht effizient.
struct cl_LF_cosh_sinh_t {
	cl_LF cosh;
	cl_LF sinh;
// Constructor:
	cl_LF_cosh_sinh_t (const cl_LF& u, const cl_LF& v) : cosh (u), sinh (v) {}
	cl_LF_cosh_sinh_t () {}
};
extern const cl_LF_cosh_sinh_t cl_coshsinh_aux (const cl_I& p, uintE lq, uintC len);

// atanhx(x) liefert zu einem Float x (betragsmäßig <1/2) atanh(x) als Float.
extern const cl_F atanhx (const cl_F& x);

// atanx(x) liefert zu einem Float x (betragsmäßig <=1) atan(x) als Float.
extern const cl_F atanx (const cl_F& x);

// sinx(x) liefert zu einem Float x (betragsmäßig <1) sin(x)^2 als Float.
// sinxbyx(x) liefert zu einem Float x (betragsmäßig <1) (sin(x)/x)^2 als Float.
extern const cl_LF sinx_naive (const cl_LF& x); // requires cl_F_extendsqrt
extern const cl_F sinxbyx_naive (const cl_F& x); // requires cl_F_extendsqrt
// (cos(x),sin(x)) für ein Long-Float x (betragsmäßig <1).
extern const cl_LF_cos_sin_t cl_cossin_ratseries (const cl_LF& x); // requires extend by 1

// sinhx(x) liefert zu einem Float x (betragsmäßig <1) sinh(x)^2 als Float.
// sinhxbyx(x) liefert zu einem Float x (betragsmäßig <1) (sinh(x)/x)^2 als Float.
extern const cl_LF sinhx_naive (const cl_LF& x); // requires cl_F_extendsqrt
extern const cl_F sinhxbyx_naive (const cl_F& x); // requires cl_F_extendsqrt
// (cosh(x),sinh(x)) für ein Long-Float x (betragsmäßig <1).
extern const cl_LF_cosh_sinh_t cl_coshsinh_ratseries (const cl_LF& x); // requires extend by 1

// cl_round_pi(x) dividiert ein Float x mit Rest durch pi.
// Beide Werte von (round x (float pi x)).
extern const cl_F_div_t cl_round_pi (const cl_F& x);

// cl_round_pi2(x) dividiert ein Float x mit Rest durch pi/2.
// Beide Werte von (round x (float pi/2 x)).
extern const cl_F_div_t cl_round_pi2 (const cl_F& x);

// cl_atan_recip(m,len) liefert arctan(1/m) mit len Digits.
extern const cl_LF cl_atan_recip (cl_I m, uintC len);

// lnx(x) liefert zu einem Float x (>=1/2, <=2) ln(x) als Float.
extern const cl_F lnx_naive (const cl_F& x); // requires cl_F_extendsqrtx
extern const cl_LF lnx_naive (const cl_LF& x); // requires cl_F_extendsqrtx
extern const cl_LF lnx_ratseries (const cl_LF& x); // requires extend by 1

// cl_atanh_recip(m,len) liefert artanh(1/m) mit len Digits.
extern const cl_LF cl_atanh_recip (cl_I m, uintC len);

// ln(2).
extern const cl_SF& cl_SF_ln2();
extern const cl_FF& cl_FF_ln2();
extern const cl_DF& cl_DF_ln2();
extern cl_LF& cl_LF_ln2(); // as long as it has ever been computed
extern const cl_LF cl_ln2 (uintC len); // computes it even further

// cl_ln2(y) liefert die Zahl ln(2) im selben Float-Format wie y.
// > y: ein Float
extern const cl_F cl_ln2 (const cl_F& y);

// cl_ln2(y) liefert die Zahl ln(2) im Float-Format f.
// > f: eine Float-Format-Spezifikation
extern const cl_F cl_ln2 (float_format_t f);

// ln(10).
extern const cl_SF& cl_SF_ln10();
extern const cl_FF& cl_FF_ln10();
extern const cl_DF& cl_DF_ln10();
extern cl_LF& cl_LF_ln10(); // as long as it has ever been computed
extern const cl_LF cl_ln10 (uintC len); // computes it even further

// cl_ln10(y) liefert die Zahl ln(10) im selben Float-Format wie y.
// > y: ein Float
extern const cl_F cl_ln10 (const cl_F& y);

// cl_ln10(y) liefert die Zahl ln(10) im Float-Format f.
// > f: eine Float-Format-Spezifikation
extern const cl_F cl_ln10 (float_format_t f);

// e = exp(1).
extern const cl_SF& cl_SF_exp1();
extern const cl_FF& cl_FF_exp1();
extern const cl_DF& cl_DF_exp1();
extern cl_LF& cl_LF_exp1(); // as long as it has ever been computed
extern const cl_LF exp1 (uintC len); // computes it even further

// expx(x) liefert zu einem Float x (betragsmäßig <1) exp(x) als Float.
extern const cl_F expx_naive (const cl_F& x); // requires cl_F_extendsqrtx
extern const cl_LF expx_naive (const cl_LF& x); // requires cl_F_extendsqrtx
extern const cl_LF expx_ratseries (const cl_LF& x); // requires extend by 1

// Eulersche Konstante.
extern const cl_SF& cl_SF_eulerconst();
extern const cl_FF& cl_FF_eulerconst();
extern const cl_DF& cl_DF_eulerconst();
extern cl_LF& cl_LF_eulerconst(); // as long as it has ever been computed
extern const cl_LF eulerconst (uintC len); // computes it even further

// Catalansche Konstante.
extern const cl_SF& cl_SF_catalanconst();
extern const cl_FF& cl_FF_catalanconst();
extern const cl_DF& cl_DF_catalanconst();
extern cl_LF& cl_LF_catalanconst(); // as long as it has ever been computed
extern const cl_LF catalanconst (uintC len); // computes it even further

// Zeta-Funktion für s>1 ganzzahlig.
extern const cl_LF zeta (int s, uintC len);
// Zeta-Funktion für s=3.
extern const cl_LF zeta3 (uintC len);

}  // namespace cln

#endif /* _CL_F_TRAN_H */
