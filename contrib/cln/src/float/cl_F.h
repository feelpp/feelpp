// cl_F internals

#ifndef _CL_F_H
#define _CL_F_H

#include "cln/number.h"
#include "base/cl_macros.h"
#include "cln/float.h"

namespace cln {

#define underflow_allowed()  (! cl_inhibit_floating_point_underflow)


// For all floating-point formats:
// Sign s, Exponent e, Mantissa mk-1,...,m0
// represents the number (-1)^s * 2^(e-_EXP_MID) * [0 . 1 mk-1 ... m0]
// e=0 represents the number 0, always with sign s=0 (and mantissa =0).
// _exp_low and _exp_high are (inclusive) bounds for e.
// Bits for   Sign s    Exponent e    Mantissa m (= k)
// SF           1           8             16
// FF           1           8             23
// DF           1          11             52
// LF           1       32 or 64      intDsize*n >= 53


// Konversionen ohne Rundung:

// cl_SF_to_FF(x) wandelt ein Short-Float x in ein Single-Float um.
extern const cl_FF cl_SF_to_FF (const cl_SF& x);

// cl_SF_to_DF(x) wandelt ein Short-Float x in ein Double-Float um.
extern const cl_DF cl_SF_to_DF (const cl_SF& x);

// cl_SF_to_LF(x,len) wandelt ein Short-Float x in ein Long-Float mit len Digits um.
// > uintC len: gewünschte Anzahl Digits, >=LF_minlen
extern const cl_LF cl_SF_to_LF (const cl_SF& x, uintC len);

// cl_FF_to_DF(x) wandelt ein Single-Float x in ein Double-Float um.
extern const cl_DF cl_FF_to_DF (const cl_FF& x);

// cl_FF_to_LF(x,len) wandelt ein Single-Float x in ein Long-Float mit len Digits um.
// > uintC len: gewünschte Anzahl Digits, >=LF_minlen
extern const cl_LF cl_FF_to_LF (const cl_FF& x, uintC len);

// cl_DF_to_LF(x,len) wandelt ein Double-Float x in ein Long-Float mit len Digits um.
// > uintC len: gewünschte Anzahl Digits, >=LF_minlen
extern const cl_LF cl_DF_to_LF (const cl_DF& x, uintC len);


// Konversionen mit Rundung:

// cl_FF_to_SF(x) wandelt ein Single-Float x in ein Short-Float um.
extern const cl_SF cl_FF_to_SF (const cl_FF& x);

// cl_DF_to_SF(x) wandelt ein Double-Float x in ein Short-Float um.
extern const cl_SF cl_DF_to_SF (const cl_DF& x);

// cl_LF_to_SF(x) wandelt ein Long-Float x in ein Short-Float um.
extern const cl_SF cl_LF_to_SF (const cl_LF& x);

// cl_DF_to_FF(x) wandelt ein Double-Float x in ein Single-Float um.
extern const cl_FF cl_DF_to_FF (const cl_DF& x);

// cl_LF_to_FF(x) wandelt ein Long-Float x in ein Single-Float um.
extern const cl_FF cl_LF_to_FF (const cl_LF& x);

// cl_LF_to_DF(x) wandelt ein Long-Float x in ein Double-Float um.
extern const cl_DF cl_LF_to_DF (const cl_LF& x);


// Runtime typing support.
extern cl_class cl_class_ffloat;
extern cl_class cl_class_dfloat;
extern cl_class cl_class_lfloat;

// Type test.
inline bool longfloatp (const cl_F& x)
{
	if (x.pointer_p())
		if (x.pointer_type() == &cl_class_lfloat)
			return true;
	return false;
}

// Macro: verteilt je nach Float-Typ eines Floats x auf 4 Statements.
// floattypecase(x, SF_statement,FF_statement,DF_statement,LF_statement);
// x sollte eine Variable sein.
#ifdef CL_WIDE_POINTERS
  #define floattypecase(x, SF_statement,FF_statement,DF_statement,LF_statement) \
    if (!(x).pointer_p())						\
      switch ((x).nonpointer_tag())					\
        { case cl_SF_tag: { SF_statement } break;			\
          case cl_FF_tag: { FF_statement } break;			\
          default: NOTREACHED						\
        }								\
      else {								\
        if ((x).pointer_type() == &cl_class_dfloat) { DF_statement }	\
        else if ((x).pointer_type() == &cl_class_lfloat) { LF_statement } \
        else NOTREACHED							\
      }
#else
  #define floattypecase(x, SF_statement,FF_statement,DF_statement,LF_statement) \
    if (!(x).pointer_p())						\
      switch ((x).nonpointer_tag())					\
        { case cl_SF_tag: { SF_statement } break;			\
          default: NOTREACHED						\
        }								\
      else {								\
        if ((x).pointer_type() == &cl_class_ffloat) { FF_statement }	\
        else if ((x).pointer_type() == &cl_class_dfloat) { DF_statement } \
        else if ((x).pointer_type() == &cl_class_lfloat) { LF_statement } \
        else NOTREACHED							\
      }
#endif

// Macro: verteilt je nach Float-Typ eines Floats x auf 4 Statements,
// die x vom jeweiligen Float-Typ benutzen dürfen.
// floatcase(x, SF_statement,FF_statement,DF_statement,LF_statement);
// x sollte eine Variable sein.
  #define floatcase(x, SF_statement,FF_statement,DF_statement,LF_statement) \
    floattypecase(x							   \
      , var cl_SF& __tmp = *(cl_SF*)&x; var cl_SF& x = __tmp; SF_statement \
      , var cl_FF& __tmp = *(cl_FF*)&x; var cl_FF& x = __tmp; FF_statement \
      , var cl_DF& __tmp = *(cl_DF*)&x; var cl_DF& x = __tmp; DF_statement \
      , var cl_LF& __tmp = *(cl_LF*)&x; var cl_LF& x = __tmp; LF_statement \
      )


// GEN_F_OP1(arg1,F_OP,ergebnis_zuweisung)
// generates the body of a float operation with one argument.
// LF_OP is executed once the argument has been converted to its exact
// float type.
#define GEN_F_OP1(arg1,F_OP,ergebnis_zuweisung)  \
{									\
	floatcase(arg1							\
	, /* SF */	ergebnis_zuweisung F_OP(arg1);			\
	, /* FF */	ergebnis_zuweisung F_OP(arg1);			\
	, /* DF */	ergebnis_zuweisung F_OP(arg1);			\
	, /* LF */	ergebnis_zuweisung F_OP(arg1);			\
	);								\
}


// GEN_F_OP2(arg1,arg2,F_OP,r,s,ergebnis_zuweisung)
// generates the body of a float operation with two arguments.
// F_OP is executed once both arguments have been converted to the same
// float format (the longer one of arg1 and arg2). The r results are then
// converted the shorter of the two float formats. (r = 0,1,2.)
// s = 0,1. s=0 means the LF operation needs two long-floats of the same size.
// s=1 means they may be of different sizes.
#define GEN_F_OP2(arg1,arg2,F_OP,r,s,ergebnis_zuweisung)  \
{									\
	floatcase(arg1							\
	, /* arg1 SF */							\
		floatcase(arg2						\
		, /* arg2 SF */						\
			ergebnis_zuweisung CONCAT(NOMAP,r)(SF,		\
			F_OP(arg1,arg2) );				\
		, /* arg2 FF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(FF,cl_FF_to_SF,\
			F_OP(cl_SF_to_FF(arg1),arg2) );			\
		, /* arg2 DF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(DF,cl_DF_to_SF,\
			F_OP(cl_SF_to_DF(arg1),arg2) );			\
		, /* arg2 LF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(LF,cl_LF_to_SF,\
			F_OP(cl_SF_to_LF(arg1,CONCAT(LFlen,s)(arg2)),arg2) ); \
		);							\
	, /* arg1 FF */							\
		floatcase(arg2						\
		, /* arg2 SF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(FF,cl_FF_to_SF,\
			F_OP(arg1,cl_SF_to_FF(arg2)) );			\
		, /* arg2 FF */						\
			ergebnis_zuweisung CONCAT(NOMAP,r)(FF,		\
			F_OP(arg1,arg2) );				\
		, /* arg2 DF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(DF,cl_DF_to_FF,\
			F_OP(cl_FF_to_DF(arg1),arg2) );			\
		, /* arg2 LF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(LF,cl_LF_to_FF,\
			F_OP(cl_FF_to_LF(arg1,CONCAT(LFlen,s)(arg2)),arg2) ); \
		);							\
	, /* arg1 DF */							\
		floatcase(arg2						\
		, /* arg2 SF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(DF,cl_DF_to_SF,\
			F_OP(arg1,cl_SF_to_DF(arg2)) );			\
		, /* arg2 FF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(DF,cl_DF_to_FF,\
			F_OP(arg1,cl_FF_to_DF(arg2)) );			\
		, /* arg2 DF */						\
			ergebnis_zuweisung CONCAT(NOMAP,r)(DF,		\
			F_OP(arg1,arg2) );				\
		, /* arg2 LF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(LF,cl_LF_to_DF,\
			F_OP(cl_DF_to_LF(arg1,CONCAT(LFlen,s)(arg2)),arg2) ); \
		);							\
	, /* arg1 LF */							\
		floatcase(arg2						\
		, /* arg2 SF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(LF,cl_LF_to_SF,\
			F_OP(arg1,cl_SF_to_LF(arg2,CONCAT(LFlen,s)(arg1))) ); \
		, /* arg2 FF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(LF,cl_LF_to_FF,\
			F_OP(arg1,cl_FF_to_LF(arg2,CONCAT(LFlen,s)(arg1))) ); \
		, /* arg2 DF */						\
			ergebnis_zuweisung CONCAT(MAP,r)(LF,cl_LF_to_DF,\
			F_OP(arg1,cl_DF_to_LF(arg2,CONCAT(LFlen,s)(arg1))) ); \
		, /* arg2 LF */						\
			GEN_LF_OP2_AUX(arg1,arg2,F_OP,r,s,ergebnis_zuweisung) \
		);							\
	);								\
}
#define GEN_LF_OP2_AUX(arg1,arg2,F_OP,r,s,ergebnis_zuweisung)  \
  CONCAT(GEN_LF_OP2_AUX,s)(arg1,arg2,F_OP,r,ergebnis_zuweisung)
#define GEN_LF_OP2_AUX0(arg1,arg2,F_OP,r,ergebnis_zuweisung)  \
  var uintC len1 = TheLfloat(arg1)->len;				\
  var uintC len2 = TheLfloat(arg2)->len;				\
  if (len1 == len2) /* gleich -> direkt ausführen */			\
    { ergebnis_zuweisung CONCAT(NOMAP,r) (LF, F_OP(arg1,arg2)); }	\
  elif (len1 > len2) /* -> arg2 auf die Länge von arg1 bringen */	\
    { ergebnis_zuweisung CONCAT(MAP,r) (LF, LF_shorten_len2,		\
      F_OP(arg1,extend(arg2,len1)) );					\
    }									\
  else /* (len1 < len2) -> arg1 auf die Länge von arg2 bringen */	\
    { ergebnis_zuweisung CONCAT(MAP,r) (LF, LF_shorten_len1,		\
      F_OP(extend(arg1,len2),arg2) );					\
    }
#define LF_shorten_len1(arg)  shorten(arg,len1)
#define LF_shorten_len2(arg)  shorten(arg,len2)
#define GEN_LF_OP2_AUX1(arg1,arg2,F_OP,r,ergebnis_zuweisung)  \
  ergebnis_zuweisung CONCAT(NOMAP,r) (LF, F_OP(arg1,arg2));

#define NOMAP0(F,EXPR)  EXPR
#define NOMAP1(F,EXPR)  EXPR
#define MAP0(F,FN,EXPR)  EXPR
#define MAP1(F,FN,EXPR)  FN(EXPR)

#define LFlen0(arg)  TheLfloat(arg)->len
#define LFlen1(arg)  LF_minlen


// cl_F_extendsqrt(x) erweitert die Genauigkeit eines Floats x um eine Stufe
// SF -> FF -> DF -> LF(4) -> LF(5) -> LF(6) -> ...
// Ein Float mit d Mantissenbits wird so zu einem Float mit
// mindestens d+sqrt(d)+2 Mantissenbits.
extern const cl_F cl_F_extendsqrt (const cl_F& x);

// cl_F_extendsqrtx(x) erweitert die Genauigkeit eines Floats x um eine Stufe
// SF -> FF -> DF -> LF(4) -> LF(5) -> LF(6) -> ...
// Ein Float mit d Mantissenbits und l Exponentenbits wird so zu einem Float
// mit mindestens d+sqrt(d)+2+(l-1) Mantissenbits.
extern const cl_F cl_F_extendsqrtx (const cl_F& x);

// cl_F_shortenrelative(x,y) tries to reduce the size of x, such that one
// wouldn't notice it when adding x to y. y must be /= 0. More precisely,
// this returns a float approximation of x, such that 1 ulp(x) < 1 ulp(y).
extern const cl_F cl_F_shortenrelative (const cl_F& x, const cl_F& y);


// Macro: dispatches according to a float_format_t value.
// floatformatcase(value, SF_statement,FF_statement,DF_statement,LF_statement)
// LF_statement darf auf `len' zugreifen, die zu `value' korrespondierende
// Mantissenlänge (gemessen in Digits).
  #define floatformatcase(value, SF_statement,FF_statement,DF_statement,LF_statement)  \
    { if ((value) <= float_format_sfloat) { SF_statement }		\
      elif ((value) <= float_format_ffloat) { FF_statement }		\
      elif ((value) <= float_format_dfloat) { DF_statement }		\
      else { var uintC len = ceiling((uintC)(value),intDsize); LF_statement } \
    }

}  // namespace cln

#endif /* _CL_F_H */
