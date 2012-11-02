// cl_R internals

#ifndef _CL_R_H
#define _CL_R_H

#include "cln/number.h"
#include "cln/real.h"

namespace cln {

extern cl_class cl_class_bignum;
extern cl_class cl_class_ratio;
extern cl_class cl_class_ffloat;
extern cl_class cl_class_dfloat;
extern cl_class cl_class_lfloat;

// Type tests.
inline bool rationalp (const cl_R& x)
{
	if (!x.pointer_p()) {
		if (x.nonpointer_tag() == cl_FN_tag)
			return true;
	} else {
		if (x.pointer_type()->flags & cl_class_flags_subclass_rational)
			return true;
	}
	return false;
}
inline bool integerp (const cl_R& x)
{
	if (!x.pointer_p()) {
		if (x.nonpointer_tag() == cl_FN_tag)
			return true;
	} else {
		if (x.pointer_type() == &cl_class_bignum)
			return true;
	}
	return false;
}
inline bool floatp (const cl_R& x)
{
	if (!x.pointer_p()) {
		switch (x.nonpointer_tag()) {
		case cl_SF_tag:
		#if defined(CL_WIDE_POINTERS)
		case cl_FF_tag:
		#endif
			return true;
		}
	} else {
		if (x.pointer_type()->flags & cl_class_flags_subclass_float)
			return true;
	}
	return false;
}

// Comparison with a fixnum.
inline bool eq (const cl_R& x, sint32 y)
{
	return x.word == cl_combine(cl_FN_tag,y);
}

inline bool exact_zerop (const cl_R& x)
{
	return eq(x,0);
}


// Macro: verteilt je nach Real-Typ eines Floats x auf 2 Statements,
// die x vom jeweiligen Real-Typ benutzen dürfen.
// realcase2(x, RA_statement,F_statement);
// x sollte eine Variable sein.
  #define realcase2(x, RA_statement,F_statement) \
    if (rationalp(x))							     \
      { var cl_RA& __tmp = *(cl_RA*)&x; var cl_RA& x = __tmp; RA_statement } \
    else								     \
      { var cl_F& __tmp = *(cl_F*)&x; var cl_F& x = __tmp; F_statement }

// Macro: verteilt je nach Real-Typ eines Floats x auf 7 Statements.
// realtypecase(x, FN_statement,BN_statement,RT_statement,SF_statement,FF_statement,DF_statement,LF_statement);
// x sollte eine Variable sein.
#ifdef CL_WIDE_POINTERS
  #define realtypecase(x, FN_statement,BN_statement,RT_statement,SF_statement,FF_statement,DF_statement,LF_statement) \
    if (!(x).pointer_p())						\
      switch ((x).nonpointer_tag())					\
        { case cl_FN_tag: { FN_statement } break;			\
          case cl_SF_tag: { SF_statement } break;			\
          case cl_FF_tag: { FF_statement } break;			\
          default: NOTREACHED						\
        }								\
      else {								\
        if ((x).pointer_type() == &cl_class_bignum) { BN_statement }	\
        else if ((x).pointer_type() == &cl_class_ratio) { RT_statement } \
        else if ((x).pointer_type() == &cl_class_dfloat) { DF_statement } \
        else if ((x).pointer_type() == &cl_class_lfloat) { LF_statement } \
        else NOTREACHED							\
      }
#else
  #define realtypecase(x, FN_statement,BN_statement,RT_statement,SF_statement,FF_statement,DF_statement,LF_statement) \
    if (!(x).pointer_p())						\
      switch ((x).nonpointer_tag())					\
        { case cl_FN_tag: { FN_statement } break;			\
          case cl_SF_tag: { SF_statement } break;			\
          default: NOTREACHED						\
        }								\
      else {								\
        if ((x).pointer_type() == &cl_class_bignum) { BN_statement }	\
        else if ((x).pointer_type() == &cl_class_ratio) { RT_statement } \
        else if ((x).pointer_type() == &cl_class_ffloat) { FF_statement } \
        else if ((x).pointer_type() == &cl_class_dfloat) { DF_statement } \
        else if ((x).pointer_type() == &cl_class_lfloat) { LF_statement } \
        else NOTREACHED							\
      }
#endif

// Macro: verteilt je nach Real-Typ eines Floats x auf 7 Statements,
// die x vom jeweiligen Real-Typ benutzen dürfen.
// realcase7(x, FN_statement,BN_statement,RT_statement,SF_statement,FF_statement,DF_statement,LF_statement);
// x sollte eine Variable sein.
  #define realcase7(x, FN_statement,BN_statement,RT_statement,SF_statement,FF_statement,DF_statement,LF_statement) \
    realtypecase(x							   \
      , var cl_FN& __tmp = *(cl_FN*)&x; var cl_FN& x = __tmp; FN_statement \
      , var cl_BN& __tmp = *(cl_BN*)&x; var cl_BN& x = __tmp; BN_statement \
      , var cl_RT& __tmp = *(cl_RT*)&x; var cl_RT& x = __tmp; RT_statement \
      , var cl_SF& __tmp = *(cl_SF*)&x; var cl_SF& x = __tmp; SF_statement \
      , var cl_FF& __tmp = *(cl_FF*)&x; var cl_FF& x = __tmp; FF_statement \
      , var cl_DF& __tmp = *(cl_DF*)&x; var cl_DF& x = __tmp; DF_statement \
      , var cl_LF& __tmp = *(cl_LF*)&x; var cl_LF& x = __tmp; LF_statement \
      )

// Macro: verteilt je nach Real-Typ eines Floats x auf 6 Statements,
// die x vom jeweiligen Real-Typ benutzen dürfen.
// realcase6(x, I_statement,RT_statement,SF_statement,FF_statement,DF_statement,LF_statement);
// x sollte eine Variable sein.
  #define realcase6(x, I_statement,RT_statement,SF_statement,FF_statement,DF_statement,LF_statement) \
    realcase7(x, I_statement,I_statement,RT_statement,SF_statement,FF_statement,DF_statement,LF_statement)

// contagion(x,y) liefert eine reelle Zahl, die so ungenau ist wie die
// ungenauere der beiden reellen Zahlen x und y.
extern const cl_R contagion (const cl_R& x, const cl_R& y);
// ?? Lieber ein uintL (0, SF_mant_len+1, FF_mant_len+1, DF_mant_len+1, intDsize*len liefern, weniger Kopieraufwand!

// GEN_R_OP1_2(arg1,R_OP,ergebnis_zuweisung)
// generates the body of a real operation with one argument.
// Distinguish two cases (rational/float) only.
#define GEN_R_OP1_2(arg1,R_OP,ergebnis_zuweisung)  \
{									\
	realcase2(arg1							\
	, ergebnis_zuweisung R_OP(arg1);				\
	, ergebnis_zuweisung R_OP(arg1);				\
	);								\
}

// GEN_R_OP1_7(arg1,R_OP,ergebnis_zuweisung)
// generates the body of a real operation with one argument.
// Full type dispatch, faster than GEN_R_OP1_2.
#define GEN_R_OP1_7(arg1,R_OP,ergebnis_zuweisung)  \
{									\
	realcase7(arg1							\
	, ergebnis_zuweisung R_OP(arg1);				\
	, ergebnis_zuweisung R_OP(arg1);				\
	, ergebnis_zuweisung R_OP(arg1);				\
	, ergebnis_zuweisung R_OP(arg1);				\
	, ergebnis_zuweisung R_OP(arg1);				\
	, ergebnis_zuweisung R_OP(arg1);				\
	, ergebnis_zuweisung R_OP(arg1);				\
	);								\
}

// GEN_R_OP2_2(arg1,arg2,R_OP,ergebnis_zuweisung)
// generates the body of a real operation with two arguments.
// Distinguish two cases (rational/float) only.
#define GEN_R_OP2_2(arg1,arg2,R_OP,ergebnis_zuweisung)  \
{									\
	realcase2(arg1							\
	,	realcase2(arg2						\
		,	/* beides rationale Zahlen */			\
			ergebnis_zuweisung R_OP(arg1,arg2);		\
		,	/* arg1 rational, arg2 Float -> arg1 in Float umwandeln */ \
			ergebnis_zuweisung R_OP(cl_float(arg1,arg2),arg2); \
		);							\
	,	realcase2(arg2						\
		,	/* arg1 Float, arg2 rational -> arg2 in Float umwandeln */ \
			ergebnis_zuweisung R_OP(arg1,cl_float(arg2,arg1)); \
		,	/* beides Floats */				\
			ergebnis_zuweisung R_OP(arg1,arg2);		\
		);							\
	);								\
}

// cl_somefloat(x,y) wandelt eine reelle Zahl x in ein Float-Format um
// (das von y, falls x rational ist) und rundet dabei nötigenfalls.
// > x: eine reelle Zahl
// > y: ein Float
// < ergebnis: x als Float
inline const cl_F cl_somefloat (const cl_R& x, const cl_F& y)
{
	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		return cl_float(x,y);
	} else {
		DeclareType(cl_F,x);
		return x;
	}
}

}  // namespace cln

#endif /* _CL_R_H */
