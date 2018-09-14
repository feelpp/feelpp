// cl_LF implementation

#ifndef _CL_LF_IMPL_H
#define _CL_LF_IMPL_H

#include "cln/number.h"
#include "float/lfloat/cl_LF.h"
#include "cln/malloc.h"
#include "base/cl_offsetof.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

extern cl_class cl_class_lfloat;

// Builds a long-float, without filling the mantissa.
// allocate_lfloat(len,expo,sign)
// > uintC len: length of mantissa (in digits)
// > uintE expo: exponent
// > cl_signean sign: sign (0 = +, -1 = -)
// The long-float is only complete when the mantissa has been filled in!
inline cl_heap_lfloat* allocate_lfloat (uintC len, uintE expo, cl_signean sign)
{
	cl_heap_lfloat* p = (cl_heap_lfloat*) malloc_hook(offsetofa(cl_heap_lfloat,data)+sizeof(uintD)*len);
	p->refcount = 1;
	p->type = &cl_class_lfloat;
	p->len = len;
	p->sign = sign;
	p->expo = expo;
	return p;
}

// Private constructor.
// ptr should be the result of some allocate_lfloat() call.
inline cl_LF::cl_LF (cl_heap_lfloat* ptr) : cl_F ((cl_private_thing) ptr) {}

// Both work, but the first definition results in less compiler-generated
// temporaries.
#if 1
  #define Lfloat  cl_heap_lfloat*
#else
  #define Lfloat  cl_LF
#endif

// Pointers to the mantissa.
#if 1
  inline const uintD* LF_MSDptr (Lfloat lf)
    { return (const uintD*) arrayMSDptr(lf->data,lf->len); }
  inline const uintD* LF_LSDptr (Lfloat lf)
    { return (const uintD*) arrayLSDptr(lf->data,lf->len); }
#endif
  inline const uintD* LF_MSDptr (const cl_LF& obj)
    { var cl_heap_lfloat* lf = TheLfloat(obj); return (const uintD*) arrayMSDptr(lf->data,lf->len); }
  inline const uintD* LF_LSDptr (const cl_LF& obj)
    { var cl_heap_lfloat* lf = TheLfloat(obj); return (const uintD*) arrayLSDptr(lf->data,lf->len); }


// Entpacken eines Long-Float:
// LF_decode(obj, zero_statement, sign=,exp=,mantMSDptr=,mantlen=,mantLSDptr=);
// zerlegt ein Long-Float obj.
// Ist obj=0.0, wird zero_statement ausgeführt.
// Sonst: cl_signean sign = Vorzeichen (0 = +, -1 = -),
//        sintE exp = Exponent (vorzeichenbehaftet),
//        UDS mantMSDptr/mantlen/mantLSDptr = Mantisse
//          (>= 2^(intDsize*mantlen-1), < 2^(intDsize*mantlen)),
//          mit mantlen>=LF_minlen.
  #define LF_decode(obj, zero_statement, sign_zuweisung,exp_zuweisung,mantMSDptr_zuweisung,mantlen_zuweisung,mantLSDptr_zuweisung)  \
    { var Lfloat _x = TheLfloat(obj);					\
      var uintE uexp = _x->expo;					\
      if (uexp==0)							\
        { unused (mantlen_zuweisung _x->len); zero_statement } /* e=0 -> Zahl 0.0 */\
        else								\
        { exp_zuweisung (sintE)(uexp - LF_exp_mid);	/* Exponent */	\
          sign_zuweisung _x->sign;			/* Vorzeichen */\
          unused (mantMSDptr_zuweisung arrayMSDptr(_x->data, (uintP)(mantlen_zuweisung _x->len))); /* Mantissen-UDS */\
          unused (mantLSDptr_zuweisung arrayLSDptr(_x->data, (uintP)(mantlen_zuweisung _x->len))); \
    }   }

// Einpacken eines Long-Float:
// encode_LF0(len) liefert ein Long-Float 0.0 mit len Digits.
// > uintC len: Anzahl der Digits
// < cl_LF ergebnis: neues Long-Float 0.0 mit len Digits
inline const cl_LF encode_LF0 (uintC len)
{
	var Lfloat erg = allocate_lfloat(len,0,0); // Exponent 0, Vorzeichen +
	DS_clear_loop(arrayMSDptr(TheLfloat(erg)->data,len),len,arrayLSDptr(TheLfloat(erg)->data,len)); // Mantisse := 0
	return erg;
}

// Einpacken eines Long-Float:
// encode_LF1s(sign,len) liefert ein Long-Float +-1.0 mit len Digits.
// > cl_signean sign: Vorzeichen
// > uintC len: Anzahl der Digits
// < cl_LF ergebnis: neues Long-Float +1.0 oder -1.0 mit len Digits
inline const cl_LF encode_LF1s (cl_signean sign, uintC len)
{
	var Lfloat erg = allocate_lfloat(len,LF_exp_mid+1,sign); // Exponent 1
	mspref(arrayMSDptr(TheLfloat(erg)->data,len),0) = bit(intDsize-1); // Mantisse := 2^(intDsize*len-1)
	DS_clear_loop(arrayMSDptr(TheLfloat(erg)->data,len) mspop 1,len-1,arrayLSDptr(TheLfloat(erg)->data,len));
	return erg;
}

// Einpacken eines Long-Float:
// encode_LF1(len) liefert ein Long-Float 1.0 mit len Digits.
// > uintC len: Anzahl der Digits
// < cl_LF ergebnis: neues Long-Float 1.0 mit len Digits
inline const cl_LF encode_LF1 (uintC len)
{
	return encode_LF1s(0,len);
}

// Einpacken eines Long-Float:
// encode_LFu(sign,uexp,mantMSDptr,mantlen) liefert ein Long-Float
// > cl_signean sign: Vorzeichen
// > uintE exp: Exponent + LF_exp_mid
// > uintD* mantMSDptr: Pointer auf eine NUDS mit gesetztem höchstem Bit
// > uintC mantlen: Anzahl der Digits, >= LF_minlen
// < cl_LF erg: neues Long-Float mit der UDS mantMSDptr/mantlen/.. als Mantisse
// Der Exponent wird nicht auf Überlauf/Unterlauf getestet.
inline const cl_LF encode_LFu (cl_signean sign, uintE uexp, const uintD* mantMSDptr, uintC mantlen)
{
	var Lfloat erg = allocate_lfloat(mantlen,uexp,sign); /* Exponent */
	copy_loop_msp(mantMSDptr,arrayMSDptr(TheLfloat(erg)->data,mantlen),mantlen); /* Mantisse übertragen */
	return erg;
}

// Einpacken eines Long-Float:
// encode_LF(sign,exp,mantMSDptr,mantlen) liefert ein Long-Float
// > cl_signean sign: Vorzeichen
// > sintE exp: Exponent
// > uintD* mantMSDptr: Pointer auf eine NUDS mit gesetztem höchstem Bit
// > uintC mantlen: Anzahl der Digits, >= LF_minlen
// < cl_LF erg: neues Long-Float mit der UDS mantMSDptr/mantlen/.. als Mantisse
// Der Exponent wird nicht auf Überlauf/Unterlauf getestet.
inline const cl_LF encode_LF (cl_signean sign, sintE exp, const uintD* mantMSDptr, uintC mantlen)
{
	return encode_LFu(sign,LF_exp_mid+(uintE)exp,mantMSDptr,mantlen);
}

// Einpacken eines Long-Float:
// encode_LF_array(sign,exp,mantarr,mantlen) liefert ein Long-Float
// > cl_signean sign: Vorzeichen
// > sintE exp: Exponent
// > uintD mantarr[]: NUDS mit gesetztem höchstem Bit
// > uintC mantlen: Anzahl der Digits, >= LF_minlen
// < cl_LF erg: neues Long-Float mit der UDS mantarr[] als Mantisse
// Der Exponent wird nicht auf Überlauf/Unterlauf getestet.
#define encode_LF_array(sign,exp,mantarr,mantlen)  \
  encode_LF(sign,exp,arrayMSDptr(mantarr,mantlen),mantlen)

}  // namespace cln

#endif /* _CL_LF_IMPL_H */
