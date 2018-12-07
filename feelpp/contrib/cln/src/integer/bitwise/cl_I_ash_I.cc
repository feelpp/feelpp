// ash().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_I ash (const cl_I& x, const cl_I& y)
{
    // Methode:
    // x = 0 -> 0 als Ergebnis
    // y = 0 -> x als Ergebnis
    // y > 0 -> y = intDsize*k + i, j=k+(1 falls i>0, 0 falls i=0).
    //          j Wörter mehr reservieren, k Nullwörter, dann übertragen,
    //          bei i>0: um i Bits links schieben (i=1 geht einfacher).
    // y < 0 -> y <= - intDsize * (Länge(A0) in Digits) -> Ergebnis = 0 oder -1.
    //          Sonst: -y = intDsize*k + i mit k<Länge(A0).
    //                  Übertrage die (Länge(A0)-k) MSDigits,
    //                  falls i>0: schiebe sie um i Bits nach rechts (i=1 geht einfacher).
	if (zerop(x))
		return 0;		// x=0 -> 0 als Ergebnis
	if (zerop(y))
		return x;		// y=0 -> x als Ergebnis
	CL_ALLOCA_STACK;
	if (!minusp(y)) {
		// y>=0
		var uintL i; // i = y mod intDsize, >=0, <intDsize
		var uintC k; // k = y div intDsize, >=0, <2^intCsize
		if (bignump(y)) {
			#if (log2_intDsize+intCsize <= cl_value_len-1)
			// y >= 2^(cl_value_len-1) >= intDsize*2^intCsize
			throw ash_exception(y);
			#else
			// y >= 2^(cl_value_len-1)
			// usable only if y < intDsize*2^intCsize
			var cl_heap_bignum* bn = TheBignum(y);
			var uintC len = bn->length;
			if (len > ceiling(log2_intDsize+intCsize+1,intDsize))
				throw ash_exception(y);
			// bn_minlength <= len <= ceiling(log2_intDsize+intCsize+1,intDsize).
			if (bn_minlength == ceiling(log2_intDsize+intCsize+1,intDsize)
			    || len == ceiling(log2_intDsize+intCsize+1,intDsize))
				if (mspref(arrayMSDptr(bn->data,len),0) >= (uintD)bit((log2_intDsize+intCsize)%intDsize))
					throw ash_exception(y);
			#if (log2_intDsize+intCsize > intDsize)
			#define IF_LENGTH(i)  \
			  if (bn_minlength <= i && i <= ceiling(log2_intDsize+intCsize+1,intDsize) && (i == ceiling(log2_intDsize+intCsize+1,intDsize) || len == i))
			IF_LENGTH(1)
				k = 0;
			else IF_LENGTH(2)
				k = get_uint1D_Dptr(arrayLSDptr(bn->data,2) lspop 1);
			else IF_LENGTH(3)
				k = get_uint2D_Dptr(arrayLSDptr(bn->data,3) lspop 1);
			else IF_LENGTH(4)
				k = get_uint3D_Dptr(arrayLSDptr(bn->data,4) lspop 1);
			else IF_LENGTH(5)
				k = get_uint4D_Dptr(arrayLSDptr(bn->data,5) lspop 1);
			else
				throw runtime_exception();
			#undef IF_LENGTH
			k = k << (intDsize-log2_intDsize);
			#else
			// log2_intDsize+intCsize <= intDsize,
			// implies len==1 or len==2 && lspref(arrayLSDptr(bn->data,len),1) == 0.
			k = 0;
			#endif
			k |= lspref(arrayLSDptr(bn->data,len),0) >> log2_intDsize;
			i = lspref(arrayLSDptr(bn->data,len),0) % intDsize;
			#endif
		} else {
			var uintV y_ = FN_to_V(y); // Wert von y, >=0, <intDsize*2^intCsize
			i = y_%intDsize;
			k = floor(y_,intDsize);
		}
		var uintD* LSDptr;
		var uintC len;
		var const uintD* x_LSDptr;
		I_to_NDS_nocopy(x, ,len=,x_LSDptr=,false,); // DS zu x bilden.
		if (k >= (uintC)(~len)) // kann len+k+1 Überlauf geben?
			{ throw ash_exception(y); } // ja -> Fehler
		num_stack_alloc_1(len+k,,LSDptr=);
		LSDptr = clear_loop_lsp(LSDptr,k); // k Nulldigits
		var uintD* MSDptr = copy_loop_lsp(x_LSDptr,LSDptr,len);
		// Nun ist MSDptr/len/LSDptr die DS zu x.
		// Oberhalb von ihr liegen k Nulldigits, unterhalb ist 1 Digit Platz.
		// MSDptr/len+k/.. ist jetzt die Gesamt-DS.
		// Noch um i Bits nach links schieben:
		if (!(i==0)) // Bei i>0
		  { // noch ein weiteres Digit dazunehmen (Vorzeichen)
		    {var uintD sign = sign_of_sintD(mspref(MSDptr,0));
		     lsprefnext(MSDptr) = sign;
		     len++;
		    }
		    // Schiebeschleife: die unteren len Digits um i Bits schieben
		    if (i==1)
		      { shift1left_loop_lsp(LSDptr,len); }
		    else
		      { shiftleft_loop_lsp(LSDptr,len,i,0); }
		  }
		return DS_to_I(MSDptr,len+k);
	} else {
		// y<0
		var uintL i; // i = (-y) mod intDsize, >=0, <intDsize
		var uintC k; // k = (-y) div intDsize, >=0, <2^intCsize
		if (bignump(y)) {
			#if (log2_intDsize+intCsize <= cl_value_len-1)
			// -y-1 >= 2^(cl_value_len-1) >= intDsize*2^intCsize
			goto sign;
			#else
			// -y-1 >= 2^(cl_value_len-1)
			// usable only if -y-1 < intDsize*2^intCsize
			// We write -y-1 = lognot(y) = k*intDsize+i and then add 1.
			var cl_heap_bignum* bn = TheBignum(y);
			var uintC len = bn->length;
			if (len > ceiling(log2_intDsize+intCsize+1,intDsize))
				goto sign;
			// bn_minlength <= len <= ceiling(log2_intDsize+intCsize+1,intDsize).
			if (bn_minlength == ceiling(log2_intDsize+intCsize+1,intDsize)
			    || len == ceiling(log2_intDsize+intCsize+1,intDsize))
				if (mspref(arrayMSDptr(bn->data,len),0) < (uintD)(-bit((log2_intDsize+intCsize)%intDsize)))
					goto sign;
			#if (log2_intDsize+intCsize > intDsize)
			#define IF_LENGTH(i)  \
			  if (bn_minlength <= i && i <= ceiling(log2_intDsize+intCsize+1,intDsize) && (i == ceiling(log2_intDsize+intCsize+1,intDsize) || len == i))
			IF_LENGTH(1)
				k = 0;
			else IF_LENGTH(2)
				k = ~get_sint1D_Dptr(arrayLSDptr(bn->data,2) lspop 1);
			else IF_LENGTH(3)
				k = ~get_sint2D_Dptr(arrayLSDptr(bn->data,3) lspop 1);
			else IF_LENGTH(4)
				k = ~get_sint3D_Dptr(arrayLSDptr(bn->data,4) lspop 1);
			else IF_LENGTH(5)
				k = ~get_sint4D_Dptr(arrayLSDptr(bn->data,5) lspop 1);
			else
				throw runtime_exception();
			#undef IF_LENGTH
			k = k << (intDsize-log2_intDsize);
			#else
			// log2_intDsize+intCsize <= intDsize,
			// implies len==1 or len==2 && lspref(arrayLSDptr(bn->data,len),1) == ~0.
			k = 0;
			#endif
			k |= (uintD)(~lspref(arrayLSDptr(bn->data,len),0)) >> log2_intDsize;
			i = (uintD)(-lspref(arrayLSDptr(bn->data,len),0)) % intDsize;
			if (i == 0)
				if (++k == 0)
					goto sign;
			#endif
		} else {
			var uintV y_ = -FN_to_V(y); // Wert von -y, >0, <intDsize*2^intCsize
			i = y_%intDsize;
			k = floor(y_,intDsize);
		}
		// DS zu x bilden:
		var uintD* MSDptr;
		var uintC len;
		I_to_NDS(x, MSDptr=,len=,); // DS zu x bilden.
		if (k>=len) goto sign; // -y >= intDsize*len -> Vorzeichen von x zurück
		len -= k; // rechte k Digits einfach streichen
		// Noch ist len>0. Um i Bits nach rechts schieben:
		if (!(i==0)) // Bei i>0:
		  { // Schiebe len Digits ab MSDptr um i Bits nach rechts:
		    if (i==1)
		      { shift1right_loop_msp(MSDptr,len,sign_of_sintD(mspref(MSDptr,0))); }
		      else
		      { shiftrightsigned_loop_msp(MSDptr,len,i); }
		  }
		return DS_to_I(MSDptr,len);
	}
sign:	// Ergebnis ist 0, falls x>=0, und -1, falls x<0:
	return (minusp(x) ? cl_I(-1) : cl_I(0));
}

}  // namespace cln
