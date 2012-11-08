// Fast integer multiplication using Nussbaumer's FFT based algorithm.
// [Donald Ervin Knuth: The Art of Computer Programming, Vol. II:
//  Seminumerical Algorithms, second edition.
//  Section 4.6.4, exercise 59, p. 503, 652-654.]
// [Henri Jean Nussbaumer, IEEE Trans. ASSP-28 (1980), 205-215.]
// Bruno Haible 4.-5.5.1996

// This algorithm has the benefit of working on entire words, not single bits,
// and involving no non-integer numbers. (The root of unity is chosen in
// an appropriate polynomial ring.)

// If at the beginning all words x_i, y_i are >= 0 and < M, then
// the intermediate X_{i,j}, Y_{i,j} are < M * N in absolute value
// (where N = number of words), hence the |Z_{i,j}| < M^2 * N^2.
// We therefore reserve 2 32-bit words for every X_{i,j} and 4 32-bit words
// for every Z_{i,j}.

#if !(intDsize==32)
#error "nussbaumer implemented only for intDsize==32"
#endif

// Define this if you want the external loops instead of inline operations.
//#define NUSS_IN_EXTERNAL_LOOPS
#define NUSS_OUT_EXTERNAL_LOOPS

// Define this if you want inline operations which access the stack directly.
// This looks like better code, but is in effect 3% slower. No idea why.
//#define NUSS_ASM_DIRECT

// Define this for (cheap) consistency checks.
//#define DEBUG_NUSS

// Define this for extensive consistency checks.
//#define DEBUG_NUSS_OPERATIONS

#if (intDsize==32)

//typedef struct { sint32 iw1; uint32 iw0; } nuss_inword;
//typedef struct { uint32 iw0; sint32 iw1; } nuss_inword;
typedef struct { uintD _iw[2]; } nuss_inword;
#if CL_DS_BIG_ENDIAN_P
  #define iw1 _iw[0]
  #define iw0 _iw[1]
#else
  #define iw0 _iw[0]
  #define iw1 _iw[1]
#endif

//typedef struct { sint32 ow3; uint32 ow2; uint32 ow1; uint32 ow0; } nuss_outword;
//typedef struct { uint32 ow0; uint32 ow1; uint32 ow2; sint32 ow3; } nuss_outword;
typedef struct { uintD _ow[4]; } nuss_outword;
#if CL_DS_BIG_ENDIAN_P
  #define ow3 _ow[0]
  #define ow2 _ow[1]
  #define ow1 _ow[2]
  #define ow0 _ow[3]
#else
  #define ow0 _ow[0]
  #define ow1 _ow[1]
  #define ow2 _ow[2]
  #define ow3 _ow[3]
#endif

// r := a + b
static inline void add (const nuss_inword& a, const nuss_inword& b, nuss_inword& r)
{
#if defined(__GNUC__) && defined(__i386__)
	var uintD dummy;
  #ifdef NUSS_ASM_DIRECT
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"addl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.iw0), "m" (b.iw0), "m" (r.iw0)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"adcl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.iw1), "m" (b.iw1), "m" (r.iw1)
		: "cc"
		);
  #else
    #if CL_DS_BIG_ENDIAN_P
	__asm__ __volatile__ (
		"movl 4(%1),%0" "\n\t"
		"addl 4(%2),%0" "\n\t"
		"movl %0,4(%3)" "\n\t"
		"movl (%1),%0"  "\n\t"
		"adcl (%2),%0"  "\n\t"
		"movl %0,(%3)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b), "r" (&r)
		: "cc"
		);
    #else
	__asm__ __volatile__ (
		"movl (%1),%0"  "\n\t"
		"addl (%2),%0"  "\n\t"
		"movl %0,(%3)"  "\n\t"
		"movl 4(%1),%0" "\n\t"
		"adcl 4(%2),%0" "\n\t"
		"movl %0,4(%3)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b), "r" (&r)
		: "cc"
		);
    #endif
  #endif
#elif defined(NUSS_IN_EXTERNAL_LOOPS)
	add_loop_lsp(arrayLSDptr(a._iw,2),arrayLSDptr(b._iw,2),arrayLSDptr(r._iw,2),2);
#else
	var uint32 tmp;

	tmp = a.iw0 + b.iw0;
	if (tmp >= a.iw0) {
		// no carry
		r.iw0 = tmp;
		r.iw1 = a.iw1 + b.iw1;
	} else {
		// carry
		r.iw0 = tmp;
		r.iw1 = a.iw1 + b.iw1 + 1;
	}
#endif
}

// r := a - b
static inline void sub (const nuss_inword& a, const nuss_inword& b, nuss_inword& r)
{
#if defined(__GNUC__) && defined(__i386__)
	var uintD dummy;
  #ifdef NUSS_ASM_DIRECT
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"subl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.iw0), "m" (b.iw0), "m" (r.iw0)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"sbbl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.iw1), "m" (b.iw1), "m" (r.iw1)
		: "cc"
		);
  #else
    #if CL_DS_BIG_ENDIAN_P
	__asm__ __volatile__ (
		"movl 4(%1),%0" "\n\t"
		"subl 4(%2),%0" "\n\t"
		"movl %0,4(%3)" "\n\t"
		"movl (%1),%0"  "\n\t"
		"sbbl (%2),%0"  "\n\t"
		"movl %0,(%3)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b), "r" (&r)
		: "cc"
		);
    #else
	__asm__ __volatile__ (
		"movl (%1),%0"  "\n\t"
		"subl (%2),%0"  "\n\t"
		"movl %0,(%3)"  "\n\t"
		"movl 4(%1),%0" "\n\t"
		"sbbl 4(%2),%0" "\n\t"
		"movl %0,4(%3)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b), "r" (&r)
		: "cc"
		);
    #endif
  #endif
#elif defined(NUSS_IN_EXTERNAL_LOOPS)
	sub_loop_lsp(arrayLSDptr(a._iw,2),arrayLSDptr(b._iw,2),arrayLSDptr(r._iw,2),2);
#else
	var uint32 tmp;

	tmp = a.iw0 - b.iw0;
	if (tmp <= a.iw0) {
		// no carry
		r.iw0 = tmp;
		r.iw1 = a.iw1 - b.iw1;
	} else {
		// carry
		r.iw0 = tmp;
		r.iw1 = a.iw1 - b.iw1 - 1;
	}
#endif
}

// r := a * b
static void mul (const nuss_inword& a, const nuss_inword& b, nuss_outword& r)
{
#ifdef NUSS_IN_EXTERNAL_LOOPS
	mulu_2loop(arrayLSDptr(a._iw,2),2, arrayLSDptr(b._iw,2),2, arrayLSDptr(r._ow,4));
	if ((sintD)mspref(arrayMSDptr(a._iw,2),0) < 0)
		subfrom_loop_lsp(arrayLSDptr(b._iw,2),arrayLSDptr(r._ow,4) lspop 2,2);
	if ((sintD)mspref(arrayMSDptr(b._iw,2),0) < 0)
		subfrom_loop_lsp(arrayLSDptr(a._iw,2),arrayLSDptr(r._ow,4) lspop 2,2);
#else
	if (a.iw1 == 0) {
		// a small positive
		if (b.iw1 == 0) {
			// a, b small positive
			mulu32(a.iw0, b.iw0, r.ow1 =, r.ow0 =);
			r.ow3 = 0; r.ow2 = 0;
			return;
		}
		else if (b.iw1 == -(uint32)1 && b.iw0 != 0) {
			// b small negative
			var uint32 hi, lo;
			mulu32(a.iw0, -b.iw0, hi=, lo=);
			r.ow0 = -lo;
			if (lo) {
				r.ow1 = ~hi;
			} else if (hi) {
				r.ow1 = -hi;
			} else /* a.iw0 == 0 */ {
				r.ow3 = 0; r.ow2 = 0; r.ow1 = 0;
				return;
			}
			r.ow3 = -(uint32)1; r.ow2 = -(uint32)1;
			return;
		}
		var uint32 hi1, lo1, hi0;
		mulu32(a.iw0, b.iw0, hi0 =, r.ow0 =);
		mulu32(a.iw0, b.iw1, hi1 =, lo1 =);
		if ((lo1 += hi0) < hi0)
			hi1++;
		// hi1|lo1|r.ow0 = a.iw0 * b(unsigned).
		r.ow1 = lo1;
		if ((sint32)b.iw1 >= 0) {
			r.ow2 = hi1;
			r.ow3 = 0;
		} else {
			// b was negative -> subtract a * 2^64
			if (a.iw0) {
				r.ow2 = hi1 - a.iw0;
				r.ow3 = -(uint32)1;
			} else /* a.iw0 == 0 */ {
				r.ow3 = 0; r.ow2 = 0;
			}
		}
		return;
	}
	else if (a.iw1 == -(uint32)1 && a.iw0 != 0) {
		// a small negative
		if (b.iw1 == 0) {
			// b small positive
			var uint32 hi, lo;
			mulu32(-a.iw0, b.iw0, hi=, lo=);
			r.ow0 = -lo;
			if (lo) {
				r.ow1 = ~hi;
			} else if (hi) {
				r.ow1 = -hi;
			} else /* b.iw0 == 0 */ {
				r.ow3 = 0; r.ow2 = 0; r.ow1 = 0;
				return;
			}
			r.ow3 = -(uint32)1; r.ow2 = -(uint32)1;
			return;
		}
		else if (b.iw1 == -(uint32)1 && b.iw0 != 0) {
			// a, b small negative
			mulu32(-a.iw0, -b.iw0, r.ow1 =, r.ow0 =);
			r.ow3 = 0; r.ow2 = 0;
			return;
		}
		var uint32 hi1, lo1, hi0, lo0;
		mulu32(-a.iw0, b.iw0, hi0 =, lo0 =);
		mulu32(-a.iw0, b.iw1, hi1 =, lo1 =);
		if ((lo1 += hi0) < hi0)
			hi1++;
		// hi1|lo1|lo0 = -a * b(unsigned).
		if (lo0) {
			lo0 = -lo0;
			lo1 = ~lo1;
			hi1 = ~hi1;
		} else if (lo1) {
			lo1 = -lo1;
			hi1 = ~hi1;
		} else
			hi1 = -hi1;
		// hi1|lo1|lo0 = a * b(unsigned).
		r.ow0 = lo0;
		r.ow1 = lo1;
		if ((sint32)b.iw1 >= 0) {
			r.ow2 = hi1;
			r.ow3 = -(uint32)1;
		} else {
			// b was negative -> subtract a * 2^64
			r.ow2 = hi1 - a.iw0;
			r.ow3 = 0;
		}
		return;
	}
	else if (b.iw1 == 0) {
		// b small positive
		var uint32 hi1, lo1, hi0;
		mulu32(b.iw0, a.iw0, hi0 =, r.ow0 =);
		mulu32(b.iw0, a.iw1, hi1 =, lo1 =);
		if ((lo1 += hi0) < hi0)
			hi1++;
		// hi1|lo1|r.ow0 = a(unsigned) * b.iw0.
		r.ow1 = lo1;
		if ((sint32)a.iw1 >= 0) {
			r.ow2 = hi1;
			r.ow3 = 0;
		} else {
			// a was negative -> subtract b * 2^64
			if (b.iw0) {
				r.ow2 = hi1 - b.iw0;
				r.ow3 = -(uint32)1;
			} else /* b.iw0 == 0 */ {
				r.ow3 = 0; r.ow2 = 0;
			}
		}
		return;
	}
	else if (b.iw1 == -(uint32)1 && b.iw0 != 0) {
		// b small negative
		var uint32 hi1, lo1, hi0, lo0;
		mulu32(-b.iw0, a.iw0, hi0 =, lo0 =);
		mulu32(-b.iw0, a.iw1, hi1 =, lo1 =);
		if ((lo1 += hi0) < hi0)
			hi1++;
		// hi1|lo1|lo0 = a(unsigned) * -b.
		if (lo0) {
			lo0 = -lo0;
			lo1 = ~lo1;
			hi1 = ~hi1;
		} else if (lo1) {
			lo1 = -lo1;
			hi1 = ~hi1;
		} else
			hi1 = -hi1;
		// hi1|lo1|lo0 = a(unsigned) * b.
		r.ow0 = lo0;
		r.ow1 = lo1;
		if ((sint32)a.iw1 >= 0) {
			r.ow2 = hi1;
			r.ow3 = -(uint32)1;
		} else {
			// a was negative -> subtract b * 2^64
			r.ow2 = hi1 - b.iw0;
			r.ow3 = 0;
		}
		return;
	}
	// This is the main and most frequent case (65% to 80%).
	var uint32 w3, w2, w1, hi, lo;
	mulu32(a.iw0, b.iw0, w1=, r.ow0=);
	mulu32(a.iw1, b.iw1, w3=, w2=);
	mulu32(a.iw0, b.iw1, hi=, lo=);
	if ((w1 += lo) < lo)
		hi++;
	if ((w2 += hi) < hi)
		w3++;
	mulu32(a.iw1, b.iw0, hi=, lo=);
	if ((w1 += lo) < lo)
		hi++;
	if ((w2 += hi) < hi)
		w3++;
	// w3|w2|w1|r.ow0 = a(unsigned) * b(unsigned).
	r.ow1 = w1;
	if ((sint32)a.iw1 < 0) {
		// a was negative -> subtract b * 2^64
		if (w2 >= b.iw0) {
			w2 -= b.iw0;
			w3 -= b.iw1;
		} else {
			// carry
			w2 -= b.iw0;
			w3 = w3 - b.iw1 - 1;
		}
	}
	if ((sint32)b.iw1 < 0) {
		// b was negative -> subtract a * 2^64
		if (w2 >= a.iw0) {
			w2 -= a.iw0;
			w3 -= a.iw1;
		} else {
			// carry
			w2 -= a.iw0;
			w3 = w3 - a.iw1 - 1;
		}
	}
	r.ow2 = w2;
	r.ow3 = w3;
	return;
#endif
}
#ifdef DEBUG_NUSS_OPERATIONS
static void mul_doublecheck (const nuss_inword& a, const nuss_inword& b, nuss_outword& r)
{
	nuss_outword or;
	mulu_2loop(arrayLSDptr(a._iw,2),2, arrayLSDptr(b._iw,2),2, arrayLSDptr(or._ow,4));
	if ((sintD)mspref(arrayMSDptr(a._iw,2),0) < 0)
		subfrom_loop_lsp(arrayLSDptr(b._iw,2),arrayLSDptr(or._ow,4) lspop 2,2);
	if ((sintD)mspref(arrayMSDptr(b._iw,2),0) < 0)
		subfrom_loop_lsp(arrayLSDptr(a._iw,2),arrayLSDptr(or._ow,4) lspop 2,2);
	mul(a,b, r);
	if (compare_loop_msp(arrayMSDptr(r._ow,4),arrayMSDptr(or._ow,4),4))
		throw runtime_exception();
}
#define mul mul_doublecheck
#endif

// r := 0
static inline void zero (nuss_outword& r)
{
	r.ow0 = 0;
	r.ow1 = 0;
	r.ow2 = 0;
	r.ow3 = 0;
}

// r := a + b
static inline void add (const nuss_outword& a, const nuss_outword& b, nuss_outword& r)
{
#if defined(__GNUC__) && defined(__i386__)
	var uintD dummy;
  #ifdef NUSS_ASM_DIRECT
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"addl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.ow0), "m" (b.ow0), "m" (r.ow0)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"adcl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.ow1), "m" (b.ow1), "m" (r.ow1)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"adcl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.ow2), "m" (b.ow2), "m" (r.ow2)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"adcl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.ow3), "m" (b.ow3), "m" (r.ow3)
		: "cc"
		);
  #else
    #if CL_DS_BIG_ENDIAN_P
	__asm__ __volatile__ (
		"movl 12(%1),%0" "\n\t"
		"addl 12(%2),%0" "\n\t"
		"movl %0,12(%3)" "\n\t"
		"movl 8(%1),%0"  "\n\t"
		"adcl 8(%2),%0"  "\n\t"
		"movl %0,8(%3)"  "\n\t"
		"movl 4(%1),%0"  "\n\t"
		"adcl 4(%2),%0"  "\n\t"
		"movl %0,4(%3)"  "\n\t"
		"movl (%1),%0"   "\n\t"
		"adcl (%2),%0"   "\n\t"
		"movl %0,(%3)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b), "r" (&r)
		: "cc"
		);
    #else
	__asm__ __volatile__ (
		"movl (%1),%0"   "\n\t"
		"addl (%2),%0"   "\n\t"
		"movl %0,(%3)"   "\n\t"
		"movl 4(%1),%0"  "\n\t"
		"adcl 4(%2),%0"  "\n\t"
		"movl %0,4(%3)"  "\n\t"
		"movl 8(%1),%0"  "\n\t"
		"adcl 8(%2),%0"  "\n\t"
		"movl %0,8(%3)"  "\n\t"
		"movl 12(%1),%0" "\n\t"
		"adcl 12(%2),%0" "\n\t"
		"movl %0,12(%3)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b), "r" (&r)
		: "cc"
		);
    #endif
  #endif
#elif defined(NUSS_OUT_EXTERNAL_LOOPS)
	add_loop_lsp(arrayLSDptr(a._ow,4),arrayLSDptr(b._ow,4),arrayLSDptr(r._ow,4),4);
#else
	var uint32 tmp;

	tmp = a.ow0 + b.ow0;
	if (tmp >= a.ow0) {
		// no carry
		r.ow0 = tmp;
		tmp = a.ow1 + b.ow1;
		if (tmp >= a.ow1) goto no_carry_1; else goto carry_1;
	} else {
		// carry
		r.ow0 = tmp;
		tmp = a.ow1 + b.ow1 + 1;
		if (tmp > a.ow1) goto no_carry_1; else goto carry_1;
	}
	if (1) {
		no_carry_1: // no carry
		r.ow1 = tmp;
		tmp = a.ow2 + b.ow2;
		if (tmp >= a.ow2) goto no_carry_2; else goto carry_2;
	} else {
		carry_1: // carry
		r.ow1 = tmp;
		tmp = a.ow2 + b.ow2 + 1;
		if (tmp > a.ow2) goto no_carry_2; else goto carry_2;
	}
	if (1) {
		no_carry_2: // no carry
		r.ow2 = tmp;
		tmp = a.ow3 + b.ow3;
	} else {
		carry_2: // carry
		r.ow2 = tmp;
		tmp = a.ow3 + b.ow3 + 1;
	}
	r.ow3 = tmp;
#endif
}

// r := a - b
static inline void sub (const nuss_outword& a, const nuss_outword& b, nuss_outword& r)
{
#if defined(__GNUC__) && defined(__i386__)
	var uintD dummy;
  #ifdef NUSS_ASM_DIRECT
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"subl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.ow0), "m" (b.ow0), "m" (r.ow0)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"sbbl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.ow1), "m" (b.ow1), "m" (r.ow1)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"sbbl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.ow2), "m" (b.ow2), "m" (r.ow2)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"sbbl %2,%0" "\n\t"
		"movl %0,%3"
		: "=&q" (dummy)
		: "m" (a.ow3), "m" (b.ow3), "m" (r.ow3)
		: "cc"
		);
  #else
    #if CL_DS_BIG_ENDIAN_P
	__asm__ __volatile__ (
		"movl 12(%1),%0" "\n\t"
		"subl 12(%2),%0" "\n\t"
		"movl %0,12(%3)" "\n\t"
		"movl 8(%1),%0"  "\n\t"
		"sbbl 8(%2),%0"  "\n\t"
		"movl %0,8(%3)"  "\n\t"
		"movl 4(%1),%0"  "\n\t"
		"sbbl 4(%2),%0"  "\n\t"
		"movl %0,4(%3)"  "\n\t"
		"movl (%1),%0"   "\n\t"
		"sbbl (%2),%0"   "\n\t"
		"movl %0,(%3)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b), "r" (&r)
		: "cc"
		);
    #else
	__asm__ __volatile__ (
		"movl (%1),%0"   "\n\t"
		"subl (%2),%0"   "\n\t"
		"movl %0,(%3)"   "\n\t"
		"movl 4(%1),%0"  "\n\t"
		"sbbl 4(%2),%0"  "\n\t"
		"movl %0,4(%3)"  "\n\t"
		"movl 8(%1),%0"  "\n\t"
		"sbbl 8(%2),%0"  "\n\t"
		"movl %0,8(%3)"  "\n\t"
		"movl 12(%1),%0" "\n\t"
		"sbbl 12(%2),%0" "\n\t"
		"movl %0,12(%3)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b), "r" (&r)
		: "cc"
		);
    #endif
  #endif
#elif defined(NUSS_OUT_EXTERNAL_LOOPS)
	sub_loop_lsp(arrayLSDptr(a._ow,4),arrayLSDptr(b._ow,4),arrayLSDptr(r._ow,4),4);
#else
	var uint32 tmp;

	tmp = a.ow0 - b.ow0;
	if (tmp <= a.ow0) {
		// no carry
		r.ow0 = tmp;
		tmp = a.ow1 - b.ow1;
		if (tmp <= a.ow1) goto no_carry_1; else goto carry_1;
	} else {
		// carry
		r.ow0 = tmp;
		tmp = a.ow1 - b.ow1 - 1;
		if (tmp < a.ow1) goto no_carry_1; else goto carry_1;
	}
	if (1) {
		no_carry_1: // no carry
		r.ow1 = tmp;
		tmp = a.ow2 - b.ow2;
		if (tmp <= a.ow2) goto no_carry_2; else goto carry_2;
	} else {
		carry_1: // carry
		r.ow1 = tmp;
		tmp = a.ow2 - b.ow2 - 1;
		if (tmp < a.ow2) goto no_carry_2; else goto carry_2;
	}
	if (1) {
		no_carry_2: // no carry
		r.ow2 = tmp;
		tmp = a.ow3 - b.ow3;
	} else {
		carry_2: // carry
		r.ow2 = tmp;
		tmp = a.ow3 - b.ow3 - 1;
	}
	r.ow3 = tmp;
#endif
}

// b := a >> 1
static inline void shift (const nuss_outword& a, nuss_outword& b)
{
#if defined(__GNUC__) && defined(__i386__) && !defined(DEBUG_NUSS)
	var uintD dummy;
  #ifdef NUSS_ASM_DIRECT
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"sarl $1,%0" "\n\t"
		"movl %0,%2"
		: "=&q" (dummy)
		: "m" (a.ow3), "m" (b.ow3)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"rcrl $1,%0" "\n\t"
		"movl %0,%2"
		: "=&q" (dummy)
		: "m" (a.ow2), "m" (b.ow2)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"rcrl $1,%0" "\n\t"
		"movl %0,%2"
		: "=&q" (dummy)
		: "m" (a.ow1), "m" (b.ow1)
		: "cc"
		);
	__asm__ __volatile__ (
		"movl %1,%0" "\n\t"
		"rcrl $1,%0" "\n\t"
		"movl %0,%2"
		: "=&q" (dummy)
		: "m" (a.ow0), "m" (b.ow0)
		: "cc"
		);
  #else
    #if CL_DS_BIG_ENDIAN_P
	__asm__ __volatile__ (
		"movl (%1),%0"   "\n\t"
		"sarl $1,%0"     "\n\t"
		"movl %0,(%2)"   "\n\t"
		"movl 4(%1),%0"  "\n\t"
		"rcrl $1,%0"     "\n\t"
		"movl %0,4(%2)"  "\n\t"
		"movl 8(%1),%0"  "\n\t"
		"rcrl $1,%0"     "\n\t"
		"movl %0,8(%2)"  "\n\t"
		"movl 12(%1),%0" "\n\t"
		"rcrl $1,%0"     "\n\t"
		"movl %0,12(%2)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b)
		: "cc"
		);
    #else
	__asm__ __volatile__ (
		"movl 12(%1),%0" "\n\t"
		"sarl $1,%0"     "\n\t"
		"movl %0,12(%2)" "\n\t"
		"movl 8(%1),%0"  "\n\t"
		"rcrl $1,%0"     "\n\t"
		"movl %0,8(%2)"  "\n\t"
		"movl 4(%1),%0"  "\n\t"
		"rcrl $1,%0"     "\n\t"
		"movl %0,4(%2)"  "\n\t"
		"movl (%1),%0"   "\n\t"
		"rcrl $1,%0"     "\n\t"
		"movl %0,(%2)"
		: "=&q" (dummy)
		: "r" (&a), "r" (&b)
		: "cc"
		);
    #endif
  #endif
#elif defined(NUSS_OUT_EXTERNAL_LOOPS)
	#ifdef DEBUG_NUSS
	if (shiftrightcopy_loop_msp(arrayMSDptr(a._ow,4),arrayMSDptr(b._ow,4),4,1,mspref(arrayMSDptr(a._ow,4),0)>>31))
		throw runtime_exception();
	#else
	shiftrightcopy_loop_msp(arrayMSDptr(a._ow,4),arrayMSDptr(b._ow,4),4,1,mspref(arrayMSDptr(a._ow,4),0)>>31);
	#endif
#else
	var uint32 tmp, carry;

	tmp = a.ow3;
	b.ow3 = (sint32)tmp >> 1;
	carry = tmp << 31;
	tmp = a.ow2;
	b.ow2 = (tmp >> 1) | carry;
	carry = tmp << 31;
	tmp = a.ow1;
	b.ow1 = (tmp >> 1) | carry;
	carry = tmp << 31;
	tmp = a.ow0;
	b.ow0 = (tmp >> 1) | carry;
	#ifdef DEBUG_NUSS
	carry = tmp << 31;
	if (carry)
		throw runtime_exception();
	#endif
#endif
}

#endif // (intDsize==32)

#if (intDsize==64)

//typedef struct { sint64 iw1; uint64 iw0; } nuss_inword;
//typedef struct { uint64 iw0; sint64 iw1; } nuss_inword;
typedef struct { uintD _iw[2]; } nuss_inword;
#if CL_DS_BIG_ENDIAN_P
  #define iw1 _iw[0]
  #define iw0 _iw[1]
#else
  #define iw0 _iw[0]
  #define iw1 _iw[1]
#endif

//typedef struct { sint64 ow2; uint64 ow1; uint64 ow0; } nuss_outword;
//typedef struct { uint64 ow0; uint64 ow1; sint64 ow2; } nuss_outword;
typedef struct { uintD _ow[3]; } nuss_outword;
#if CL_DS_BIG_ENDIAN_P
  #define ow2 _ow[0]
  #define ow1 _ow[1]
  #define ow0 _ow[2]
#else
  #define ow0 _ow[0]
  #define ow1 _ow[1]
  #define ow2 _ow[2]
#endif

// r := a + b
static inline void add (const nuss_inword& a, const nuss_inword& b, nuss_inword& r)
{
#ifdef NUSS_IN_EXTERNAL_LOOPS
	add_loop_lsp(arrayLSDptr(a._iw,2),arrayLSDptr(b._iw,2),arrayLSDptr(r._iw,2),2);
#else
	var uint64 tmp;

	tmp = a.iw0 + b.iw0;
	if (tmp >= a.iw0) {
		// no carry
		r.iw0 = tmp;
		r.iw1 = a.iw1 + b.iw1;
	} else {
		// carry
		r.iw0 = tmp;
		r.iw1 = a.iw1 + b.iw1 + 1;
	}
#endif
}

// r := a - b
static inline void sub (const nuss_inword& a, const nuss_inword& b, nuss_inword& r)
{
#ifdef NUSS_IN_EXTERNAL_LOOPS
	sub_loop_lsp(arrayLSDptr(a._iw,2),arrayLSDptr(b._iw,2),arrayLSDptr(r._iw,2),2);
#else
	var uint64 tmp;

	tmp = a.iw0 - b.iw0;
	if (tmp <= a.iw0) {
		// no carry
		r.iw0 = tmp;
		r.iw1 = a.iw1 - b.iw1;
	} else {
		// carry
		r.iw0 = tmp;
		r.iw1 = a.iw1 - b.iw1 - 1;
	}
#endif
}

// r := 0
static inline void zero (nuss_outword& r)
{
	r.ow0 = 0;
	r.ow1 = 0;
	r.ow2 = 0;
}

// r := a + b
static inline void add (const nuss_outword& a, const nuss_outword& b, nuss_outword& r)
{
#ifdef NUSS_OUT_EXTERNAL_LOOPS
	add_loop_lsp(arrayLSDptr(a._ow,3),arrayLSDptr(b._ow,3),arrayLSDptr(r._ow,3),3);
#else
	var uint64 tmp;

	tmp = a.ow0 + b.ow0;
	if (tmp >= a.ow0) {
		// no carry
		r.ow0 = tmp;
		tmp = a.ow1 + b.ow1;
		if (tmp >= a.ow1) goto no_carry_1; else goto carry_1;
	} else {
		// carry
		r.ow0 = tmp;
		tmp = a.ow1 + b.ow1 + 1;
		if (tmp > a.ow1) goto no_carry_1; else goto carry_1;
	}
	if (1) {
		no_carry_1: // no carry
		r.ow1 = tmp;
		tmp = a.ow2 + b.ow2;
	} else {
		carry_1: // carry
		r.ow1 = tmp;
		tmp = a.ow2 + b.ow2 + 1;
	}
	r.ow2 = tmp;
#endif
}

// r := a - b
static inline void sub (const nuss_outword& a, const nuss_outword& b, nuss_outword& r)
{
#ifdef NUSS_OUT_EXTERNAL_LOOPS
	sub_loop_lsp(arrayLSDptr(a._ow,3),arrayLSDptr(b._ow,3),arrayLSDptr(r._ow,3),3);
#else
	var uint64 tmp;

	tmp = a.ow0 - b.ow0;
	if (tmp <= a.ow0) {
		// no carry
		r.ow0 = tmp;
		tmp = a.ow1 - b.ow1;
		if (tmp <= a.ow1) goto no_carry_1; else goto carry_1;
	} else {
		// carry
		r.ow0 = tmp;
		tmp = a.ow1 - b.ow1 - 1;
		if (tmp < a.ow1) goto no_carry_1; else goto carry_1;
	}
	if (1) {
		no_carry_1: // no carry
		r.ow1 = tmp;
		tmp = a.ow2 - b.ow2;
	} else {
		carry_1: // carry
		r.ow1 = tmp;
		tmp = a.ow2 - b.ow2 - 1;
	}
	r.ow2 = tmp;
#endif
}

// b := a >> 1
static inline void shift (const nuss_outword& a, nuss_outword& b)
{
#ifdef NUSS_OUT_EXTERNAL_LOOPS
	#ifdef DEBUG_NUSS
	if (shiftrightcopy_loop_msp(arrayMSDptr(a._ow,3),arrayMSDptr(b._ow,3),3,1,mspref(arrayMSDptr(a._ow,3),0)>>63))
		throw runtime_exception();
	#else
	shiftrightcopy_loop_msp(arrayMSDptr(a._ow,3),arrayMSDptr(b._ow,3),3,1,mspref(arrayMSDptr(a._ow,3),0)>>63);
	#endif
#else
	var uint64 tmp, carry;

	tmp = a.ow2;
	b.ow2 = (sint64)tmp >> 1;
	carry = tmp << 63;
	tmp = a.ow1;
	b.ow1 = (tmp >> 1) | carry;
	carry = tmp << 63;
	tmp = a.ow0;
	b.ow0 = (tmp >> 1) | carry;
	#ifdef DEBUG_NUSS
	carry = tmp << 63;
	if (carry)
		throw runtime_exception();
	#endif
#endif
}

#endif // (intDsize==64)

// This is a recursive implementation.
// TODO: Write a non-recursive one.

#ifndef _BIT_REVERSE
#define _BIT_REVERSE
// Reverse an n-bit number x. n>0.
static uintC bit_reverse (uintL n, uintC x)
{
	var uintC y = 0;
	do {
		y <<= 1;
		y |= (x & 1);
		x >>= 1;
	} while (!(--n == 0));
	return y;
}
#endif

// Threshold for recursion base in mulu_nuss_negacyclic().
// Time of a multiplication with len1=len2=10000 on Linux i486:
//                     normal    asm-optimized
//   threshold1 = 1:   40.1 sec  25.5 sec
//   threshold1 = 2:   28.6 sec  18.3 sec
//   threshold1 = 3:   25.6 sec  16.6 sec
//   threshold1 = 4:   25.7 sec  17.6 sec
//   threshold1 = 5:   26.1 sec  18.0 sec
const uintC cl_nuss_threshold1 = 3;

// Threshold for recursion base in mulu_nuss_cyclic().
const uintC cl_nuss_threshold2 = 1;

// Computes z[k] := sum(i+j==k mod N, x[i]*y[j]*(-1)^((i+j-k)/N))
// for all k=0..N-1.
static void mulu_nuss_negacyclic (const uintL n, const uintC N, // N = 2^n
                                  const nuss_inword * x, // N words
                                  const nuss_inword * y, // N words
                                  nuss_outword * z       // N words result
                                 )
{
	#if 0 // always n > 0
	if (n == 0) {
		// z[0] := x0 y0
		mul(x[0],y[0], z[0]);
		return;
	}
	#endif
	if (n <= cl_nuss_threshold1) {
		if (n == 1) {
			// z[0] := x0 (y0 + y1) - (x0 + x1) y1
			// z[1] := x0 (y0 + y1) + (x1 - x0) y0
			var nuss_inword x_sum;
			var nuss_inword y_sum;
			var nuss_outword first, second;
			add(x[0],x[1], x_sum);
			add(y[0],y[1], y_sum);
			mul(x[0],y_sum, first);
			mul(x_sum,y[1], second); sub(first,second, z[0]);
			sub(x[1],x[0], x_sum);
			mul(x_sum,y[0], second); add(first,second, z[1]);
			return;
		}
		// 1 < n <= cl_nuss_threshold1.
		#if 0 // straightforward, but slow
		var uintC k;
		for (k = 0; k < N; k++) {
			var uintC i;
			var nuss_outword accu;
			mul(x[0],y[k], accu);
			for (i = 1; i <= k; i++) {
				var nuss_outword temp;
				mul(x[i],y[k-i], temp);
				add(accu,temp, accu);
			}
			for (i = k+1; i < N; i++) {
				var nuss_outword temp;
				mul(x[i],y[N-i+k], temp);
				sub(accu,temp, accu);
			}
			z[k] = accu;
		}
		#else
		var const uintC M = (uintC)1 << (n-1); // M = N/2
		var uintC i, j, k;
		for (k = 0; k < N; k++)
			zero(z[k]);
		for (i = 0; i < M; i++) {
			var uintC iM = i+M;
			for (j = 0; j < M-i; j++) {
				var uintC jM = j+M;
				// z[i+j]   += x[i] (y[j] + y[j+M]) - (x[i] + x[i+M]) y[j+M]
				// z[i+j+M] += x[i] (y[j] + y[j+M]) + (x[i+M] - x[i]) y[j]
				var nuss_inword x_sum;
				var nuss_inword y_sum;
				var nuss_outword first, second, temp;
				add(x[i],x[iM], x_sum);
				add(y[j],y[jM], y_sum);
				mul(x[i],y_sum, first);
				mul(x_sum,y[jM], second); sub(first,second, temp); add(z[i+j],temp, z[i+j]);
				sub(x[iM],x[i], x_sum);
				mul(x_sum,y[j], second); add(first,second, temp); add(z[i+j+M],temp, z[i+j+M]);
			}
			for (j = M-i; j < M; j++) {
				var uintC jM = j+M;
				// z[i+j]   += x[i] (y[j] + y[j+M]) - (x[i] + x[i+M]) y[j+M]
				// z[i+j-M] -= x[i] (y[j] + y[j+M]) + (x[i+M] - x[i]) y[j]
				var nuss_inword x_sum;
				var nuss_inword y_sum;
				var nuss_outword first, second, temp;
				add(x[i],x[iM], x_sum);
				add(y[j],y[jM], y_sum);
				mul(x[i],y_sum, first);
				mul(x_sum,y[jM], second); sub(first,second, temp); add(z[i+j],temp, z[i+j]);
				sub(x[iM],x[i], x_sum);
				mul(x_sum,y[j], second); add(first,second, temp); sub(z[i+j-M],temp, z[i+j-M]);
			}
		}
		#endif
		return;
	}
	// Recursive FFT.
	var const uintL m = n >> 1; // floor(n/2)
	var const uintL r = n - m;  // ceiling(n/2)
	var const uintC M = (uintC)1 << m; // M = 2^m
	var const uintC R = (uintC)1 << r; // R = 2^r
	CL_ALLOCA_STACK;
	var nuss_inword* const auX = cl_alloc_array(nuss_inword,2*N);
	var nuss_inword* const auY = cl_alloc_array(nuss_inword,2*N);
	var nuss_outword* const auZ = cl_alloc_array(nuss_outword,2*N);
	#define X(i,j) auX[((i)<<r)+(j)] /* 0 <= i < 2*M, 0 <= j < R */
	#define Y(i,j) auY[((i)<<r)+(j)] /* 0 <= i < 2*M, 0 <= j < R */
	#define Z(i,j) auZ[((i)<<r)+(j)] /* 0 <= i < 2*M, 0 <= j < R */
	var nuss_inword* const tmp1 = cl_alloc_array(nuss_inword,R);
	var nuss_inword* const tmp2 = cl_alloc_array(nuss_inword,R);
	var nuss_outword* const tmpZ = cl_alloc_array(nuss_outword,R);
	var bool squaring = (x == y);
	var uintC i, j;
	// Initialize polynomials X(i) and Y(i).
	for (i = 0; i < M; i++) {
		{
			for (j = 0; j < R; j++)
				X(i,j) = x[(j<<m) + i];
		}
		if (!squaring) {
			for (j = 0; j < R; j++)
				Y(i,j) = y[(j<<m) + i];
		}
	}
	// For i = M..2*M-1, the polynomials are implicitly 0.
	// Do an FFT of length 2*M on X.
	{
		var sintL l;
		// Level l = m:
		for (i = 0; i < M; i++)
			for (j = 0; j < R; j++)
				X(i+M,j) = X(i,j);
		// Level l = m-1..0:
		for (l = m-1; l>=0; l--) {
			var const uintC smax = (uintC)1 << (m-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(m-l,s) << (l + r-m);
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Butterfly: replace (X(i1),X(i2)) by
					// (X(i1) + w^exp*X(i2), X(i1) - w^exp*X(i2)).
					for (j = 0; j < exp; j++) {
						// note that w^R = -1
						sub(X(i1,j),X(i2,j-exp+R), tmp1[j]);
						add(X(i1,j),X(i2,j-exp+R), tmp2[j]);
					}
					for (j = exp; j < R; j++) {
						add(X(i1,j),X(i2,j-exp), tmp1[j]);
						sub(X(i1,j),X(i2,j-exp), tmp2[j]);
					}
					for (j = 0; j < R; j++) {
						X(i1,j) = tmp1[j];
						X(i2,j) = tmp2[j];
					}
				}
			}
		}
	}
	// Do an FFT of length 2*M on Y.
	if (!squaring) {
		var sintL l;
		// Level l = m:
		for (i = 0; i < M; i++)
			for (j = 0; j < R; j++)
				Y(i+M,j) = Y(i,j);
		// Level l = m-1..0:
		for (l = m-1; l>=0; l--) {
			var const uintC smax = (uintC)1 << (m-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(m-l,s) << (l + r-m);
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Butterfly: replace (Y(i1),Y(i2)) by
					// (Y(i1) + w^exp*Y(i2), Y(i1) - w^exp*Y(i2)).
					for (j = 0; j < exp; j++) {
						// note that w^R = -1
						sub(Y(i1,j),Y(i2,j-exp+R), tmp1[j]);
						add(Y(i1,j),Y(i2,j-exp+R), tmp2[j]);
					}
					for (j = exp; j < R; j++) {
						add(Y(i1,j),Y(i2,j-exp), tmp1[j]);
						sub(Y(i1,j),Y(i2,j-exp), tmp2[j]);
					}
					for (j = 0; j < R; j++) {
						Y(i1,j) = tmp1[j];
						Y(i2,j) = tmp2[j];
					}
				}
			}
		}
	}
	// Recursively compute the negacyclic product X(i)*Y(i) for all i.
	if (!squaring) {
		for (i = 0; i < 2*M; i++)
			mulu_nuss_negacyclic(r,R, &X(i,0), &Y(i,0), &Z(i,0));
	} else {
		for (i = 0; i < 2*M; i++)
			mulu_nuss_negacyclic(r,R, &X(i,0), &X(i,0), &Z(i,0));
	}
	// Undo an FFT of length 2*M on Z.
	{
		var uintL l;
		// Level l = 0..m-1:
		for (l = 0; l < m; l++) {
			var const uintC smax = (uintC)1 << (m-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(m-l,s) << (l + r-m);
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (Z(i1),Z(i2)) by
					// ((Z(i1)+Z(i2))/2, (Z(i1)-Z(i2))/(2*w^exp)).
					for (j = 0; j < exp; j++)
						// note that w^R = -1
						sub(Z(i2,j),Z(i1,j), tmpZ[j-exp+R]);
					for (j = exp; j < R; j++)
						sub(Z(i1,j),Z(i2,j), tmpZ[j-exp]);
					for (j = 0; j < R; j++) {
						var nuss_outword sum;
						add(Z(i1,j),Z(i2,j), sum);
						shift(sum, Z(i1,j));
						shift(tmpZ[j], Z(i2,j));
					}
				}
			}
		}
		// Level l=m:
		for (i = 0; i < M; i++) {
			var uintC i1 = i;
			var uintC i2 = i1 + M;
			// Inverse Butterfly: replace (Z(i1),Z(i2)) by
			// ((Z(i1)+Z(i2))/2, (Z(i1)-Z(i2))/2).
			for (j = 0; j < R; j++) {
				var nuss_outword sum;
				var nuss_outword diff;
				add(Z(i1,j),Z(i2,j), sum);
				sub(Z(i1,j),Z(i2,j), diff);
				shift(sum, Z(i1,j));
				shift(diff, Z(i2,j));
			}
		}
	}
	// Reduce to length M.
	for (i = 0; i < M; i++) {
		sub(Z(i,0),Z(i+M,R-1), z[i]);
		for (j = 1; j < R; j++)
			add(Z(i,j),Z(i+M,j-1), z[(j<<m)+i]);
	}
	#undef Z
	#undef Y
	#undef X
}

// Computes z[k] := sum(i+j==k mod N, x[i]*y[j])
// for all k=0..N-1.
static void mulu_nuss_cyclic (const uintL n, const uintC N, // N = 2^n
                              nuss_inword * x, // N words, modified!
                              nuss_inword * y, // N words, modified!
                              nuss_outword * z // N words result
                             )
{
	unused N;
	#if 0 // always n > 0
	if (n == 0) {
		// z[0] := x0 y0
		mul(x[0],y[0], z[0]);
		return;
	}
	#endif
	if (n == 1) {
		// z[0] := ((x0 + x1) (y0 + y1) + (x0 - x1) (y0 - y1)) / 2
		// z[1] := ((x0 + x1) (y0 + y1) - (x0 - x1) (y0 - y1)) / 2
		var nuss_inword x_sum;
		var nuss_inword y_sum;
		var nuss_inword x_diff;
		var nuss_inword y_diff;
		var nuss_outword first, second;
		add(x[0],x[1], x_sum);
		add(y[0],y[1], y_sum);
		sub(x[0],x[1], x_diff);
		sub(y[0],y[1], y_diff);
		mul(x_sum,y_sum, first);
		mul(x_diff,y_diff, second);
		add(first,second, z[0]); shift(z[0], z[0]);
		sub(first,second, z[1]); shift(z[1], z[1]);
		return;
	}
	#if 0 // useless code because cl_nuss_threshold2 == 1
	if (n <= cl_nuss_threshold2) {
		#if 0 // straightforward, but slow
		var uintC k;
		for (k = 0; k < N; k++) {
			var uintC i;
			var nuss_outword accu;
			mul(x[0],y[k], accu);
			for (i = 1; i <= k; i++) {
				var nuss_outword temp;
				mul(x[i],y[k-i], temp);
				add(accu,temp, accu);
			}
			for (i = k+1; i < N; i++) {
				var nuss_outword temp;
				mul(x[i],y[N-i+k], temp);
				add(accu,temp, accu);
			}
			z[k] = accu;
		}
		#else
		var const uintC M = (uintC)1 << (n-1); // M = N/2
		var uintC i, j, k;
		for (k = 0; k < N; k++)
			zero(z[k]);
		for (i = 0; i < M; i++) {
			var uintC iM = i+M;
			for (j = 0; j < M; j++) {
				var uintC jM = j+M;
				// z[i+j]   += ((x[i] + x[i+M]) (y[j] + y[j+M]) + (x[i] - x[i+M]) (y[j] - y[j+M])) / 2
				// z[i+j+M] += ((x[i] + x[i+M]) (y[j] + y[j+M]) - (x[i] - x[i+M]) (y[j] - y[j+M])) / 2
				var nuss_inword x_sum;
				var nuss_inword y_sum;
				var nuss_inword x_diff;
				var nuss_inword y_diff;
				var nuss_outword first, second, temp;
				add(x[i],x[iM], x_sum);
				add(y[j],y[jM], y_sum);
				sub(x[i],x[iM], x_diff);
				sub(y[j],y[jM], y_diff);
				mul(x_sum,y_sum, first);
				mul(x_diff,y_diff, second);
				add(first,second, temp); add(z[i+j],temp, z[i+j]);
				var uintC ijM = (i+j+M) & (N-1);
				sub(first,second, temp); add(z[ijM],temp, z[ijM]);
			}
		}
		for (k = 0; k < N; k++)
			shift(z[k], z[k]);
		#endif
		return;
	}
	#endif
	var const uintL m = n-1;
	var const uintC M = (uintC)1 << m; // M = 2^m = N/2
	var uintC i;
	// Chinese remainder theorem: u^N-1 = (u^M-1)*(u^M+1)
	for (i = 0; i < M; i++) {
		// Butterfly: replace (x(i),x(i+M))
		// by (x(i)+x(i+M),x(i)-x(i+M)).
		var nuss_inword tmp;
		sub(x[i],x[i+M], tmp);
		add(x[i],x[i+M], x[i]);
		x[i+M] = tmp;
	}
	if (!(x == y)) // squaring?
	for (i = 0; i < M; i++) {
		// Butterfly: replace (y(i),y(i+M))
		// by (y(i)+y(i+M),y(i)-y(i+M)).
		var nuss_inword tmp;
		sub(y[i],y[i+M], tmp);
		add(y[i],y[i+M], y[i]);
		y[i+M] = tmp;
	}
	// Recurse.
	mulu_nuss_cyclic(m,M, &x[0], &y[0], &z[0]);
	mulu_nuss_negacyclic(m,M, &x[M], &y[M], &z[M]);
	for (i = 0; i < M; i++) {
		// Inverse Butterfly: replace (z(i),z(i+M))
		// by ((z(i)+z(i+M))/2,(z(i)-z(i+M))/2).
		var nuss_outword sum;
		var nuss_outword diff;
		add(z[i],z[i+M], sum);
		sub(z[i],z[i+M], diff);
		shift(sum, z[i]);
		shift(diff, z[i+M]);
	}
}

static void mulu_nussbaumer (const uintD* sourceptr1, uintC len1,
                             const uintD* sourceptr2, uintC len2,
                             uintD* destptr)
// Es ist 2 <= len1 <= len2.
{
	// Methode:
	// source1 ist ein Stück der Länge N1, source2 ein oder mehrere Stücke
	// der Länge N2, mit N1+N2 <= N, wobei N Zweierpotenz ist.
	// sum(i=0..N-1, x_i b^i) * sum(i=0..N-1, y_i b^i) wird errechnet,
	// indem man die beiden Polynome
	// sum(i=0..N-1, x_i T^i), sum(i=0..N-1, y_i T^i)
	// multipliziert, und zwar durch Fourier-Transformation (s.o.).
	var uint32 n;
	integerlengthC(len1-1, n=); // 2^(n-1) < len1 <= 2^n
	var uintC len = (uintC)1 << n; // kleinste Zweierpotenz >= len1
	// Wählt man N = len, so hat man ceiling(len2/(len-len1+1)) * FFT(len).
	// Wählt man N = 2*len, so hat man ceiling(len2/(2*len-len1+1)) * FFT(2*len).
	// Wir wählen das billigere von beiden:
	// Bei ceiling(len2/(len-len1+1)) <= 2 * ceiling(len2/(2*len-len1+1))
	// nimmt man N = len, bei ....... > ........ dagegen N = 2*len.
	// (Wahl von N = 4*len oder mehr bringt nur in Extremfällen etwas.)
	if (len2 > 2 * (len-len1+1) * (len2 <= (2*len-len1+1) ? 1 : ceiling(len2,(2*len-len1+1)))) {
		n = n+1;
		len = len << 1;
	}
	var const uintC N = len; // N = 2^n
	CL_ALLOCA_STACK;
	var nuss_inword* const x = cl_alloc_array(nuss_inword,N);
	var nuss_inword* const y = cl_alloc_array(nuss_inword,N);
	var nuss_outword* const z = cl_alloc_array(nuss_outword,N);
	var uintD* const tmpprod = cl_alloc_array(uintD,len1+1);
	var uintP i;
	var uintC destlen = len1+len2;
	clear_loop_lsp(destptr,destlen);
	do {
		var uintC len2p; // length of a piece of source2
		len2p = N - len1 + 1;
		if (len2p > len2)
			len2p = len2;
		// len2p = min(N-len1+1,len2).
		if (len2p == 1) {
			// cheap case
			var uintD* tmpptr = arrayLSDptr(tmpprod,len1+1);
			mulu_loop_lsp(lspref(sourceptr2,0),sourceptr1,tmpptr,len1);
			if (addto_loop_lsp(tmpptr,destptr,len1+1))
				if (inc_loop_lsp(destptr lspop (len1+1),destlen-(len1+1)))
					throw runtime_exception();
		} else {
			var uintC destlenp = len1 + len2p - 1;
			// destlenp = min(N,destlen-1).
			var bool squaring = ((sourceptr1 == sourceptr2) && (len1 == len2p));
			// Fill factor x.
			{
				for (i = 0; i < len1; i++) {
					x[i].iw0 = lspref(sourceptr1,i);
					x[i].iw1 = 0;
				}
				for (i = len1; i < N; i++) {
					x[i].iw0 = 0;
					x[i].iw1 = 0;
				}
			}
			// Fill factor y.
			if (!squaring) {
				for (i = 0; i < len2p; i++) {
					y[i].iw0 = lspref(sourceptr2,i);
					y[i].iw1 = 0;
				}
				for (i = len2p; i < N; i++) {
					y[i].iw0 = 0;
					y[i].iw1 = 0;
				}
			}
			// Multiply.
			if (!squaring)
				mulu_nuss_cyclic(n,N, &x[0], &y[0], &z[0]);
			else
				mulu_nuss_cyclic(n,N, &x[0], &x[0], &z[0]);
			#ifdef DEBUG_NUSS
			// Check result.
			for (i = 0; i < N; i++)
				if (!(z[i].ow3 == 0))
					throw runtime_exception();
			#endif
			// Add result to destptr[-destlen..-1]:
			{
				var uintD* ptr = destptr;
				// ac2|ac1|ac0 are an accumulator.
				var uint32 ac0 = 0;
				var uint32 ac1 = 0;
				var uint32 ac2 = 0;
				var uint32 tmp;
				for (i = 0; i < destlenp; i++) {
					// Add z[i] to the accumulator.
					tmp = z[i].ow0;
					if ((ac0 += tmp) < tmp) {
						if (++ac1 == 0)
							++ac2;
					}
					tmp = z[i].ow1;
					if ((ac1 += tmp) < tmp)
						++ac2;
					tmp = z[i].ow2;
					ac2 += tmp;
					// Add the accumulator's least significant word to destptr:
					tmp = lspref(ptr,0);
					if ((ac0 += tmp) < tmp) {
						if (++ac1 == 0)
							++ac2;
					}
					lspref(ptr,0) = ac0;
					lsshrink(ptr);
					ac0 = ac1;
					ac1 = ac2;
					ac2 = 0;
				}
				// ac2 = 0.
				if (ac1 > 0) {
					if (!((i += 2) <= destlen))
						throw runtime_exception();
					tmp = lspref(ptr,0);
					if ((ac0 += tmp) < tmp)
						++ac1;
					lspref(ptr,0) = ac0;
					lsshrink(ptr);
					tmp = lspref(ptr,0);
					ac1 += tmp;
					lspref(ptr,0) = ac1;
					lsshrink(ptr);
					if (ac1 < tmp)
						if (inc_loop_lsp(ptr,destlen-i))
							throw runtime_exception();
				} else if (ac0 > 0) {
					if (!((i += 1) <= destlen))
						throw runtime_exception();
					tmp = lspref(ptr,0);
					ac0 += tmp;
					lspref(ptr,0) = ac0;
					lsshrink(ptr);
					if (ac0 < tmp)
						if (inc_loop_lsp(ptr,destlen-i))
							throw runtime_exception();
				}
			}
			#ifdef DEBUG_NUSS
			// If destlenp < N, check that the remaining z[i] are 0.
			for (i = destlenp; i < N; i++)
				if (z[i].ow2 > 0 || z[i].ow1 > 0 || z[i].ow0 > 0)
					throw runtime_exception();
			#endif
		}
		// Decrement len2.
		destptr = destptr lspop len2p;
		destlen -= len2p;
		sourceptr2 = sourceptr2 lspop len2p;
		len2 -= len2p;
	} while (len2 > 0);
}

#undef iw0
#undef iw1
#undef ow0
#undef ow1
#undef ow2
#undef ow3
