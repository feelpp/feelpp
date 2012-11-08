// Fast integer multiplication using FFT in a modular ring.
// Bruno Haible 14.5.,16.5.1996

// FFT in the complex domain has the drawback that it needs careful round-off
// error analysis. So here we choose another field of characteristic 0: Q_p.
// Since Q_p contains exactly the (p-1)th roots of unity, we choose
// p == 1 mod N and have the Nth roots of unity (N = 2^n) in Q_p and
// even in Z_p. Actually, we compute in Z/(p^m Z).

// All operations the FFT algorithm needs is addition, subtraction,
// multiplication, multiplication by the Nth root and unity and division
// by N. Hence we can use the domain Z/(p^m Z) even if p is not a prime!

// We use the Schönhage-Strassen choice of the modulus: p = 2^R+1. This
// has two big advantages: Multiplication and division by 2 (which is a
// (2R)th root of unity) or a power of 2 is just a shift and an subtraction.
// And multiplication mod p is just a normal multiplication, followed by
// a subtraction.
// In order to exploit the (2R)th root of unity for FFT, we choose R = 2^r,
// and do an FFT of size M with M = 2^m and M | 2R.

// Say we want to compute the product of two integers with N1 and N2 bits,
// respectively. We choose N >= N1+N2 and K, R, M with
//      ceiling(N1/K)+ceiling(N2/K)-1 <= M  (i.e. roughly N <= K*M),
//      2*K+ceiling(log2(M)) <= R,
//      R = 2^r, M = 2^m, M | 2R.
// We then split each of the factors in M K-bit chunks each, and do
// an FFT mod p = 2^R+1. We then recover the convolution of the chunks
// from the FFT product (the first inequality ensures that this is possible).
// The second inequality ensures that we have no overflow, i.e. the
// convolution result is valid in Z, not only in Z/pZ.

// The computation time (bit complexity) will be proportional to
//    Mul(N) = O(M log(M) * O(2R)) + M * Mul(R+1).
// Hence we try to choose R as small as possible.
// Roughly, R >= 2*K, R >= M/2, hence R^2 >= K*M >= N.

// For example, when N1 = N2 = 1000000:
// Choosing R = 1024, M = 2048, K = 506, ceiling(N1/K) = ceiling(N2/K) = 1977,
// M >= 3953, doesn't work.
// Choosing R = 2048, M = 4096, K = 1018, ceiling(N1/K) = ceiling(N2/K) = 983,
// M >= 1965, works.
// Actually, we will also want  intDsize | K,  so that splitting into chunks
// and putting together the result can be done without shifts. So
// choose R = 2048, M = 4096, K = 992, ceiling(N1/K) = ceiling(N2/K) = 1009.
// We see that M = 2048 suffices.

// In contrast to Nussbaumer multiplication, here we can use the standard
// Karatsuba algorithm for multiplication mod p = 2^R+1. We don't have to
// recurse until N=1.


// Define this for (cheap) consistency checks.
//#define DEBUG_FFTM


// Operations modulo p = 2^R+1, each chunk represented as chlen words
// (chlen = floor(R/intDsize)+1).

static inline void assign (const uintC R, const uintC chlen,
                           const uintD* a, uintD* r)
{
	unused R;
	copy_loop_lsp(a,r,chlen);
}

// r := (a + b) mod p
static void addm (const uintC R, const uintC chlen,
                  const uintD* a, const uintD* b, uintD* r)
{
	unused R;
	// r := a+b.
	add_loop_lsp(a,b, r, chlen);
#if 0
	if (lspref(r,chlen-1) < ((uintD)1 << (R % intDsize)))
		return;
	if (lspref(r,chlen-1) == ((uintD)1 << (R % intDsize)))
		if (!DS_test_loop(r lspop (chlen-1),chlen-1,r))
			return;
	// r >= p, so subtract r := r-p.
	lspref(r,chlen-1) -= ((uintD)1 << (R % intDsize));
	dec_loop_lsp(r,chlen);
#else
	if (lspref(r,chlen-1) < 1)
		return;
	if (lspref(r,chlen-1) == 1)
		if (!DS_test_loop(r lspop (chlen-1),chlen-1,r))
			return;
	// r >= p, so subtract r := r-p.
	lspref(r,chlen-1) -= 1;
	dec_loop_lsp(r,chlen);
#endif
}

// r := (a - b) mod p
static void subm (const uintC R, const uintC chlen,
                  const uintD* a, const uintD* b, uintD* r)
{
	unused R;
	// r := a-b.
	sub_loop_lsp(a,b, r, chlen);
#if 0
	if ((sintD)lspref(r,chlen-1) >= 0)
		return;
	// r < 0, so add r := r+p.
	lspref(r,chlen-1) += ((uintD)1 << (R % intDsize));
	inc_loop_lsp(r,chlen);
#else
	if ((sintD)lspref(r,chlen-1) >= 0)
		return;
	// r < 0, so add r := r+p.
	lspref(r,chlen-1) += 1;
	inc_loop_lsp(r,chlen);
#endif
}

// r := (a << s) mod p (0 <= s < R).
// Assume that a and r don't overlap.
static void shiftleftm (const uintC R, const uintC chlen,
                        const uintD* a, uintL s, uintD* r)
{
	// Write a = 2^(R-s)*b + c, then
	// a << s = 2^R*b + (c << s) = (c << s) - b.
#if 0
	if (chlen == 1) {
		// R < intDsize.
		var uintD b = lspref(a,0) >> (R-s);
		var uintD c = lspref(a,0) & (((uintD)1 << (R-s)) - 1);
		c = c << s;
		c -= b;
		if ((sintD)c < 0)
			c += ((uintD)1 << R) + 1;
		lspref(r,0) = c;
		return;
	}
#endif
	// Here R >= intDsize, hence intDsize | R.
	if ((s % intDsize) == 0) {
		var uintP lenb = s/intDsize;
		var uintP lenc = (R-s)/intDsize;
		// chlen = 1 + lenb + lenc.
		lspref(r,lenb+lenc) = 0;
		copy_loop_lsp(a,r lspop lenb,lenc);
		copy_loop_lsp(a lspop lenc,r,lenb);
		if ((lspref(a,lenb+lenc) > 0) || neg_loop_lsp(r,lenb)) // -b gives carry?
			if (dec_loop_lsp(r lspop lenb,lenc))
				// add p = 2^R+1 to compensate with carry
				inc_loop_lsp(r,chlen);
	} else {
		var uintP lenb = floor(s,intDsize);
		var uintP lenc = floor(R-s,intDsize)+1;
		// chlen = 1 + lenb + lenc.
		s = s % intDsize;
		lspref(r,lenb+lenc) = 0;
		var uintD b0 = shiftleftcopy_loop_lsp(a,r lspop lenb,lenc,s);
		var uintD bov;
		if (lenb == 0)
			bov = b0;
		else {
			bov = shiftleftcopy_loop_lsp(a lspop lenc,r,lenb,s);
			lspref(r,0) |= b0;
		}
		bov |= lspref(a,lenb+lenc) << s;
		if (neg_loop_lsp(r,lenb))
			bov++;
		if (lspref(r,lenb) >= bov)
			lspref(r,lenb) -= bov;
		else {
			lspref(r,lenb) -= bov;
			if (dec_loop_lsp(r lspop (lenb+1),lenc-1))
				// add p = 2^R+1 to compensate with carry
				inc_loop_lsp(r,chlen);
		}
	}
}

// r := (a * b) mod p
static void mulm (const uintC R, const uintC chlen,
                  const uintD* a, const uintD* b, uintD* r)
{
	unused R;
	// The leading digits are very likely to be 0.
	var uintP a_len = chlen;
	if (lspref(a,a_len-1) == 0)
		do {
			a_len--;
		} while ((a_len > 0) && (lspref(a,a_len-1) == 0));
	if (a_len == 0) {
		clear_loop_lsp(r,chlen);
		return;
	}
	var uintP b_len = chlen;
	if (lspref(b,b_len-1) == 0)
		do {
			b_len--;
		} while ((b_len > 0) && (lspref(b,b_len-1) == 0));
	if (b_len == 0) {
		clear_loop_lsp(r,chlen);
		return;
	}
	CL_SMALL_ALLOCA_STACK;
	var uintD* tmp = cl_small_alloc_array(uintD,2*chlen);
	cl_UDS_mul(a,a_len, b,b_len, arrayLSDptr(tmp,2*chlen));
	DS_clear_loop(arrayMSDptr(tmp,2*chlen),2*chlen-(a_len+b_len),arrayLSDptr(tmp,2*chlen) lspop (a_len+b_len));
	// To divide c (0 <= c < p^2) by p = 2^R+1,
	// we set q := floor(c/2^R) and r := c - q*p = (c mod 2^R) - q.
	// If this becomes negative, set r := r + p (at most twice).
	// (This works because  floor(c/p) <= q <= floor(c/p)+2.)
	// (Actually, here, 0 <= c <= (p-1)^2, hence
	// floor(c/p) <= q <= floor(c/p)+1, so we have
	// to set r := r + p at most once!)
#if 0
	if (chlen == 1) {
		// R < intDsize.
		var uintD r0 = (arrayLSref(tmp,2,0) & (((uintD)1 << R) - 1))
		             - ((arrayLSref(tmp,2,1) << (intDsize-R)) | (arrayLSref(tmp,2,0) >> R));
		if ((sintD)r0 < 0)
			r0 += ((uintD)1 << R) + 1;
		lspref(r,0) = r0;
		return;
	}
#endif
	// Here R >= intDsize, hence intDsize | R.
	// R/intDsize = chlen-1.
	// arrayLSref(tmp,2*chlen,2*chlen-1) = 0, arrayLSref(tmp,2*chlen,2*chlen-2) <= 1.
	lspref(r,chlen-1) = 0;
	if (sub_loop_lsp(arrayLSDptr(tmp,2*chlen),arrayLSDptr(tmp,2*chlen) lspop (chlen-1),r,chlen-1) || arrayLSref(tmp,2*chlen,2*chlen-2))
		// add p = 2^R+1 to compensate with carry
		inc_loop_lsp(r,chlen);
}

// b := (a / 2) mod p
static void shiftm (const uintC R, const uintC chlen,
                    const uintD* a, uintD* b)
{
	unused R;
	shiftrightcopy_loop_msp(a lspop chlen,b lspop chlen,chlen,1,0);
	if (lspref(a,0) & 1) {
		// ((a + p) >> 1) = (a >> 1) + (p>>1) + 1.
#if 0
		if (chlen == 1)
			// R < intDsize.
			lspref(b,0) |= ((uintD)1 << (R-1));
		else
#endif
			// intDsize | R.
			lspref(b,chlen-2) |= ((uintD)1 << (intDsize-1));
		inc_loop_lsp(b,chlen);
	}
}


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

static void mulu_fftm (const uintL r, const uintC R, // R = 2^r
                       const uintL m, const uintC M, // M = 2^m
                       const uintC k, // K = intDsize*k
                       const uintD* sourceptr1, uintC len1,
                       const uintD* sourceptr2, uintC len2,
                       uintD* destptr)
// Assume:
//      ceiling(len1/k)+ceiling(len2/k)-1 <= M,
//      2*K+m <= R,
//      R = 2^r, M = 2^m, M | 2R.
//      m > 0.
{
	var const uintC chlen = floor(R,intDsize)+1; // chunk length (in words)
	CL_ALLOCA_STACK;
	var uintD* const arrX = cl_alloc_array(uintD,chlen<<m);
	var uintD* const arrY = cl_alloc_array(uintD,chlen<<m);
	#ifdef DEBUG_FFTM
	var uintD* const arrZ = cl_alloc_array(uintD,chlen<<m);
	#else
	var uintD* const arrZ = arrX; // put Z in place of X - saves memory
	#endif
	#define X(i) arrayLSDptr(&arrX[chlen*(i)],chlen)
	#define Y(i) arrayLSDptr(&arrY[chlen*(i)],chlen)
	#define Z(i) (arrayLSDptr(&arrZ[chlen*(i)],chlen))
	var uintD* tmp;
	var uintD* sum;
	var uintD* diff;
	num_stack_alloc(chlen,,tmp=);
	num_stack_alloc(chlen,,sum=);
	num_stack_alloc(chlen,,diff=);
	var bool squaring = ((sourceptr1 == sourceptr2) && (len1 == len2));
	var uintC i;
	// Initialize factors X(i) and Y(i).
	{
		var const uintD* sptr = sourceptr1;
		var uintC slen = len1;
		for (i = 0; i < M; i++) {
			var uintD* ptr = X(i);
			if (slen >= k) {
				copy_loop_lsp(sptr,ptr,k);
				clear_loop_lsp(ptr lspop k,chlen-k);
				sptr = sptr lspop k;
				slen -= k;
			} else {
				copy_loop_lsp(sptr,ptr,slen);
				clear_loop_lsp(ptr lspop slen,chlen-slen);
				i++;
				break;
			}
		}
		// X(i) := ... := X(M-1) := 0
		clear_loop_up(&arrX[chlen*i],chlen*(M-i));
	}
	if (!squaring) {
		var const uintD* sptr = sourceptr2;
		var uintC slen = len2;
		for (i = 0; i < M; i++) {
			var uintD* ptr = Y(i);
			if (slen >= k) {
				copy_loop_lsp(sptr,ptr,k);
				clear_loop_lsp(ptr lspop k,chlen-k);
				sptr = sptr lspop k;
				slen -= k;
			} else {
				copy_loop_lsp(sptr,ptr,slen);
				clear_loop_lsp(ptr lspop slen,chlen-slen);
				i++;
				break;
			}
		}
		// Y(i) := ... := Y(M-1) := 0
		clear_loop_up(&arrY[chlen*i],chlen*(M-i));
	}
	// Do an FFT of length M on X. w = 2^(2R/M) = 2^(2^(r+1-m)).
	{
		var sintL l;
		/* l = m-1 */ {
			var const uintC tmax = M>>1; // tmax = 2^(m-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Butterfly: replace (X(i1),X(i2)) by
				// (X(i1) + X(i2), X(i1) - X(i2)).
				assign(R,chlen, X(i2), tmp);
				subm(R,chlen, X(i1),tmp, X(i2));
				addm(R,chlen, X(i1),tmp, X(i1));
			}
		}
		for (l = m-2; l>=0; l--) {
			var const uintC smax = (uintC)1 << (m-1-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				// w^exp = 2^(exp << (r+1-m)).
				var uintC exp = bit_reverse(m-1-l,s) << (r-(m-1-l));
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Butterfly: replace (X(i1),X(i2)) by
					// (X(i1) + w^exp*X(i2), X(i1) - w^exp*X(i2)).
					shiftleftm(R,chlen, X(i2),exp, tmp);
					subm(R,chlen, X(i1),tmp, X(i2));
					addm(R,chlen, X(i1),tmp, X(i1));
				}
			}
		}
	}
	// Do an FFT of length M on Y. w = 2^(2R/M) = 2^(2^(r+1-m)).
	if (!squaring) {
		var sintL l;
		/* l = m-1 */ {
			var const uintC tmax = M>>1; // tmax = 2^(m-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Butterfly: replace (Y(i1),Y(i2)) by
				// (Y(i1) + Y(i2), Y(i1) - Y(i2)).
				assign(R,chlen, Y(i2), tmp);
				subm(R,chlen, Y(i1),tmp, Y(i2));
				addm(R,chlen, Y(i1),tmp, Y(i1));
			}
		}
		for (l = m-2; l>=0; l--) {
			var const uintC smax = (uintC)1 << (m-1-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				// w^exp = 2^(exp << (r+1-m)).
				var uintC exp = bit_reverse(m-1-l,s) << (r-(m-1-l));
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Butterfly: replace (Y(i1),Y(i2)) by
					// (Y(i1) + w^exp*Y(i2), Y(i1) - w^exp*Y(i2)).
					shiftleftm(R,chlen, Y(i2),exp, tmp);
					subm(R,chlen, Y(i1),tmp, Y(i2));
					addm(R,chlen, Y(i1),tmp, Y(i1));
				}
			}
		}
	}
	// Multiply the transformed vectors into Z.
	if (!squaring) {
		for (i = 0; i < M; i++)
			mulm(R,chlen, X(i),Y(i),Z(i));
	} else {
		for (i = 0; i < M; i++)
			mulm(R,chlen, X(i),X(i),Z(i));
	}
	// Undo an FFT of length M on Z. w = 2^(2R/M) = 2^(2^(r+1-m)).
	{
		var uintL l;
		for (l = 0; l < m-1; l++) {
			var const uintC smax = (uintC)1 << (m-1-l);
			var const uintC tmax = (uintC)1 << l;
			/* s = 0, exp = 0 */ {
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (Z(i1),Z(i2)) by
					// ((Z(i1)+Z(i2))/2, (Z(i1)-Z(i2))/(2*w^exp)),
					// with exp <-- 0.
					addm(R,chlen, Z(i1),Z(i2), sum);
					subm(R,chlen, Z(i1),Z(i2), diff);
					shiftm(R,chlen, sum, Z(i1));
					shiftm(R,chlen, diff, Z(i2));
				}
			}
			for (var uintC s = 1; s < smax; s++) {
				// w^exp = 2^(exp << (r+1-m)).
				var uintC exp = bit_reverse(m-1-l,s) << (r-(m-1-l));
				exp = R - exp; // negate exp (use w^-1 instead of w)
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (Z(i1),Z(i2)) by
					// ((Z(i1)+Z(i2))/2, (Z(i1)-Z(i2))/(2*w^exp)),
					// with exp <-- (M/2 - exp).
					addm(R,chlen, Z(i1),Z(i2), sum);
					subm(R,chlen, Z(i2),Z(i1), diff); // note that w^(M/2) = 2^R = -1
					shiftm(R,chlen, sum, Z(i1));
					shiftleftm(R,chlen, diff,exp-1, Z(i2));
				}
			}
		}
		/* l = m-1 */ {
			var const uintC tmax = M>>1; // tmax = 2^(m-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Inverse Butterfly: replace (Z(i1),Z(i2)) by
				// ((Z(i1)+Z(i2))/2, (Z(i1)-Z(i2))/2).
				addm(R,chlen, Z(i1),Z(i2), sum);
				subm(R,chlen, Z(i1),Z(i2), diff);
				shiftm(R,chlen, sum, Z(i1));
				shiftm(R,chlen, diff, Z(i2));
			}
		}
	}
	var uintC zchlen = 2*k + ceiling(m,intDsize);
	#ifdef DEBUG_FFTM
	// Check that every Z(i) has at most 2*K+m bits.
	{
		var uintC zerodigits = chlen - zchlen;
		for (i = 0; i < M; i++)
			if (DS_test_loop(Z(i) lspop chlen,zerodigits,Z(i) lspop zchlen))
				throw runtime_exception();
	}
	#endif
	// Put together result.
	var uintC destlen = len1+len2;
	clear_loop_lsp(destptr,destlen);
	for (i = 0; i < M; i++, destptr = destptr lspop k, destlen -= k) {
		if (zchlen <= destlen) {
			if (addto_loop_lsp(Z(i),destptr,zchlen))
				if (inc_loop_lsp(destptr lspop zchlen,destlen-zchlen))
					throw runtime_exception();
		} else {
			#ifdef DEBUG_FFTM
			if (DS_test_loop(Z(i) lspop zchlen,zchlen-destlen,Z(i) lspop destlen))
				throw runtime_exception();
			#endif
			if (addto_loop_lsp(Z(i),destptr,destlen))
				throw runtime_exception();
		}
		if (destlen <= k) {
			i++;
			break;
		}
	}
	#ifdef DEBUG_FFTM
	// Check that Z(i)..Z(M-1) are all zero.
	if (test_loop_up(&arrZ[chlen*i],chlen*(M-i)))
		throw runtime_exception();
	#endif
	#undef diff
	#undef sum
	#undef tmp
	#undef Z
	#undef Y
	#undef X
}

// The running time of mulu_fftm() is roughly
//      O(M log(M) * O(2R)) + M * R^(1+c),  where c = log3/log2 - 1 = 0.585...
// Try to minimize this given the constraints
//      ceiling(len1/k)+ceiling(len2/k)-1 <= M,
//      K = intDsize*k, 2*K+m <= R,
//      R = 2^r, M = 2^m, M | 2R.
//      m > 0.
// Necessary conditions:
//      len1+len2 <= k*(M+1),  intDsize*(len1+len2) <= K*(M+1) <= (R-1)/2 * (2*R+1) < R^2.
//      2*intDsize+1 <= R, log2_intDsize+1 < r.
// So we start with len1 <= len2,
//    r := max(log2_intDsize+2,ceiling(ceiling(log2(intDsize*2*len1))/2)), R := 2^r.
//    try
//      kmax := floor((R-(r+1))/(2*intDsize)), Kmax := intDsize*kmax,
//      m := max(1,ceiling(log2(2*ceiling(len1/kmax)-1))), M := 2^m,
//      if m > r+1 retry with r <- r+1.
//    [Now we are sure that we can at least multiply len1 and len1 digits using these
//     values of r and m, symbolically (r,m) OKFOR (len1,len1).]
//    [Normally, we will have m=r+1 or m=r.]
//    For (len1,len2), we might want to split the second integer into pieces.
//    If (r,m) OKFOR (len1,len2)
//      If (r-1,m) OKFOR (len1,ceiling(len2/2))
//        then use (r-1,m) and two pieces
//        else use (r,m) and one piece
//    else
//      q1 := number of pieces len2 needs to be splitted into to be OKFOR (r,m),
//      If m<r+1
//        q2 := number of pieces len2 needs to be splitted into to be OKFOR (r,m+1),
//        If 2*q2 <= q1
//          then use (r,m+1) and q2 pieces
//          else use (r,m) and q1 pieces
//      else m=r+1
//        q2 := number of pieces len2 needs to be splitted into to be OKFOR (r+1,m),
//        If 3*q2 <= q1
//          then use (r+1,m) and q2 pieces
//          else use (r,m) and q1 pieces

// Because we always choose r >= log2_intDsize+2, R >= 4*intDsize, so chlen >= 5.
// To avoid infinite recursion, mulu_fft_modm() must only be called with len1 > 5.

static bool okfor (uintL r, uintL m, uintC len1, uintC len2)
{
	var uintC R = (uintC)1 << r;
	var uintC M = (uintC)1 << m;
	var uintC k = floor(R-m,2*intDsize);
	return (ceiling(len1,k)+ceiling(len2,k) <= M+1);
}

static uintC numpieces (uintL r, uintL m, uintC len1, uintC len2)
{
	var uintC R = (uintC)1 << r;
	var uintC M = (uintC)1 << m;
	var uintC k = floor(R-m,2*intDsize);
	var uintC piecelen2 = (M+1-ceiling(len1,k))*k;
	#ifdef DEBUG_FFTM
	if ((sintC)piecelen2 <= 0)
		throw runtime_exception();
	#endif
	return ceiling(len2,piecelen2);
}

static void mulu_fft_modm (const uintD* sourceptr1, uintC len1,
                           const uintD* sourceptr2, uintC len2,
                           uintD* destptr)
  // Called only with 6 <= len1 <= len2.
  {
	var uint32 n;
	integerlengthC(len1-1, n=); // 2^(n-1) < len1 <= 2^n
	var uintL r;
	var uintL m;
	r = ceiling(log2_intDsize+1+n,2);
	if (r < log2_intDsize+2)
		r = log2_intDsize+2;
	retry: {
		var uintC k = floor(((uintC)1 << r) - (r+1), 2*intDsize);
		var uintC M = 2*ceiling(len1,k)-1;
		integerlengthC(M, m=);
		if (m == 0)
			m = 1;
		if (m > r+1) {
			r++;
			goto retry;
		}
	}
	#ifdef DEBUG_FFTM
	if (!(m > 0 && m <= r+1 && okfor(r,m,len1,len1)))
		throw runtime_exception();
	#endif
	if (okfor(r,m,len1,len2)) {
		if ((m <= r) && (r > log2_intDsize+2) && okfor(r-1,m,len1,ceiling(len2,2)))
			if (!(sourceptr1 == sourceptr2 && len1 == len2)) // when squaring, keep one piece
				r--;
	} else {
		var uintC q1 = numpieces(r,m,len1,len2);
		if (m <= r) {
			var uintC q2 = numpieces(r,m+1,len1,len2);
			if (2*q2 <= q1)
				m++;
		} else {
			var uintC q2 = numpieces(r+1,m,len1,len2);
			if (3*q2 <= q1)
				r++;
		}
	}
	var uintC R = (uintC)1 << r;
	var uintC M = (uintC)1 << m;
	var uintC k = floor(R-m,2*intDsize);
	var uintC piecelen2 = (M+1-ceiling(len1,k))*k;
	if (piecelen2 >= len2) {
		// One piece only.
		mulu_fftm(r,R, m,M, k, sourceptr1,len1, sourceptr2,len2, destptr);
		return;
	}
	CL_ALLOCA_STACK;
	var uintD* tmpptr;
	num_stack_alloc(len1+piecelen2,,tmpptr=);
	var uintC destlen = len1+len2;
	clear_loop_lsp(destptr,destlen);
	do {
		var uintC len2p; // length of a piece of source2
		len2p = piecelen2;
		if (len2p > len2)
			len2p = len2;
		// len2p = min(piecelen2,len2).
		var uintC destlenp = len1 + len2p;
		// destlenp = min(len1+piecelen2,destlen).
		// Use tmpptr[-destlenp..-1].
		if (len2p == 1) {
			// cheap case
			mulu_loop_lsp(lspref(sourceptr2,0),sourceptr1,tmpptr,len1);
		} else if (2*len2p < piecelen2) {
			// semi-cheap case
			cl_UDS_mul(sourceptr1,len1, sourceptr2,len2p, tmpptr);
		} else {
			mulu_fftm(r,R, m,M, k, sourceptr1,len1, sourceptr2,len2p, tmpptr);
		}
		if (addto_loop_lsp(tmpptr,destptr,destlenp))
			if (inc_loop_lsp(destptr lspop destlenp,destlen-destlenp))
				throw runtime_exception();
		// Decrement len2.
		destptr = destptr lspop len2p;
		destlen -= len2p;
		sourceptr2 = sourceptr2 lspop len2p;
		len2 -= len2p;
	} while (len2 > 0);
}
