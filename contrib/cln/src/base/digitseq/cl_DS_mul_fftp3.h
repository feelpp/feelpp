// Fast integer multiplication using FFT in a modular ring.
// Bruno Haible 5.5.1996, 30.6.1996

// FFT in the complex domain has the drawback that it needs careful round-off
// error analysis. So here we choose another field of characteristic 0: Q_p.
// Since Q_p contains exactly the (p-1)th roots of unity, we choose
// p == 1 mod N and have the Nth roots of unity (N = 2^n) in Q_p and
// even in Z_p. Actually, we compute in Z/(p^m Z).

// All operations the FFT algorithm needs is addition, subtraction,
// multiplication, multiplication by the Nth root of unity and division
// by N. Hence we can use the domain Z/(p^m Z) even if p is not a prime!

// We want to compute the convolution of N 32-bit words. The resulting
// words are < (2^32)^2 * N. To avoid computing with numbers greater than
// 32 bits, we compute in Z/pZ for three different primes p in parallel,
// i.e. we compute in the ring (Z / p1 Z) x (Z / p2 Z) x (Z / p3 Z). We choose
// p1 = 3*2^30+1, p2 = 13*2^28+1, p3 = 29*2^27+1 (or 15*2^27+1 or 17*2^27+1).
// Because of p1*p2*p3 >= (2^32)^2 * N, the chinese remainder theorem will
// faithfully combine 3 32-bit words to a word < (2^32)^2 * N.

#if !(intDsize==32)
#error "fft mod p implemented only for intDsize==32"
#endif

static const uint32 p1 = 1+(3<<30); // = 3221225473
static const uint32 p2 = 1+(13<<28); // = 3489660929
static const uint32 p3 = 1+(29<<27); // = 3892314113

typedef struct {
	uint32 w1; // remainder mod p1
	uint32 w2; // remainder mod p2
	uint32 w3; // remainder mod p3
} fftp3_word;

static const fftp3_word fftp3_roots_of_1 [27+1] =
  // roots_of_1[n] is a (2^n)th root of unity in our ring.
  // (Also roots_of_1[n-1] = roots_of_1[n]^2, but we don't need this.)
  {
    {          1,          1,          1 },
    { 3221225472, 3489660928, 3892314112 },
    { 1013946479, 1647819299,  800380159 },
    { 1031213943, 1728043888, 1502037594 },
    {  694614138,  156262243, 1093602721 },
    {  347220834,  408915340,  491290336 },
    {  680684264,  452506952,  846570852 },
    { 1109768284, 1230864251,  390870396 },
    {  602134989,   11914870, 1791906422 },
    { 1080308101,  336213294,  158126993 },
    {  381653707, 1548648704, 1432108380 },
    {  902453688,  429650884,  798472051 },
    { 1559299664,  775532293, 1877725713 },
    {  254499731,  727160889, 1192318337 },
    { 1376063215, 1557302953, 1642774092 },
    { 1284040478,  937059094, 1876917422 },
    {  336664489, 1411926644,  682311165 },
    {  894491787,  534027329,  473693773 },
    {  795860341, 1178663675, 1313928891 },
    {   23880336, 1707047452,   93147496 },
    {  790585193,  892284267, 1647947905 },
    {  877386874, 1729337527, 1233672227 },
    { 1510644826,  333000282,  296593948 },
    {  353060343,  901807544,  659274384 },
    {  716717815, 1281544649, 1457308949 },
    { 1020271667, 1713714919,  627726344 },
    {  139914905,  950720020, 1119863241 },
    {  709308748,  675675166,  538726428 }
  };

// Define this for (cheap) consistency checks.
//#define DEBUG_FFTP3

// Define this for extensive consistency checks.
//#define DEBUG_FFTP3_OPERATIONS

// Define the algorithm of the backward FFT:
// Either FORWARD (a normal FFT followed by a permutation)
// or     RECIPROOT (an FFT with reciprocal root of unity)
// or     CLEVER (an FFT with reciprocal root of unity but clever computation
//                of the reciprocals).
// Drawback of FORWARD: the permutation pass.
// Drawback of RECIPROOT: need all the powers of the root, not only half of them.
#define FORWARD   42
#define RECIPROOT 43
#define CLEVER    44
#define FFTP3_BACKWARD CLEVER

#ifdef DEBUG_FFTP3_OPERATIONS
#define check_fftp3_word(x)  if ((x.w1 >= p1) || (x.w2 >= p2) || (x.w3 >= p3)) throw runtime_exception()
#else
#define check_fftp3_word(x)
#endif

// r := 0 mod p
static inline void zerop3 (fftp3_word& r)
{
	r.w1 = 0;
	r.w2 = 0;
	r.w3 = 0;
}

// r := x mod p
static inline void setp3 (uint32 x, fftp3_word& r)
{
	if (p1 >= ((uint32)1 << 31))
		r.w1 = (x >= p1 ? x - p1 : x);
	else
		divu_3232_3232(x,p1, ,r.w1=);
	if (p2 >= ((uint32)1 << 31))
		r.w2 = (x >= p2 ? x - p2 : x);
	else
		divu_3232_3232(x,p2, ,r.w2=);
	if (p3 >= ((uint32)1 << 31))
		r.w3 = (x >= p3 ? x - p3 : x);
	else
		divu_3232_3232(x,p3, ,r.w3=);
}

// Chinese remainder theorem:
// (Z / p1 Z) x (Z / p2 Z) x (Z / p3 Z) == Z / p1*p2*p3 Z = Z / P Z.
// Return r as an integer >= 0, < p1*p2*p3, as 3-digit-sequence res.
static void combinep3 (const fftp3_word& r, uintD* resLSDptr)
{
	check_fftp3_word(r);
	// Compute e1 * r.w1 + e2 * r.w2 + e3 * r.w3 where the idempotents are
	// found as: xgcd(pi,p/pi) = 1 = ui*pi + vi*P/pi, ei = 1 - ui*pi.
	// e1 = 35002755423056150739595925972
	// e2 = 14584479687667766215746868453
	// e3 = 37919651490985126265126719818
	// Since e1+e2+e3 > 2*P, we prefer to compute with the negated
	// idempotents, their sum is < P:
	// -e1 =  8750687877798370870638831149
	// -e2 = 29168963613186755394487888668
	// -e3 =  5833791809869395345108037303
	// We will have 0 <= -e1 * r.w1 + -e2 * r.w2 + -e3 * r.w3 <
	// < -e1 * p1 + -e2 * p2 + -e3 * p3 < 2^32 * p1*p2*p3 < 2^128.
	// The sum of the products fits in 4 digits, we divide by p1*p2*p3
	// as a 3-digit sequence and finally negate the remainder.
	#if CL_DS_BIG_ENDIAN_P
	var const uintD p123 [3] = { 0x8D600002, 0x06800002, 0x78000001 };
	var const uintD e1 [3] = { 0x1C46663C, 0x647FFF9D, 0x7E66662D };
	var const uintD e2 [3] = { 0x5E40004D, 0xEDAAAB66, 0xEAAAAB1C };
	var const uintD e3 [3] = { 0x12D99977, 0xB45554FE, 0x0EEEEEB7 };
	#else
	var const uintD p123 [3] = { 0x78000001, 0x06800002, 0x8D600002 };
	var const uintD e1 [3] = { 0x7E66662D, 0x647FFF9D, 0x1C46663C };
	var const uintD e2 [3] = { 0xEAAAAB1C, 0xEDAAAB66, 0x5E40004D };
	var const uintD e3 [3] = { 0x0EEEEEB7, 0xB45554FE, 0x12D99977 };
	#endif
	var uintD sum [4];
	var uintD* const sumLSDptr = arrayLSDptr(sum,4);
	mulu_loop_lsp(r.w1,arrayLSDptr(e1,3), sumLSDptr,3);
	lspref(sumLSDptr,3) += muluadd_loop_lsp(r.w2,arrayLSDptr(e2,3), sumLSDptr,3);
	lspref(sumLSDptr,3) += muluadd_loop_lsp(r.w3,arrayLSDptr(e3,3), sumLSDptr,3);
	#if 0
	{CL_ALLOCA_STACK;
	 var DS q;
	 var DS r;
	 UDS_divide(arrayMSDptr(sum,4),4,arrayLSDptr(sum,4),
	            arrayMSDptr(p123,3),3,arrayLSDptr(p123,3),
	            &q,&r
	           );
	 ASSERT(q.len <= 1)
	 ASSERT(r.len <= 3)
	 copy_loop_lsp(r.LSDptr,arrayLSDptr(sum,4),r.len);
	 DS_clear_loop(arrayMSDptr(sum,4) mspop 1,3-r.len,arrayLSDptr(sum,4) lspop r.len);
	}
	#else
	// Division wie UDS_divide mit a_len=4, b_len=3.
	{
		var uintD q_stern;
		var uintD c1;
		#if HAVE_DD
		  divuD(highlowDD(lspref(sumLSDptr,3),lspref(sumLSDptr,2)),lspref(arrayLSDptr(p123,3),2), q_stern=,c1=);
		  { var uintDD c2 = highlowDD(c1,lspref(sumLSDptr,1));
		    var uintDD c3 = muluD(lspref(arrayLSDptr(p123,3),1),q_stern);
		    if (c3 > c2)
		      { q_stern = q_stern-1;
		        if (c3-c2 > highlowDD(lspref(arrayLSDptr(p123,3),2),lspref(arrayLSDptr(p123,3),1)))
		          { q_stern = q_stern-1; }
		  }   }
		#else
		  divuD(lspref(sumLSDptr,3),lspref(sumLSDptr,2),lspref(arrayLSDptr(p123,3),2), q_stern=,c1=);
		  { var uintD c2lo = lspref(sumLSDptr,1);
		    var uintD c3hi;
		    var uintD c3lo;
		    muluD(lspref(arrayLSDptr(p123,3),1),q_stern, c3hi=,c3lo=);
		    if ((c3hi > c1) || ((c3hi == c1) && (c3lo > c2lo)))
		      { q_stern = q_stern-1;
		        c3hi -= c1; if (c3lo < c2lo) { c3hi--; }; c3lo -= c2lo;
		        if ((c3hi > lspref(arrayLSDptr(p123,3),2)) || ((c3hi == lspref(arrayLSDptr(p123,3),2)) && (c3lo > lspref(arrayLSDptr(p123,3),1))))
		          { q_stern = q_stern-1; }
                   }   }
		#endif
		if (!(q_stern==0))
		  { var uintD carry = mulusub_loop_lsp(q_stern,arrayLSDptr(p123,3),sumLSDptr,3);
		    if (carry > lspref(sumLSDptr,3))
		      { q_stern = q_stern-1;
		        addto_loop_lsp(arrayLSDptr(p123,3),sumLSDptr,3);
		  }   }
	}
	#endif
	if (lspref(sumLSDptr,0)==0 && lspref(sumLSDptr,1)==0 && lspref(sumLSDptr,2)==0) {
		clear_loop_lsp(resLSDptr,3);
	} else {
		sub_loop_lsp(arrayLSDptr(p123,3),sumLSDptr,resLSDptr,3);
	}
}

// r := (a + b) mod p
static inline void addp3 (const fftp3_word& a, const fftp3_word& b, fftp3_word& r)
{
	var uint32 x;

	check_fftp3_word(a); check_fftp3_word(b);
	// Add single 32-bit words mod pi.
	if (((x = (a.w1 + b.w1)) < b.w1) || (x >= p1))
		x -= p1;
	r.w1 = x;
	if (((x = (a.w2 + b.w2)) < b.w2) || (x >= p2))
		x -= p2;
	r.w2 = x;
	if (((x = (a.w3 + b.w3)) < b.w3) || (x >= p3))
		x -= p3;
	r.w3 = x;
	check_fftp3_word(r);
}

// r := (a - b) mod p
static inline void subp3 (const fftp3_word& a, const fftp3_word& b, fftp3_word& r)
{
	check_fftp3_word(a); check_fftp3_word(b);
	// Subtract single 32-bit words mod pi.
	r.w1 = (a.w1 < b.w1 ? a.w1-b.w1+p1 : a.w1-b.w1);
	r.w2 = (a.w2 < b.w2 ? a.w2-b.w2+p2 : a.w2-b.w2);
	r.w3 = (a.w3 < b.w3 ? a.w3-b.w3+p3 : a.w3-b.w3);
	check_fftp3_word(r);
}

// r := (a * b) mod p
static void mulp3 (const fftp3_word& a, const fftp3_word& b, fftp3_word& res)
{
	check_fftp3_word(a); check_fftp3_word(b);
	// To divide c (0 <= c < p^2) by p = m*2^n+1,
	// we set q := floor(floor(c/2^n)/m) and
	// r := c - q*p = (floor(c/2^n) mod m)*2^n + (c mod 2^n) - q.
	// If this becomes negative, set r := r + p (at most twice).
	// (This works because  floor(c/p) <= q <= floor(c/p)+2.)
	// (Actually, here, 0 <= c <= (p-1)^2, hence
	// floor(c/p) <= q <= floor(c/p)+1, so we have
	// to set r := r + p at most once!)
	#if 1
	#define mul_mod_p(aw,bw,result_zuweisung,p,n,m)  \
	{								\
		var uint32 hi;						\
		var uint32 lo;						\
		mulu32(aw,bw, hi=,lo=);					\
		divu_6432_3232(hi,lo,p, ,result_zuweisung);		\
	}
	#else
	#define mul_mod_p(aw,bw,result_zuweisung,p,n,m)  \
	{								\
		var uint32 hi;						\
		var uint32 lo;						\
		mulu32(aw,bw, hi=,lo=);					\
		var uint32 q;						\
		var uint32 r;						\
		divu_6432_3232(hi>>n,(hi<<(32-n))|(lo>>n), m, q=,r=);	\
		r = (r << n) | (lo & (((uint32)1<<n)-1));		\
		if (r >= q) { r = r-q; } else { r = r-q+p; }		\
		result_zuweisung r;					\
	}
	#endif
	// p1 = 3*2^30+1, n = 30, m = 3
	mul_mod_p(a.w1,b.w1,res.w1=,p1,30,3);
	// p2 = 13*2^28+1, n = 28, m = 13
	mul_mod_p(a.w2,b.w2,res.w2=,p2,28,13);
	// p3 = 29*2^27+1, n = 27, m = 29
	mul_mod_p(a.w3,b.w3,res.w3=,p3,27,29);
	#undef mul_mod_p
	check_fftp3_word(res);
}
#ifdef DEBUG_FFTP3_OPERATIONS
static void mulp3_doublecheck (const fftp3_word& a, const fftp3_word& b, fftp3_word& r)
{
	fftp3_word zero, ma, mb, or;
	zerop3(zero);
	subp3(zero,a, ma);
	subp3(zero,b, mb);
	mulp3(ma,mb, or);
	mulp3(a,b, r);
	if (!((r.w1 == or.w1) && (r.w2 == or.w2) && (r.w3 == or.w3)))
		throw runtime_exception();
}
#define mulp3 mulp3_doublecheck
#endif /* DEBUG_FFTP3_OPERATIONS */

// b := (a / 2) mod p
static inline void shiftp3 (const fftp3_word& a, fftp3_word& b)
{
	check_fftp3_word(a);
	b.w1 = (a.w1 & 1 ? (a.w1 >> 1) + (p1 >> 1) + 1 : (a.w1 >> 1));
	b.w2 = (a.w2 & 1 ? (a.w2 >> 1) + (p2 >> 1) + 1 : (a.w2 >> 1));
	b.w3 = (a.w3 & 1 ? (a.w3 >> 1) + (p3 >> 1) + 1 : (a.w3 >> 1));
	check_fftp3_word(b);
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

// Compute an convolution mod p using FFT: z[0..N-1] := x[0..N-1] * y[0..N-1].
static void fftp3_convolution (const uintL n, const uintC N, // N = 2^n
                               fftp3_word * x, // N words
                               fftp3_word * y, // N words
                               fftp3_word * z  // N words result
                              )
{
	CL_ALLOCA_STACK;
	#if (FFTP3_BACKWARD == RECIPROOT) || defined(DEBUG_FFTP3)
	var fftp3_word* const w = cl_alloc_array(fftp3_word,N);
	#else
	var fftp3_word* const w = cl_alloc_array(fftp3_word,(N>>1)+1);
	#endif
	var uintC i;
	// Initialize w[i] to w^i, w a primitive N-th root of unity.
	w[0] = fftp3_roots_of_1[0];
	w[1] = fftp3_roots_of_1[n];
	#if (FFTP3_BACKWARD == RECIPROOT) || defined(DEBUG_FFTP3)
	for (i = 2; i < N; i++)
		mulp3(w[i-1],fftp3_roots_of_1[n], w[i]);
	#else // need only half of the roots
	for (i = 2; i < N>>1; i++)
		mulp3(w[i-1],fftp3_roots_of_1[n], w[i]);
	#endif
	#ifdef DEBUG_FFTP3
	// Check that w is really a primitive N-th root of unity.
	{
		var fftp3_word w_N;
		mulp3(w[N-1],fftp3_roots_of_1[n], w_N);
		if (!(w_N.w1 == 1 && w_N.w2 == 1 && w_N.w3 == 1))
			throw runtime_exception();
		w_N = w[N>>1];
		if (!(w_N.w1 == p1-1 && w_N.w2 == p2-1 && w_N.w3 == p3-1))
			throw runtime_exception();
	}
	#endif
	var bool squaring = (x == y);
	// Do an FFT of length N on x.
	{
		var sintL l;
		/* l = n-1 */ {
			var const uintC tmax = N>>1; // tmax = 2^(n-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Butterfly: replace (x(i1),x(i2)) by
				// (x(i1) + x(i2), x(i1) - x(i2)).
				var fftp3_word tmp;
				tmp = x[i2];
				subp3(x[i1],tmp, x[i2]);
				addp3(x[i1],tmp, x[i1]);
			}
		}
		for (l = n-2; l>=0; l--) {
			var const uintC smax = (uintC)1 << (n-1-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Butterfly: replace (x(i1),x(i2)) by
					// (x(i1) + w^exp*x(i2), x(i1) - w^exp*x(i2)).
					var fftp3_word tmp;
					mulp3(x[i2],w[exp], tmp);
					subp3(x[i1],tmp, x[i2]);
					addp3(x[i1],tmp, x[i1]);
				}
			}
		}
	}
	// Do an FFT of length N on y.
	if (!squaring) {
		var sintL l;
		/* l = n-1 */ {
			var uintC const tmax = N>>1; // tmax = 2^(n-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Butterfly: replace (y(i1),y(i2)) by
				// (y(i1) + y(i2), y(i1) - y(i2)).
				var fftp3_word tmp;
				tmp = y[i2];
				subp3(y[i1],tmp, y[i2]);
				addp3(y[i1],tmp, y[i1]);
			}
		}
		for (l = n-2; l>=0; l--) {
			var const uintC smax = (uintC)1 << (n-1-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Butterfly: replace (y(i1),y(i2)) by
					// (y(i1) + w^exp*y(i2), y(i1) - w^exp*y(i2)).
					var fftp3_word tmp;
					mulp3(y[i2],w[exp], tmp);
					subp3(y[i1],tmp, y[i2]);
					addp3(y[i1],tmp, y[i1]);
				}
			}
		}
	}
	// Multiply the transformed vectors into z.
	for (i = 0; i < N; i++)
		mulp3(x[i],y[i], z[i]);
	// Undo an FFT of length N on z.
	{
		var uintL l;
		for (l = 0; l < n-1; l++) {
			var const uintC smax = (uintC)1 << (n-1-l);
			var const uintC tmax = (uintC)1 << l;
			#if FFTP3_BACKWARD != CLEVER
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				#if FFTP3_BACKWARD == RECIPROOT
				if (exp > 0)
					exp = N - exp; // negate exp (use w^-1 instead of w)
				#endif
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (z(i1),z(i2)) by
					// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/(2*w^exp)).
					var fftp3_word sum;
					var fftp3_word diff;
					addp3(z[i1],z[i2], sum);
					subp3(z[i1],z[i2], diff);
					shiftp3(sum, z[i1]);
					mulp3(diff,w[exp], diff); shiftp3(diff, z[i2]);
				}
			}
			#else // FFTP3_BACKWARD == CLEVER: clever handling of negative exponents
			/* s = 0, exp = 0 */ {
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (z(i1),z(i2)) by
					// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/(2*w^exp)),
					// with exp <-- 0.
					var fftp3_word sum;
					var fftp3_word diff;
					addp3(z[i1],z[i2], sum);
					subp3(z[i1],z[i2], diff);
					shiftp3(sum, z[i1]);
					shiftp3(diff, z[i2]);
				}
			}
			for (var uintC s = 1; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				exp = (N>>1) - exp; // negate exp (use w^-1 instead of w)
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (z(i1),z(i2)) by
					// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/(2*w^exp)),
					// with exp <-- (N/2 - exp).
					var fftp3_word sum;
					var fftp3_word diff;
					addp3(z[i1],z[i2], sum);
					subp3(z[i2],z[i1], diff); // note that w^(N/2) = -1
					shiftp3(sum, z[i1]);
					mulp3(diff,w[exp], diff); shiftp3(diff, z[i2]);
				}
			}
			#endif
		}
		/* l = n-1 */ {
			var const uintC tmax = N>>1; // tmax = 2^(n-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Inverse Butterfly: replace (z(i1),z(i2)) by
				// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/2).
				var fftp3_word sum;
				var fftp3_word diff;
				addp3(z[i1],z[i2], sum);
				subp3(z[i1],z[i2], diff);
				shiftp3(sum, z[i1]);
				shiftp3(diff, z[i2]);
			}
		}
	}
	#if FFTP3_BACKWARD == FORWARD
	// Swap z[i] and z[N-i] for 0 < i < N/2.
	for (i = (N>>1)-1; i > 0; i--) {
		var fftp3_word tmp = z[i];
		z[i] = z[N-i];
		z[N-i] = tmp;
	}
	#endif
}

static void mulu_fft_modp3 (const uintD* sourceptr1, uintC len1,
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
	var fftp3_word* const x = cl_alloc_array(fftp3_word,N);
	var fftp3_word* const y = cl_alloc_array(fftp3_word,N);
	#ifdef DEBUG_FFTP3
	var fftp3_word* const z = cl_alloc_array(fftp3_word,N);
	#else
	var fftp3_word* const z = x; // put z in place of x - saves memory
	#endif
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
				for (i = 0; i < len1; i++)
					setp3(lspref(sourceptr1,i), x[i]);
				for (i = len1; i < N; i++)
					zerop3(x[i]);
			}
			// Fill factor y.
			if (!squaring) {
				for (i = 0; i < len2p; i++)
					setp3(lspref(sourceptr2,i), y[i]);
				for (i = len2p; i < N; i++)
					zerop3(y[i]);
			}
			// Multiply.
			if (!squaring)
				fftp3_convolution(n,N, &x[0], &y[0], &z[0]);
			else
				fftp3_convolution(n,N, &x[0], &x[0], &z[0]);
			// Add result to destptr[-destlen..-1]:
			{
				var uintD* ptr = destptr;
				// ac2|ac1|ac0 are an accumulator.
				var uint32 ac0 = 0;
				var uint32 ac1 = 0;
				var uint32 ac2 = 0;
				var uint32 tmp;
				for (i = 0; i < destlenp; i++) {
					// Convert z[i] to a 3-digit number.
					var uintD z_i[3];
					combinep3(z[i],arrayLSDptr(z_i,3));
					#ifdef DEBUG_FFTP3
					if (!(arrayLSref(z_i,3,2) < N))
						throw runtime_exception();
					#endif
					// Add z[i] to the accumulator.
					tmp = arrayLSref(z_i,3,0);
					if ((ac0 += tmp) < tmp) {
						if (++ac1 == 0)
							++ac2;
					}
					tmp = arrayLSref(z_i,3,1);
					if ((ac1 += tmp) < tmp)
						++ac2;
					tmp = arrayLSref(z_i,3,2);
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
			#ifdef DEBUG_FFTP3
			// If destlenp < N, check that the remaining z[i] are 0.
			for (i = destlenp; i < N; i++)
				if (z[i].w1 > 0 || z[i].w2 > 0 || z[i].w3 > 0)
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
