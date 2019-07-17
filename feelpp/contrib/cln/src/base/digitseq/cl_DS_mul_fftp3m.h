// Fast integer multiplication using FFT in a modular ring.
// Bruno Haible 5.5.1996, 30.6.1996, 20.8.1996

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
// p1 = 3*2^30+1, p2 = 15*2^27+1, p3 = 7*2^26+1.
// Because of p1*p2*p3 >= 2^91 >= (2^32)^2 * N, the chinese remainder theorem
// will faithfully combine 3 32-bit words to a word < (2^32)^2 * N.

// Furthermore we use Montgomery's modular multiplication trick
// [Peter L. Montgomery: Modular multiplication without trial division,
//  Mathematics of Computation 44 (1985), 519-521.]
//
// Assume we want to compute modulo M, M odd. V and N will be chosen
// so that V*N==1 mod M and that (a,b) --> a*b*V mod M can be more easily
// computed than (a,b) --> a*b mod M. Then, we have a ring isomorphism
//   (Z/MZ, +, * mod M)  \isomorph  (Z/MZ, +, (a,b) --> a*b*V mod M)
//   x mod M             -------->  x*N mod M
// It is thus preferrable to use x*N mod M as a "representation" of x mod M,
// especially for computations which involve at least several multiplications.
//
// The precise algorithm to compute a*b*V mod M, given a and b, and the choice
// of N and V depend on M and on the hardware. The general idea is this:
// Choose N = 2^n, so that division by N is easy. Recall that V == N^-1 mod M.
// 1. Given a and b as m-bit numbers (M <= 2^m), compute a*b in full
//    precision.
// 2. Write a*b = c*N+d, i.e. split it into components c and d.
// 3. Now a*b*V = c*N*V+d*V == c+d*V mod M.
// 4. Instead of computing d*V mod M
//    a. by full multiplication and then division mod M, or
//    b. by left shifts: repeated application of
//          x := 2*x+(0 or 1); if (x >= M) { x := x-M; }
//    we compute
//    c. by right shifts (recall that d*V == d*2^-n mod M): repeated application
//       of   if (x odd) { x := (x+M)/2; } else { x := x/2; }
// Usually one will choose N = 2^m, so that c and d have both m bits.
// Several variations are possible: In step 4 one can implement the right
// shifts in hardware. Or (for example when N = 2^160 and working on a
// 32-bit machine) one can do 32 shift steps at the same time:
// Choose M' == M^-1 mod 2^32 and compute n/32 times
//       x := (x - ((x mod 2^32) * M' mod 2^32) * M) / 2^32.
//
// Here, we deal with moduli M = p_i = j*2^k+1. These form of primes comes
// in because we need 2^n-th roots of unity mod M. But is also comes handy
// for Montgomery multiplication: Instead of choosing N = 2^32 (which makes
// up for very easy splitting in step 1) and V = j^2*2^(2*k-32), we better
// choose N = 2^k and V = -j. The algorithm now goes like this (recall that
// M is an m-bit number and j is an (m-k)-bit number):
// 1. Compute a*b in full precision, as a 2*m <= 64 bit number.
// 2. Split a*b = c*N+d, with c an (2m-k)-bit number and d an k-bit number.
// 3. a*b*V == c+d*V mod M.
// 4. Compute c mod M by splitting off the leading (m-k+1) bits of c and
//    using table lookup; the remainder (c mod 2^(m-1)) is already reduced
//    mod M.
//    Compute d*|V| the standard way; |V| has only few bits. d*|V| is
//    already reduced mod M, because d*|V| < j*2^k < M.

// In order to get best performance, we carefully choose the primes so that
// a. the table of size 2^(m-k+1) doesn't get too large,
// b. multiplication by V is easy.
// Here is a list of the interesting primes < 2^32:
//
//                       U*M+V*N = 1
//   prime       bits    N=2^n  V    U    2m<=n+32 ?
//     M           m       n
//
//   3*2^30+1     32      30    -3   1        n       (*)
//
//  13*2^28+1     32      28   -13   1        n
//
//  15*2^27+1     31      27   -15   1        n       (*)
//  17*2^27+1     32      27   -17   1        n
//  29*2^27+1     32      27   -29   1        n
//
//   7*2^26+1     29      26    -7   1        y       (*)
//  27*2^26+1     31      26   -27   1        n
//  37*2^26+1     32      26   -37   1        n
//  43*2^26+1     32      26   -43   1        n
//
//   5*2^25+1     28      25    -5   1        y
//  33*2^25+1     31      25   -33   1        n
//  51*2^25+1     31      25   -51   1        n
//  63*2^25+1     31      25   -63   1        n
//  81*2^25+1     32      25   -81   1        n
// 125*2^25+1     32      25  -125   1        n
//
//  45*2^24+1     30      24   -45   1        n
//  73*2^24+1     31      24   -73   1        n
// 127*2^24+1     31      24  -127   1        n
// 151*2^24+1     32      24  -151   1        n
// 157*2^24+1     32      24  -157   1        n
// 171*2^24+1     32      24  -171   1        n
// 193*2^24+1     32      24  -193   1        n
// 235*2^24+1     32      24  -235   1        n
// 243*2^24+1     32      24  -243   1        n
//
//  45*2^23+1     29      23   -45   1        n
// ...
//
// The inequality 2m<=n+32 would mean that c fits in a 32-bit word, but that's
// actually irrelevant because we can fetch the most significant bits of c
// before actually computing c.
// We choose the primes marked with an asterisk.


#if !(intDsize==32)
#error "fft mod p implemented only for intDsize==32"
#endif

// Avoid clash with fftp3
#define p1 fftp3m_p1
#define p2 fftp3m_p2
#define p3 fftp3m_p3
#define n1 fftp3m_n1
#define n2 fftp3m_n2
#define n3 fftp3m_n3

static const uint32 p1 = 1+(3<<30); // = 3221225473
static const uint32 p2 = 1+(15<<27); // = 2013265921
static const uint32 p3 = 1+(7<<26); // = 469762049
static const uint32 n1 = 30; // Montgomery: represent x mod p1 as x*2^n1 mod p1
static const uint32 n2 = 27; // Montgomery: represent x mod p2 as x*2^n2 mod p2
static const uint32 n3 = 26; // Montgomery: represent x mod p3 as x*2^n3 mod p3

typedef struct {
	uint32 w1; // remainder mod p1
	uint32 w2; // remainder mod p2
	uint32 w3; // remainder mod p3
} fftp3m_word;

static const fftp3m_word fftp3m_roots_of_1 [26+1] =
  // roots_of_1[n] is a (2^n)th root of unity in our ring.
  // (Also roots_of_1[n-1] = roots_of_1[n]^2, but we don't need this.)
  {
    #if 0 // in standard representation
    {          1,          1,          1 },
    { 3221225472, 2013265920,  469762048 },
    { 1013946479,  284861408,   19610091 },
    { 1031213943,  211723194,   26623616 },
    {  694614138,   78945800,  111570435 },
    {  347220834,  772607190,  135956445 },
    {  680684264,  288289890,  181505383 },
    { 1109768284,  112574482,  145518049 },
    {  602134989,  928726468,  109721424 },
    { 1080308101,  875419223,    2847903 },
    {  381653707,  510575142,  110273149 },
    {  902453688,  193023072,   65701394 },
    { 1559299664,  313561437,  181642641 },
    {  254499731,  121307056,   82315502 },
    { 1376063215,   20899142,  142137197 },
    { 1284040478,  956809618,  207661045 },
    {  336664489,  317295870,  194405005 },
    {  894491787,  785393806,    2821902 },
    {  795860341,  738526384,  230963948 },
    {   23880336,  956561758,   59211404 },
    {  790585193,  352904935,   95374542 },
    {  877386874,  836313293,  153165757 },
    { 1510644826,  971592443,   74027009 },
    {  353060343,  692611595,   24417505 },
    {  716717815,  791167605,   26032760 },
    { 1020271667,  751686895,  150976424 },
    {  139914905,  477826617,   71902965 }
    #else // in Montgomery representation
    { 1073741824,  134217728,   67108864 },
    { 2147483649, 1879048193,  402653185 },
    { 1809501489, 1054751064,  265634015 },
    { 2877487492, 1193844673,  331740947 },
    { 2989687427,  665825587,  252496823 },
    { 3105485195, 1961758775,  114795379 },
    { 1920588894, 1994046595,  175397252 },
    {  703819063,  932019131,  314756028 },
    { 3020513810, 1682915367,   51434375 },
    {  713639124, 1015380543,  133810885 },
    {  946523922, 1576574394,  454008742 },
    { 2920407577, 1597744532,  191940679 },
    { 1627717094, 1589708641,  309595372 },
    { 2062650405,  126130591,  189567235 },
    {  615054086,  267042180,  382347871 },
    { 1719470156, 1681043157,  238769593 },
    {  961520328, 1992112863,  240663313 },
    { 2923061544,   81858141,  402250056 },
    {  808455044,  487635820,  302549471 },
    { 3213265361, 1681059681,  461303277 },
    { 1883955251, 1318650285,  254810522 },
    {  781279533, 1017987605,  179445770 },
    {  570193549, 1008968995,  459186762 },
    { 3103538692,  624914534,  466273834 },
    {  834835886, 1960521414,  331825355 },
    { 1807393093, 1292064821,  246867396 },
    { 2100845347, 1578757629,   56837012 }
    #endif
  };

// Define this for (cheap) consistency checks.
//#define DEBUG_FFTP3M

// Define this for extensive consistency checks.
//#define DEBUG_FFTP3M_OPERATIONS

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
#define FFTP3M_BACKWARD CLEVER

#ifdef DEBUG_FFTP3M_OPERATIONS
#define check_fftp3m_word(x)  if ((x.w1 >= p1) || (x.w2 >= p2) || (x.w3 >= p3)) throw runtime_exception()
#else
#define check_fftp3m_word(x)
#endif

// r := 0 mod p
static inline void zerop3m (fftp3m_word& r)
{
	r.w1 = 0;
	r.w2 = 0;
	r.w3 = 0;
}

// r := x mod p
static inline void setp3m (uint32 x, fftp3m_word& r)
{
	var uint32 hi;
	var uint32 lo;
	hi = x >> (32-n1); lo = x << n1; divu_6432_3232(hi,lo,p1, ,r.w1=);
	hi = x >> (32-n2); lo = x << n2; divu_6432_3232(hi,lo,p2, ,r.w2=);
	hi = x >> (32-n3); lo = x << n3; divu_6432_3232(hi,lo,p3, ,r.w3=);
}

// Chinese remainder theorem:
// (Z / p1 Z) x (Z / p2 Z) x (Z / p3 Z) == Z / p1*p2*p3 Z = Z / P Z.
// Return r as an integer >= 0, < p1*p2*p3, as 3-digit-sequence res.
// This routine also does the "de-Montgomerizing".
static void combinep3m (const fftp3m_word& r, uintD* resLSDptr)
{
	check_fftp3m_word(r);
	// Compute e1 * v1 * r.w1 + e2 * v2 * r.w2 + e3 * v3 * r.w3 where
	// vi == 2^-ni mod pi, and the idempotents ei are found as:
	// xgcd(pi,p/pi) = 1 = ui*pi + vi*P/pi, ei = 1 - ui*pi.
	// e1 = 1709008312966733882383995583
	// e2 = 2781580629833601225216537109
	// e3 = 1602397205945693664242711343
	// e1*v1 = 965961209845827124691257285
	// e2*v2 = 927193593718183024654651603
	// e3*v3 = 969191855872201893987508667
	// We will have 0 <= e1*v1 * r.w1 + e2*v2 * r.w2 + e3*v3 * r.w3 <
	// < e1*v1 * p1 + e2*v2 * p2 + e3*v3 * p3 < 3 * 2^32 * p1*p2*p3 < 2^128.
	// The sum of the products fits in 4 digits, we divide by p1*p2*p3
	// as a 3-digit sequence, thus getting the remainder.
	#if 0
	#if CL_DS_BIG_ENDIAN_P
	var const uintD p123 [3] = { 0x09D80000, 0x7C200001, 0x54000001 };
	var const uintD e1v1 [3] = { 0x031F063E, 0x1CD1F37E, 0x20E0C7C5 };
	var const uintD e2v2 [3] = { 0x02FEF4E1, 0x6E62C875, 0x788590D3 };
	var const uintD e3v3 [3] = { 0x0321B25B, 0xC8DB371B, 0xF0E861BB };
	#else
	var const uintD p123 [3] = { 0x54000001, 0x7C200001, 0x09D80000 };
	var const uintD e1v1 [3] = { 0x20E0C7C5, 0x1CD1F37E, 0x031F063E };
	var const uintD e2v2 [3] = { 0x788590D3, 0x6E62C875, 0x02FEF4E1 };
	var const uintD e3v3 [3] = { 0xF0E861BB, 0xC8DB371B, 0x0321B25B };
	#endif
	#else
	// The final division step requires a shift left by 4 bits in order
	// to normalize p1*p2*p3. We combine this shift left with the
	// multiplications. Note that since e1v1 + e2v2 + e3v3 < p1*p2*p3,
	// there is no risk of overflow.
	#if CL_DS_BIG_ENDIAN_P
	var const uintD p123 [3] = { 0x9D800007, 0xC2000015, 0x40000010 };
	var const uintD e1v1 [3] = { 0x31F063E1, 0xCD1F37E2, 0x0E0C7C50 };
	var const uintD e2v2 [3] = { 0x2FEF4E16, 0xE62C8757, 0x88590D30 };
	var const uintD e3v3 [3] = { 0x321B25BC, 0x8DB371BF, 0x0E861BB0 };
	#else
	var const uintD p123 [3] = { 0x40000010, 0xC2000015, 0x9D800007 };
	var const uintD e1v1 [3] = { 0x0E0C7C50, 0xCD1F37E2, 0x31F063E1 };
	var const uintD e2v2 [3] = { 0x88590D30, 0xE62C8757, 0x2FEF4E16 };
	var const uintD e3v3 [3] = { 0x0E861BB0, 0x8DB371BF, 0x321B25BC };
	#endif
	#endif
	var uintD sum [4];
	var uintD* const sumLSDptr = arrayLSDptr(sum,4);
	mulu_loop_lsp(r.w1,arrayLSDptr(e1v1,3), sumLSDptr,3);
	lspref(sumLSDptr,3) += muluadd_loop_lsp(r.w2,arrayLSDptr(e2v2,3), sumLSDptr,3);
	lspref(sumLSDptr,3) += muluadd_loop_lsp(r.w3,arrayLSDptr(e3v3,3), sumLSDptr,3);
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
	#ifdef DEBUG_FFTP3M_OPERATIONS
	if (compare_loop_msp(sumLSDptr lspop 3,arrayMSDptr(p123,3),3) >= 0)
		throw runtime_exception();
	#endif
	// Renormalize the division's remainder: shift right by 4 bits.
	shiftrightcopy_loop_msp(sumLSDptr lspop 3,resLSDptr lspop 3,3,4,0);
}

// r := (a + b) mod p
static inline void addp3m (const fftp3m_word& a, const fftp3m_word& b, fftp3m_word& r)
{
	var uint32 x;

	check_fftp3m_word(a); check_fftp3m_word(b);
	// Add single 32-bit words mod pi.
	if (((x = (a.w1 + b.w1)) < b.w1) || (x >= p1))
		x -= p1;
	r.w1 = x;
	if ((x = (a.w2 + b.w2)) >= p2) // x doesn't overflow since p2 <= 2^31
		x -= p2;
	r.w2 = x;
	if ((x = (a.w3 + b.w3)) >= p3) // x doesn't overflow since p3 <= 2^31
		x -= p3;
	r.w3 = x;
	check_fftp3m_word(r);
}

// r := (a - b) mod p
static inline void subp3m (const fftp3m_word& a, const fftp3m_word& b, fftp3m_word& r)
{
	check_fftp3m_word(a); check_fftp3m_word(b);
	// Subtract single 32-bit words mod pi.
	r.w1 = (a.w1 < b.w1 ? a.w1-b.w1+p1 : a.w1-b.w1);
	r.w2 = (a.w2 < b.w2 ? a.w2-b.w2+p2 : a.w2-b.w2);
	r.w3 = (a.w3 < b.w3 ? a.w3-b.w3+p3 : a.w3-b.w3);
	check_fftp3m_word(r);
}

// r := (a * b) mod p
static void mulp3m (const fftp3m_word& a, const fftp3m_word& b, fftp3m_word& res)
{
	check_fftp3m_word(a); check_fftp3m_word(b);
	// Multiplication à la Montgomery:
	#define mul_mod_p(aw,bw,result_zuweisung,p,m,n,j,js,table)  \
	{	/* table[i] == i*2^(m-1) mod p for 0 <= i < 2^(m-n+1) */\
		var uint32 hi;						\
		var uint32 lo;						\
		mulu32(aw,bw, hi=,lo=);					\
		/* hi has 2m-32 bits */					\
		var const int l = (m-1)-(32-n);				\
		var uint32 r = table[hi>>l];				\
		hi = ((hi << (32-l)) >> (n-l)) | (lo >> n);		\
		/* hi = c mod 2^(m-1), has m-1 bits */			\
		lo = lo & (bit(n)-1);					\
		/* lo = d, has n bits */				\
		lo = (lo << js) - lo;					\
		/* lo = d*|V|, has m bits */				\
		/* Finally compute (r + hi - lo) mod p. */		\
		if (m < 32) {						\
			r += hi;					\
			if (r >= p)					\
				{ r = r - p; }				\
		} else {						\
			if (((r += hi) < hi) || (r >= p))		\
				{ r = r - p; }				\
		}							\
		r = (r < lo ? r-lo+p : r-lo);				\
		/* ifdef DEBUG_FFTP3M_OPERATIONS *			\
		var uint32 tmp;						\
		mulu32(aw,bw, hi=,lo=);					\
		divu_6432_3232(hi,lo,p, ,tmp=);				\
		mulu32(tmp,j, hi=, lo=);				\
		divu_6432_3232(hi,lo,p, ,tmp=);				\
		if (tmp != 0) { tmp = p-tmp; }				\
		if (tmp != r)						\
			throw runtime_exception();					\
		 * endif DEBUG_FFTP3M_OPERATIONS */			\
		result_zuweisung r;					\
	}
	// p1 = 3*2^30+1, n1 = 30, j1 = 3 = 2^2-1
	static uint32 table1 [8] =
	  {          0, 2147483648, 1073741823, 3221225471,
	    2147483646, 1073741821, 3221225469, 2147483644
	  };
	mul_mod_p(a.w1,b.w1,res.w1=,p1,32,30,3,2,table1);
	// p2 = 15*2^27+1, n2 = 27, j2 = 15 = 2^4-1
	static uint32 table2 [32] =
	  {          0, 1073741824,  134217727, 1207959551,
	     268435454, 1342177278,  402653181, 1476395005,
	     536870908, 1610612732,  671088635, 1744830459,
	     805306362, 1879048186,  939524089, 2013265913,
	    1073741816,  134217719, 1207959543,  268435446,
	    1342177270,  402653173, 1476394997,  536870900,
	    1610612724,  671088627, 1744830451,  805306354,
	    1879048178,  939524081, 2013265905, 1073741808
	  };
	mul_mod_p(a.w2,b.w2,res.w2=,p2,31,27,15,4,table2);
	// p3 = 7*2^26+1, n3 = 26, j3 = 7 = 2^3-1
	static uint32 table3 [16] =
	  {          0,  268435456,   67108863,  335544319,
	     134217726,  402653182,  201326589,  469762045,
	     268435452,   67108859,  335544315,  134217722,
	     402653178,  201326585,  469762041,  268435448
	  };
	mul_mod_p(a.w3,b.w3,res.w3=,p3,29,26,7,3,table3);
	#undef mul_mod_p
	check_fftp3m_word(res);
}
#ifdef DEBUG_FFTP3M_OPERATIONS
static void mulp3m_doublecheck (const fftp3m_word& a, const fftp3m_word& b, fftp3m_word& r)
{
	fftp3m_word zero, ma, mb, or;
	zerop3m(zero);
	subp3m(zero,a, ma);
	subp3m(zero,b, mb);
	mulp3m(ma,mb, or);
	mulp3m(a,b, r);
	if (!((r.w1 == or.w1) && (r.w2 == or.w2) && (r.w3 == or.w3)))
		throw runtime_exception();
}
#define mulp3m mulp3m_doublecheck
#endif /* DEBUG_FFTP3M_OPERATIONS */

// b := (a / 2) mod p
static inline void shiftp3m (const fftp3m_word& a, fftp3m_word& b)
{
	check_fftp3m_word(a);
	b.w1 = (a.w1 & 1 ? (a.w1 >> 1) + (p1 >> 1) + 1 : (a.w1 >> 1));
	b.w2 = (a.w2 & 1 ? (a.w2 >> 1) + (p2 >> 1) + 1 : (a.w2 >> 1));
	b.w3 = (a.w3 & 1 ? (a.w3 >> 1) + (p3 >> 1) + 1 : (a.w3 >> 1));
	check_fftp3m_word(b);
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
static void fftp3m_convolution (const uintL n, const uintC N, // N = 2^n
                                fftp3m_word * x, // N words
                                fftp3m_word * y, // N words
                                fftp3m_word * z  // N words result
                               )
{
	CL_ALLOCA_STACK;
	#if (FFTP3M_BACKWARD == RECIPROOT) || defined(DEBUG_FFTP3M)
	var fftp3m_word* const w = cl_alloc_array(fftp3m_word,N);
	#else
	var fftp3m_word* const w = cl_alloc_array(fftp3m_word,(N>>1)+1);
	#endif
	var uintC i;
	// Initialize w[i] to w^i, w a primitive N-th root of unity.
	w[0] = fftp3m_roots_of_1[0];
	w[1] = fftp3m_roots_of_1[n];
	#if (FFTP3M_BACKWARD == RECIPROOT) || defined(DEBUG_FFTP3M)
	for (i = 2; i < N; i++)
		mulp3m(w[i-1],fftp3m_roots_of_1[n], w[i]);
	#else // need only half of the roots
	for (i = 2; i < N>>1; i++)
		mulp3m(w[i-1],fftp3m_roots_of_1[n], w[i]);
	#endif
	#ifdef DEBUG_FFTP3M
	// Check that w is really a primitive N-th root of unity.
	{
		var fftp3m_word w_N;
		mulp3m(w[N-1],fftp3m_roots_of_1[n], w_N);
		if (!(   w_N.w1 == (uint32)1<<n1
		      && w_N.w2 == (uint32)1<<n2
		      && w_N.w3 == (uint32)1<<n3))
			throw runtime_exception();
		w_N = w[N>>1];
		if (!(   w_N.w1 == p1-((uint32)1<<n1)
		      && w_N.w2 == p2-((uint32)1<<n2)
		      && w_N.w3 == p3-((uint32)1<<n3)))
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
				var fftp3m_word tmp;
				tmp = x[i2];
				subp3m(x[i1],tmp, x[i2]);
				addp3m(x[i1],tmp, x[i1]);
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
					var fftp3m_word tmp;
					mulp3m(x[i2],w[exp], tmp);
					subp3m(x[i1],tmp, x[i2]);
					addp3m(x[i1],tmp, x[i1]);
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
				var fftp3m_word tmp;
				tmp = y[i2];
				subp3m(y[i1],tmp, y[i2]);
				addp3m(y[i1],tmp, y[i1]);
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
					var fftp3m_word tmp;
					mulp3m(y[i2],w[exp], tmp);
					subp3m(y[i1],tmp, y[i2]);
					addp3m(y[i1],tmp, y[i1]);
				}
			}
		}
	}
	// Multiply the transformed vectors into z.
	for (i = 0; i < N; i++)
		mulp3m(x[i],y[i], z[i]);
	// Undo an FFT of length N on z.
	{
		var uintL l;
		for (l = 0; l < n-1; l++) {
			var const uintC smax = (uintC)1 << (n-1-l);
			var const uintC tmax = (uintC)1 << l;
			#if FFTP3M_BACKWARD != CLEVER
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				#if FFTP3M_BACKWARD == RECIPROOT
				if (exp > 0)
					exp = N - exp; // negate exp (use w^-1 instead of w)
				#endif
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (z(i1),z(i2)) by
					// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/(2*w^exp)).
					var fftp3m_word sum;
					var fftp3m_word diff;
					addp3m(z[i1],z[i2], sum);
					subp3m(z[i1],z[i2], diff);
					shiftp3m(sum, z[i1]);
					mulp3m(diff,w[exp], diff); shiftp3m(diff, z[i2]);
				}
			}
			#else // FFTP3M_BACKWARD == CLEVER: clever handling of negative exponents
			/* s = 0, exp = 0 */ {
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (z(i1),z(i2)) by
					// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/(2*w^exp)),
					// with exp <-- 0.
					var fftp3m_word sum;
					var fftp3m_word diff;
					addp3m(z[i1],z[i2], sum);
					subp3m(z[i1],z[i2], diff);
					shiftp3m(sum, z[i1]);
					shiftp3m(diff, z[i2]);
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
					var fftp3m_word sum;
					var fftp3m_word diff;
					addp3m(z[i1],z[i2], sum);
					subp3m(z[i2],z[i1], diff); // note that w^(N/2) = -1
					shiftp3m(sum, z[i1]);
					mulp3m(diff,w[exp], diff); shiftp3m(diff, z[i2]);
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
				var fftp3m_word sum;
				var fftp3m_word diff;
				addp3m(z[i1],z[i2], sum);
				subp3m(z[i1],z[i2], diff);
				shiftp3m(sum, z[i1]);
				shiftp3m(diff, z[i2]);
			}
		}
	}
	#if FFTP3M_BACKWARD == FORWARD
	// Swap z[i] and z[N-i] for 0 < i < N/2.
	for (i = (N>>1)-1; i > 0; i--) {
		var fftp3m_word tmp = z[i];
		z[i] = z[N-i];
		z[N-i] = tmp;
	}
	#endif
}

static void mulu_fft_modp3m (const uintD* sourceptr1, uintC len1,
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
	var fftp3m_word* const x = cl_alloc_array(fftp3m_word,N);
	var fftp3m_word* const y = cl_alloc_array(fftp3m_word,N);
	#ifdef DEBUG_FFTP3M
	var fftp3m_word* const z = cl_alloc_array(fftp3m_word,N);
	#else
	var fftp3m_word* const z = x; // put z in place of x - saves memory
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
					setp3m(lspref(sourceptr1,i), x[i]);
				for (i = len1; i < N; i++)
					zerop3m(x[i]);
			}
			// Fill factor y.
			if (!squaring) {
				for (i = 0; i < len2p; i++)
					setp3m(lspref(sourceptr2,i), y[i]);
				for (i = len2p; i < N; i++)
					zerop3m(y[i]);
			}
			// Multiply.
			if (!squaring)
				fftp3m_convolution(n,N, &x[0], &y[0], &z[0]);
			else
				fftp3m_convolution(n,N, &x[0], &x[0], &z[0]);
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
					combinep3m(z[i],arrayLSDptr(z_i,3));
					#ifdef DEBUG_FFTP3M
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
			#ifdef DEBUG_FFTP3M
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

#undef n3
#undef n2
#undef n1
#undef p3
#undef p2
#undef p1
