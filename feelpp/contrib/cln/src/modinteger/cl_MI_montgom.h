// m > 1 odd, Montgomery representation

namespace cln {

// We use Montgomery's modular multiplication trick
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

// Here, we choose to use Montgomery representation only if |V| can be chosen
// to be very small, and in that case we compute d*V mod M using standard
// multiplication and division.
// So we choose N = 2^n with 0 < n <= m (the larger n, the better) and hope
// that it will yield V with |V| < 2^k. We thus replace a division of
// 2m bits by m bits (cost: approx. m^2) by a multiplication of n bits with
// k bits (cost: approx. n*k) and a division of max(2m-n,n+k) bits by m bits
// (cost: approx. (max(2m-n,n+k)-m)*m). Of course, U*M+V*N=1 implies (roughly)
// n+k >= m. It is worth it when
//         m^2 > n*k + (n+k-m)*m   and   m^2 > n*k + (m-n)*m
// <==>  3*m^2 > (n+m)*(k+m)       and     m > k
// <==   3/2*m > k+m                (assume n to be near m)
// <==>    m/2 > k .
//
// How to find N and V:
// U*M+V*N=1 means that U = (M mod 2^n)^-1 = U_m mod 2^n, where
// U_m := (M mod 2^m)^-1 (2-adic reciprocal). |V| < 2^(m/2) is more or less
// equivalent to |V*N| < 2^(n+m/2) <==> |U|*M < 2^(n+m/2) <==> |U| < n-m/2
// <==> the most significant m/2 bits of |U| are all equal. So we search
// for a bit string of at least m/2+1 equal bits in U_m, which has m bits.
// Very easy: take the middle bit of U_m, look how many bits adjacent to it
// (at the left and at the right) have the same value. Say these are the
// bits n-1,...,n-l. (Choose n and l as large as possible. k = m-l + O(1).)
// If l < m/2, forget it. Else fix n and compute V = (1-U*M)/2^n.
//
// It is now clear that only very few moduli M will allow such a good
// choice of N and V, but in these cases the Montgomery multiplication
// reduces the multiplication complexity by a large constant factor.


class cl_heap_modint_ring_montgom : public cl_heap_modint_ring {
	SUBCLASS_cl_heap_modint_ring()
public:
	// Constructor.
	cl_heap_modint_ring_montgom (const cl_I& M, uintL m, uintL n, const cl_I& V);
	// Destructor.
	~cl_heap_modint_ring_montgom () {}
	// Additional information.
	uintL m; // M = 2^m
	uintL n; // N = 2^n, n <= m
	cl_I V;
};

static void cl_modint_ring_montgom_destructor (cl_heap* pointer)
{
	(*(cl_heap_modint_ring_montgom*)pointer).~cl_heap_modint_ring_montgom();
}

cl_class cl_class_modint_ring_montgom = {
	cl_modint_ring_montgom_destructor,
	cl_class_flags_modint_ring
};

// Assuming 0 <= x < 2^(2m), return  V*x mod M.
static inline const cl_I montgom_redc (cl_heap_modint_ring_montgom* R, const cl_I& x)
{
	return mod((x >> R->n) + (R->V * ldb(x,cl_byte(R->n,0))), R->modulus);
}

static const _cl_MI montgom_canonhom (cl_heap_modint_ring* _R, const cl_I& x)
{
	var cl_heap_modint_ring_montgom* R = (cl_heap_modint_ring_montgom*)_R;
	return _cl_MI(R, mod(x << R->n, R->modulus));
}

static const cl_I montgom_retract (cl_heap_modint_ring* _R, const _cl_MI& x)
{
	var cl_heap_modint_ring_montgom* R = (cl_heap_modint_ring_montgom*)_R;
	return montgom_redc(R,x.rep);
}

static const _cl_MI montgom_one (cl_heap_modint_ring* _R)
{
	var cl_heap_modint_ring_montgom* R = (cl_heap_modint_ring_montgom*)_R;
	var cl_I zr = (cl_I)1 << R->n;
	return _cl_MI(R, R->n == R->m ? zr - R->modulus : zr);
}

static const _cl_MI montgom_mul (cl_heap_modint_ring* _R, const _cl_MI& x, const _cl_MI& y)
{
	var cl_heap_modint_ring_montgom* R = (cl_heap_modint_ring_montgom*)_R;
	return _cl_MI(R, montgom_redc(R,x.rep * y.rep));
}

static const _cl_MI montgom_square (cl_heap_modint_ring* _R, const _cl_MI& x)
{
	var cl_heap_modint_ring_montgom* R = (cl_heap_modint_ring_montgom*)_R;
	return _cl_MI(R, montgom_redc(R,square(x.rep)));
}

static const cl_MI_x montgom_recip (cl_heap_modint_ring* _R, const _cl_MI& x)
{
	var cl_heap_modint_ring_montgom* R = (cl_heap_modint_ring_montgom*)_R;
	var const cl_I& xr = x.rep;
	var cl_I u, v;
	var cl_I g = xgcd(xr,R->modulus,&u,&v);
	// g = gcd(x,M) = x*u+M*v
	if (eq(g,1))
		return cl_MI(R, mod((minusp(u) ? u + R->modulus : u) << (2*R->n), R->modulus));
	if (zerop(xr))
		throw division_by_0_exception();
	return cl_notify_composite(R,xr);
}

static const cl_MI_x montgom_div (cl_heap_modint_ring* _R, const _cl_MI& x, const _cl_MI& y)
{
	var cl_heap_modint_ring_montgom* R = (cl_heap_modint_ring_montgom*)_R;
	var const cl_I& yr = y.rep;
	var cl_I u, v;
	var cl_I g = xgcd(yr,R->modulus,&u,&v);
	// g = gcd(y,M) = y*u+M*v
	if (eq(g,1))
		return cl_MI(R, mod((x.rep * (minusp(u) ? u + R->modulus : u)) << R->n, R->modulus));
	if (zerop(yr))
		throw division_by_0_exception();
	return cl_notify_composite(R,yr);
}

#define montgom_addops std_addops
static cl_modint_mulops montgom_mulops = {
	montgom_one,
	montgom_canonhom,
	montgom_mul,
	montgom_square,
	std_expt_pos,
	montgom_recip,
	montgom_div,
	std_expt,
	std_reduce_modulo,
	montgom_retract
};

// Constructor.
inline cl_heap_modint_ring_montgom::cl_heap_modint_ring_montgom (const cl_I& M, uintL _m, uintL _n, const cl_I& _V)
	: cl_heap_modint_ring (M, &std_setops, &montgom_addops, &montgom_mulops),
	  m (_m), n (_n), V (_V)
{
	type = &cl_class_modint_ring_montgom;
}

static cl_heap_modint_ring* try_make_modint_ring_montgom (const cl_I& M)
{
	if (!oddp(M))
		return NULL;
	var uintC m = integer_length(M);
	CL_ALLOCA_STACK;
	var uintC len;
	var const uintD* M_LSDptr;
	I_to_NDS_nocopy(M, ,len=,M_LSDptr=,false,);
	if (lspref(M_LSDptr,len-1)==0) { len--; } // normalize
	// Compute U as 2-adic inverse of M.
	var uintD* U_LSDptr;
	num_stack_alloc(len,,U_LSDptr=);
	recip2adic(len,M_LSDptr,U_LSDptr);
	// Look at U's bits.
	#define U_bit(i) (lspref(U_LSDptr,floor(i,intDsize)) & ((uintD)1 << ((i)%intDsize)))
	var uintC i_min;
	var uintC i_max;
	var uintC i = floor(m,2);
	var bool negative;
	if (U_bit(i)) {
		for (; --i > 0; )
			if (!U_bit(i)) break;
		i_min = i+1;
		i = floor(m,2);
		for (; ++i < m; )
			if (!U_bit(i)) break;
		i_max = i;
		negative = true;
	} else {
		for (; --i > 0; )
			if (U_bit(i)) break;
		i_min = i+1;
		i = floor(m,2);
		for (; ++i < m; )
			if (U_bit(i)) break;
		i_max = i;
		negative = false;
	}
	#undef U_bit
	// OK, all the bits i_max-1..i_min of U are equal.
	if (i_max - i_min <= floor(m,2))
		return NULL;
	var uintC n = i_max;
	// Turn U (mod 2^n) into a signed integer.
	if (n % intDsize) {
		if (negative)
			lspref(U_LSDptr,floor(n,intDsize)) |= (uintD)(-1) << (n % intDsize);
		else
			lspref(U_LSDptr,floor(n,intDsize)) &= ((uintD)1 << (n % intDsize)) - 1;
	}
	var uintC U_len = ceiling(n,intDsize);
	var cl_I U = DS_to_I(U_LSDptr lspop U_len,U_len);
	var cl_I V_N = 1 - U*M;
	if (ldb_test(V_N,cl_byte(n,0)))
		throw runtime_exception();
	var cl_I V = V_N >> n;
	return new cl_heap_modint_ring_montgom(M,m,n,V);
}

}  // namespace cln
