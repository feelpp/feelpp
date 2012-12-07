// m > 0, m = 2^m1

namespace cln {

class cl_heap_modint_ring_pow2 : public cl_heap_modint_ring {
	SUBCLASS_cl_heap_modint_ring()
public:
	// Constructor.
	cl_heap_modint_ring_pow2 (const cl_I& m, uintC m1); // m = 2^m1
	// Destructor.
	~cl_heap_modint_ring_pow2 () {}
	// Additional information.
	uintC m1;
};

static inline const cl_I pow2_reduce_modulo (cl_heap_modint_ring* _R, const cl_I& x)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	return ldb(x,cl_byte(R->m1,0));
}

static const _cl_MI pow2_canonhom (cl_heap_modint_ring* R, const cl_I& x)
{
	return _cl_MI(R, pow2_reduce_modulo(R,x));
}

static const _cl_MI pow2_plus (cl_heap_modint_ring* _R, const _cl_MI& x, const _cl_MI& y)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	var cl_I zr = x.rep + y.rep;
	return _cl_MI(R, ldb(zr,cl_byte(R->m1,0)));
}

static const _cl_MI pow2_minus (cl_heap_modint_ring* _R, const _cl_MI& x, const _cl_MI& y)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	var cl_I zr = x.rep - y.rep;
	return _cl_MI(R, ldb(zr,cl_byte(R->m1,0)));
}

static const _cl_MI pow2_uminus (cl_heap_modint_ring* _R, const _cl_MI& x)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	var cl_I zr = - x.rep;
	return _cl_MI(R, ldb(zr,cl_byte(R->m1,0)));
}

static const _cl_MI pow2_one (cl_heap_modint_ring* _R)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	return _cl_MI(R, R->m1==0 ? 0 : 1);
}

static const _cl_MI pow2_mul (cl_heap_modint_ring* _R, const _cl_MI& x, const _cl_MI& y)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	var cl_I zr = x.rep * y.rep;
	return _cl_MI(R, ldb(zr,cl_byte(R->m1,0)));
}

static const _cl_MI pow2_square (cl_heap_modint_ring* _R, const _cl_MI& x)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	var cl_I zr = square(x.rep);
	return _cl_MI(R, ldb(zr,cl_byte(R->m1,0)));
}

// Timing comparison with std_recip, on a i486 33 MHz running Linux:
// timeMIpow2recip N inverts an (N*32)-bit number.
//   N  std_recip pow2_recip
//   10   0.0030  0.00017
//   25   0.011   0.00068
//   50   0.035   0.0024
//  100   0.124   0.0090
//  250   0.71    0.055
//  500   2.76    0.193
// 1000  11.0     0.61
// 2500  68.7     2.2
// 5000 283       5.0
static const cl_MI_x pow2_recip (cl_heap_modint_ring* _R, const _cl_MI& x)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	var const cl_I& xr = x.rep;
	if (!oddp(xr)) {
		if (R->m1 == 0)
			return cl_MI(R, 0);
		if (zerop(xr))
			throw division_by_0_exception();
		else
			return cl_notify_composite(R,xr);
	} else
		return cl_MI(R, cl_recip2adic(R->m1,xr));
}

// Timing comparison with std_div, on a i486 33 MHz running Linux:
// timeMIpow2div N divides two (N*32)-bit numbers.
//   N   std_div  pow2_div
//   10   0.0035  0.00017
//   25   0.0136  0.00068
//   50   0.043   0.0024
//  100   0.151   0.0090
//  250   0.85    0.054
//  500   3.3     0.21
// 1000  12.3     0.86
static const cl_MI_x pow2_div (cl_heap_modint_ring* _R, const _cl_MI& x, const _cl_MI& y)
{
	var cl_heap_modint_ring_pow2* R = (cl_heap_modint_ring_pow2*)_R;
	var const cl_I& yr = y.rep;
	if (!oddp(yr)) {
		if (R->m1 == 0)
			return cl_MI(R, 0);
		if (zerop(yr))
			throw division_by_0_exception();
		else
			return cl_notify_composite(R,yr);
	} else
		return cl_MI(R, cl_div2adic(R->m1,x.rep,yr));
}

static cl_modint_addops pow2_addops = {
	std_zero,
	std_zerop,
	pow2_plus,
	pow2_minus,
	pow2_uminus
};
static cl_modint_mulops pow2_mulops = {
	pow2_one,
	pow2_canonhom,
	pow2_mul,
	pow2_square,
	std_expt_pos,
	pow2_recip,
	pow2_div,
	std_expt,
	pow2_reduce_modulo,
	std_retract
};

static void cl_modint_ring_pow2_destructor (cl_heap* pointer)
{
	(*(cl_heap_modint_ring_pow2*)pointer).~cl_heap_modint_ring_pow2();
}

cl_class cl_class_modint_ring_pow2 = {
	cl_modint_ring_pow2_destructor,
	cl_class_flags_modint_ring
};

// Constructor.
inline cl_heap_modint_ring_pow2::cl_heap_modint_ring_pow2 (const cl_I& m, uintC _m1)
	: cl_heap_modint_ring (m, &std_setops, &pow2_addops, &pow2_mulops), m1 (_m1)
{
	type = &cl_class_modint_ring_pow2;
}

}  // namespace cln
