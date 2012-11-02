// 1 < m < 2^(cl_value_len-1), standard representation
// Assuming (cl_value_len <= 32).

namespace cln {

static const _cl_MI fix29_plus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 zr = FN_to_UV(x.rep) + FN_to_UV(y.rep);
	if (zr >= FN_to_UV(R->modulus)) { zr = zr - FN_to_UV(R->modulus); }
	return _cl_MI(R, L_to_FN(zr));
}

static const _cl_MI fix29_minus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 xr = FN_to_UV(x.rep);
	var uint32 yr = FN_to_UV(y.rep);
	var sint32 zr = xr - yr;
	if (zr < 0) { zr = zr + FN_to_UV(R->modulus); }
	return _cl_MI(R, L_to_FN(zr));
}

static const _cl_MI fix29_uminus (cl_heap_modint_ring* R, const _cl_MI& x)
{
	var uint32 xr = FN_to_UV(x.rep);
	var uint32 zr = (xr==0 ? 0 : FN_to_UV(R->modulus)-xr);
	return _cl_MI(R, L_to_FN(zr));
}

static const _cl_MI fix29_mul (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 xr = FN_to_UV(x.rep);
	var uint32 yr = FN_to_UV(y.rep);
	var uint32 zrhi;
	var uint32 zrlo;
	mulu32(xr,yr,zrhi=,zrlo=);
	var uint32 zr;
	divu_6432_3232(zrhi,zrlo,FN_to_UV(R->modulus),,zr=);
	return _cl_MI(R, L_to_FN(zr));
}

static const _cl_MI fix29_square (cl_heap_modint_ring* R, const _cl_MI& x)
{
	var uint32 xr = FN_to_UV(x.rep);
	var uint32 zrhi;
	var uint32 zrlo;
	mulu32(xr,xr,zrhi=,zrlo=);
	var uint32 zr;
	divu_6432_3232(zrhi,zrlo,FN_to_UV(R->modulus),,zr=);
	return _cl_MI(R, L_to_FN(zr));
}

static cl_modint_addops fix29_addops = {
	std_zero,
	std_zerop,
	fix29_plus,
	fix29_minus,
	fix29_uminus
};
static cl_modint_mulops fix29_mulops = {
	std_one,
	std_canonhom,
	fix29_mul,
	fix29_square,
	std_expt_pos,
	std_recip,
	std_div,
	std_expt,
	std_reduce_modulo,
	std_retract
};

class cl_heap_modint_ring_fix29 : public cl_heap_modint_ring {
	SUBCLASS_cl_heap_modint_ring()
public:
	// Constructor.
	cl_heap_modint_ring_fix29 (const cl_I& m);
	// Destructor.
	~cl_heap_modint_ring_fix29 () {}
};

static void cl_modint_ring_fix29_destructor (cl_heap* pointer)
{
	(*(cl_heap_modint_ring_fix29*)pointer).~cl_heap_modint_ring_fix29();
}

cl_class cl_class_modint_ring_fix29 = {
	cl_modint_ring_fix29_destructor,
	cl_class_flags_modint_ring
};

// Constructor.
inline cl_heap_modint_ring_fix29::cl_heap_modint_ring_fix29(const cl_I& m)
	: cl_heap_modint_ring (m, &std_setops, &fix29_addops, &fix29_mulops)
{
	type = &cl_class_modint_ring_fix29;
}

}  // namespace cln
