// 1 < m < 2^32, standard representation

namespace cln {

static const _cl_MI int32_plus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 xr = cl_I_to_UL(x.rep);
	var uint32 yr = cl_I_to_UL(y.rep);
	var uint32 zr = xr + yr;
	var uint32 m = cl_I_to_UL(R->modulus);
	if ((zr < xr) || (zr >= m)) { zr = zr - m; }
	return _cl_MI(R, UL_to_I(zr));
}

static const _cl_MI int32_minus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 xr = cl_I_to_UL(x.rep);
	var uint32 yr = cl_I_to_UL(y.rep);
	var sint32 zr = (xr >= yr ? xr - yr : xr - yr + cl_I_to_UL(R->modulus));
	return _cl_MI(R, UL_to_I(zr));
}

static const _cl_MI int32_uminus (cl_heap_modint_ring* R, const _cl_MI& x)
{
	var uint32 xr = cl_I_to_UL(x.rep);
	var uint32 zr = (xr==0 ? 0 : cl_I_to_UL(R->modulus)-xr);
	return _cl_MI(R, UL_to_I(zr));
}

static const _cl_MI int32_mul (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 xr = cl_I_to_UL(x.rep);
	var uint32 yr = cl_I_to_UL(y.rep);
	var uint32 zrhi;
	var uint32 zrlo;
	mulu32(xr,yr,zrhi=,zrlo=);
	var uint32 zr;
	divu_6432_3232(zrhi,zrlo,cl_I_to_UL(R->modulus),,zr=);
	return _cl_MI(R, UL_to_I(zr));
}

static const _cl_MI int32_square (cl_heap_modint_ring* R, const _cl_MI& x)
{
	var uint32 xr = cl_I_to_UL(x.rep);
	var uint32 zrhi;
	var uint32 zrlo;
	mulu32(xr,xr,zrhi=,zrlo=);
	var uint32 zr;
	divu_6432_3232(zrhi,zrlo,cl_I_to_UL(R->modulus),,zr=);
	return _cl_MI(R, UL_to_I(zr));
}

static cl_modint_addops int32_addops = {
	std_zero,
	std_zerop,
	int32_plus,
	int32_minus,
	int32_uminus
};
static cl_modint_mulops int32_mulops = {
	std_one,
	std_canonhom,
	int32_mul,
	int32_square,
	std_expt_pos,
	std_recip,
	std_div,
	std_expt,
	std_reduce_modulo,
	std_retract
};

class cl_heap_modint_ring_int32 : public cl_heap_modint_ring {
	SUBCLASS_cl_heap_modint_ring()
public:
	// Constructor.
	cl_heap_modint_ring_int32 (const cl_I& m);
	// Destructor.
	~cl_heap_modint_ring_int32 () {}
};

static void cl_modint_ring_int32_destructor (cl_heap* pointer)
{
	(*(cl_heap_modint_ring_int32*)pointer).~cl_heap_modint_ring_int32();
}

cl_class cl_class_modint_ring_int32 = {
	cl_modint_ring_int32_destructor,
	cl_class_flags_modint_ring
};

// Constructor.
inline cl_heap_modint_ring_int32::cl_heap_modint_ring_int32(const cl_I& m)
	: cl_heap_modint_ring (m, &std_setops, &int32_addops, &int32_mulops)
{
	type = &cl_class_modint_ring_int32;
}

}  // namespace cln
