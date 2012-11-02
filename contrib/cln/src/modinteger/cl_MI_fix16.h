// 1 < m < 2^16, standard representation

namespace cln {

static const _cl_MI fix16_plus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 zr = FN_to_UV(x.rep) + FN_to_UV(y.rep);
	if (zr >= FN_to_UV(R->modulus)) { zr = zr - FN_to_UV(R->modulus); }
	return _cl_MI(R, L_to_FN(zr));
}

static const _cl_MI fix16_minus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 xr = FN_to_UV(x.rep);
	var uint32 yr = FN_to_UV(y.rep);
	var sint32 zr = xr - yr;
	if (zr < 0) { zr = zr + FN_to_UV(R->modulus); }
	return _cl_MI(R, L_to_FN(zr));
}

static const _cl_MI fix16_uminus (cl_heap_modint_ring* R, const _cl_MI& x)
{
	var uint32 xr = FN_to_UV(x.rep);
	var uint32 zr = (xr==0 ? 0 : FN_to_UV(R->modulus)-xr);
	return _cl_MI(R, L_to_FN(zr));
}

static const _cl_MI fix16_mul (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var uint32 xr = FN_to_UV(x.rep);
	var uint32 yr = FN_to_UV(y.rep);
	var uint32 zr = mulu16(xr,yr);
	divu_3216_1616(zr,FN_to_UV(R->modulus),,zr=);
	return _cl_MI(R, L_to_FN(zr));
}

static const _cl_MI fix16_square (cl_heap_modint_ring* R, const _cl_MI& x)
{
	var uint32 xr = FN_to_UV(x.rep);
	var uint32 zr = mulu16(xr,xr);
	divu_3216_1616(zr,FN_to_UV(R->modulus),,zr=);
	return _cl_MI(R, L_to_FN(zr));
}

static cl_modint_addops fix16_addops = {
	std_zero,
	std_zerop,
	fix16_plus,
	fix16_minus,
	fix16_uminus
};
static cl_modint_mulops fix16_mulops = {
	std_one,
	std_canonhom,
	fix16_mul,
	fix16_square,
	std_expt_pos,
	std_recip,
	std_div,
	std_expt,
	std_reduce_modulo,
	std_retract
};

class cl_heap_modint_ring_fix16 : public cl_heap_modint_ring {
	SUBCLASS_cl_heap_modint_ring()
public:
	// Constructor.
	cl_heap_modint_ring_fix16 (const cl_I& m);
	// Destructor.
	~cl_heap_modint_ring_fix16 () {}
};

static void cl_modint_ring_fix16_destructor (cl_heap* pointer)
{
	(*(cl_heap_modint_ring_fix16*)pointer).~cl_heap_modint_ring_fix16();
}

cl_class cl_class_modint_ring_fix16 = {
	cl_modint_ring_fix16_destructor,
	cl_class_flags_modint_ring
};

// Constructor.
inline cl_heap_modint_ring_fix16::cl_heap_modint_ring_fix16(const cl_I& m)
	: cl_heap_modint_ring (m, &std_setops, &fix16_addops, &fix16_mulops)
{
	type = &cl_class_modint_ring_fix16;
}

}  // namespace cln
