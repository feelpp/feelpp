// m = 0 : Z/mZ \isomorph Z

namespace cln {

static void int_fprint (cl_heap_modint_ring* R, std::ostream& stream, const _cl_MI &x)
{
	fprint(stream,R->_retract(x));
}

static const cl_I int_reduce_modulo (cl_heap_modint_ring* R, const cl_I& x)
{
	unused R;
	return x; // reducing modulo 0 does nothing
}

// This is the only case where canonhom is injective.
static const _cl_MI int_canonhom (cl_heap_modint_ring* R, const cl_I& x)
{
	return _cl_MI(R, x);
}

// This is the only case where retract is surjective.
static const cl_I int_retract (cl_heap_modint_ring* R, const _cl_MI& x)
{
	unused R;
	return x.rep;
}

// This is the only case where random yields an error.
static const _cl_MI int_random (cl_heap_modint_ring* R, random_state& randomstate)
{
	unused R;
	unused randomstate;
	throw runtime_exception("Z / 0 Z not a finite set - no equidistributed random function.");
#if ((defined(__sparc__) || defined(__sparc64__)) && !defined(__GNUC__)) // Sun CC wants a return value
	return _cl_MI(R, 0);
#endif
}

static const _cl_MI int_zero (cl_heap_modint_ring* R)
{
	return _cl_MI(R, 0);
}

static bool int_zerop (cl_heap_modint_ring* R, const _cl_MI& x)
{
	unused R;
	return zerop(x.rep);
}

static const _cl_MI int_plus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	return _cl_MI(R, x.rep + y.rep);
}

static const _cl_MI int_minus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	return _cl_MI(R, x.rep - y.rep);
}

static const _cl_MI int_uminus (cl_heap_modint_ring* R, const _cl_MI& x)
{
	return _cl_MI(R, - x.rep);
}

static const _cl_MI int_one (cl_heap_modint_ring* R)
{
	return _cl_MI(R, 1);
}

static const _cl_MI int_mul (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	return _cl_MI(R, x.rep * y.rep);
}

static const _cl_MI int_square (cl_heap_modint_ring* R, const _cl_MI& x)
{
	return _cl_MI(R, square(x.rep));
}

static const cl_MI_x int_recip (cl_heap_modint_ring* R, const _cl_MI& x)
{
	var const cl_I& xr = x.rep;
	if (eq(xr,1) || eq(xr,-1)) { return cl_MI(R,x); }
	if (zerop(xr)) { throw division_by_0_exception(); }
	return cl_notify_composite(R,xr);
}

static const cl_MI_x int_div (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var const cl_I& yr = y.rep;
	if (eq(yr,1)) { return cl_MI(R,x.rep); }
	if (eq(yr,-1)) { return cl_MI(R,-x.rep); }
	if (zerop(yr)) { throw division_by_0_exception(); }
	return cl_notify_composite(R,yr);
}

static const _cl_MI int_expt_pos (cl_heap_modint_ring* R, const _cl_MI& x, const cl_I& y)
{
	return _cl_MI(R, expt_pos(x.rep,y));
}

static const cl_MI_x int_expt (cl_heap_modint_ring* R, const _cl_MI& x, const cl_I& y)
{
	if (eq(x.rep,1)) { return cl_MI(R,1); }
	if (eq(x.rep,-1)) { return cl_MI(R,evenp(y)?1:-1); }
	if (!minusp(y)) {
		if (zerop(y))
			return cl_MI(R,1);
		else
			return cl_MI(R,expt_pos(x.rep,y));
	}
	// y < 0, x nonunit.
	if (zerop(x.rep)) { throw division_by_0_exception(); }
	return cl_notify_composite(R,x.rep);
}

static cl_modint_setops int_setops = {
	int_fprint,
	modint_equal,
	int_random
};
static cl_modint_addops int_addops = {
	int_zero,
	int_zerop,
	int_plus,
	int_minus,
	int_uminus
};
static cl_modint_mulops int_mulops = {
	int_one,
	int_canonhom,
	int_mul,
	int_square,
	int_expt_pos,
	int_recip,
	int_div,
	int_expt,
	int_reduce_modulo,
	int_retract
};

class cl_heap_modint_ring_int : public cl_heap_modint_ring {
	SUBCLASS_cl_heap_modint_ring()
public:
	// Constructor.
	cl_heap_modint_ring_int () : cl_heap_modint_ring (0, &int_setops, &int_addops, &int_mulops) {}
	// Virtual destructor.
	~cl_heap_modint_ring_int () {}
};

}  // namespace cln
