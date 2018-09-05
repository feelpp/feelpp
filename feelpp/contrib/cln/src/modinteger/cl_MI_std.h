// m > 1, standard representation, no tricks

namespace cln {

static void std_fprint (cl_heap_modint_ring* R, std::ostream& stream, const _cl_MI &x)
{
	fprint(stream,R->_retract(x));
	fprint(stream," mod ");
	fprint(stream,R->modulus);
}

static const cl_I std_reduce_modulo (cl_heap_modint_ring* R, const cl_I& x)
{
	return mod(x,R->modulus);
}

static const _cl_MI std_canonhom (cl_heap_modint_ring* R, const cl_I& x)
{
	return _cl_MI(R, mod(x,R->modulus));
}

static const cl_I std_retract (cl_heap_modint_ring* R, const _cl_MI& x)
{
	unused R;
	return x.rep;
}

static const _cl_MI std_random (cl_heap_modint_ring* R, random_state& randomstate)
{
	return _cl_MI(R, random_I(randomstate,R->modulus));
}

static const _cl_MI std_zero (cl_heap_modint_ring* R)
{
	return _cl_MI(R, 0);
}

static bool std_zerop (cl_heap_modint_ring* R, const _cl_MI& x)
{
	unused R;
	return zerop(x.rep);
}

static const _cl_MI std_plus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var cl_I zr = x.rep + y.rep;
	return _cl_MI(R, (zr >= R->modulus ? zr - R->modulus : zr));
}

static const _cl_MI std_minus (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var cl_I zr = x.rep - y.rep;
	return _cl_MI(R, (minusp(zr) ? zr + R->modulus : zr));
}

static const _cl_MI std_uminus (cl_heap_modint_ring* R, const _cl_MI& x)
{
	return _cl_MI(R, (zerop(x.rep) ? (cl_I)0 : R->modulus - x.rep));
}

static const _cl_MI std_one (cl_heap_modint_ring* R)
{
	return _cl_MI(R, 1);
}

static const _cl_MI std_mul (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	return _cl_MI(R, mod(x.rep * y.rep, R->modulus));
}

static const _cl_MI std_square (cl_heap_modint_ring* R, const _cl_MI& x)
{
	return _cl_MI(R, mod(square(x.rep), R->modulus));
}

static const cl_MI_x std_recip (cl_heap_modint_ring* R, const _cl_MI& x)
{
	var const cl_I& xr = x.rep;
	var cl_I u, v;
	var cl_I g = xgcd(xr,R->modulus,&u,&v);
	// g = gcd(x,m) = x*u+m*v
	if (eq(g,1))
		return cl_MI(R, (minusp(u) ? u + R->modulus : u));
	if (zerop(xr))
		throw division_by_0_exception();
	return cl_notify_composite(R,xr);
}

static const cl_MI_x std_div (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y)
{
	var const cl_I& yr = y.rep;
	var cl_I u, v;
	var cl_I g = xgcd(yr,R->modulus,&u,&v);
	// g = gcd(y,m) = y*u+m*v
	if (eq(g,1))
		return cl_MI(R, mod(x.rep * (minusp(u) ? u + R->modulus : u), R->modulus));
	if (zerop(yr))
		throw division_by_0_exception();
	return cl_notify_composite(R,yr);
}

static uint8 const ord2_table[256] = // maps i -> ord2(i) for i>0
{
 63,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  5,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  6,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  5,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  7,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  5,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  6,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  5,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0
};

static uint8 const odd_table[256] = // maps i -> i/2^ord2(i) for i>0
{
  0,   1,   1,   3,   1,   5,   3,   7,   1,   9,   5,  11,   3,  13,   7,  15,
  1,  17,   9,  19,   5,  21,  11,  23,   3,  25,  13,  27,   7,  29,  15,  31,
  1,  33,  17,  35,   9,  37,  19,  39,   5,  41,  21,  43,  11,  45,  23,  47,
  3,  49,  25,  51,  13,  53,  27,  55,   7,  57,  29,  59,  15,  61,  31,  63,
  1,  65,  33,  67,  17,  69,  35,  71,   9,  73,  37,  75,  19,  77,  39,  79,
  5,  81,  41,  83,  21,  85,  43,  87,  11,  89,  45,  91,  23,  93,  47,  95,
  3,  97,  49,  99,  25, 101,  51, 103,  13, 105,  53, 107,  27, 109,  55, 111,
  7, 113,  57, 115,  29, 117,  59, 119,  15, 121,  61, 123,  31, 125,  63, 127,
  1, 129,  65, 131,  33, 133,  67, 135,  17, 137,  69, 139,  35, 141,  71, 143,
  9, 145,  73, 147,  37, 149,  75, 151,  19, 153,  77, 155,  39, 157,  79, 159,
  5, 161,  81, 163,  41, 165,  83, 167,  21, 169,  85, 171,  43, 173,  87, 175,
 11, 177,  89, 179,  45, 181,  91, 183,  23, 185,  93, 187,  47, 189,  95, 191,
  3, 193,  97, 195,  49, 197,  99, 199,  25, 201, 101, 203,  51, 205, 103, 207,
 13, 209, 105, 211,  53, 213, 107, 215,  27, 217, 109, 219,  55, 221, 111, 223,
  7, 225, 113, 227,  57, 229, 115, 231,  29, 233, 117, 235,  59, 237, 119, 239,
 15, 241, 121, 243,  61, 245, 123, 247,  31, 249, 125, 251,  63, 253, 127, 255
};

static const _cl_MI std_expt_pos (cl_heap_modint_ring* R, const _cl_MI& x, const cl_I& y)
{
#if 0
	// Methode:
	// Right-Left Binary, [Cohen, Algorithm 1.2.1.]
	//   a:=x, b:=y.
	//   Solange b gerade, setze a:=a*a, b:=b/2. (Invariante: a^b = x^y.)
	//   c:=a.
	//   Solange b:=floor(b/2) >0 ist,
	//     setze a:=a*a, und falls b ungerade, setze c:=a*c.
	//   Liefere c.
	var _cl_MI a = x;
	var cl_I b = y;
	while (!oddp(b)) { a = R->_square(a); b = b >> 1; } // a^b = x^y
	var _cl_MI c = a;
	until (eq(b,1)) {
		b = b >> 1;
		a = R->_square(a);
		// a^b*c = x^y
		if (oddp(b))
			c = R->_mul(a,c);
	}
	return c;
#else
	// Methode:
	// Left-Right base 2^k, [Cohen, Algorithm 1.2.4.]
	// Good values of k, depending on the size nn of the exponent n:
	//   k = 1 for nn <= 8
	//   k = 2 for nn <= 24
	//   k = 3 for nn <= 69.8
	//   ... k for nn <= k*(k+1)*2^(2k)/(2^(k+1)-k-2)
	var cl_I n = y;
	var uintC nn = integer_length(n);
	// n has nn bits.
	if (nn <= 8) {
		// k = 1, normal Left-Right Binary algorithm.
		var uintL _n = FN_to_UV(n);
		var _cl_MI a = x;
		for (var int i = nn-2; i >= 0; i--) {
			a = R->_square(a);
			if (_n & bit(i))
				a = R->_mul(a,x);
		}
		return cl_MI(R,a);
	} else {
		// General Left-Right base 2^k algorithm.
		CL_ALLOCA_STACK;
		var uintL k;
		     if (nn <= 24)	k = 2;
		else if (nn <= 69)	k = 3;
		else if (nn <= 196)	k = 4;
		else if (nn <= 538)	k = 5;
		else if (nn <= 1433)	k = 6;
		else if (nn <= 3714)	k = 7;
		else if (nn <= 9399)	k = 8;
		else if (nn <= 23290)	k = 9;
		else if (nn <= 56651)	k = 10;
		else if (nn <= 135598)	k = 11;
		else if (nn <= 320034)	k = 12;
		else if (nn <= 746155)	k = 13;
		else if (nn <= 1721160)	k = 14;
		else if (nn <= 3933180)	k = 15;
		else /* if (nn <= 8914120) */ k = 16;
		var uintC nnk = ceiling(nn,k); // number of base-2^k digits in n
		var uint16* n_digits = cl_alloc_array(uint16,nnk);
		// Split n into base-2^k digits.
		{
			var const uintD* n_LSDptr;
			var const uintD* n_MSDptr;
			I_to_NDS_nocopy(n, n_MSDptr=,,n_LSDptr=,false,);
			var const uintL k_mask = bit(k)-1;
			var uintD carry = 0;
			var unsigned int carrybits = 0;
			for (var uintC i = 0; i < nnk; i++) {
				if (carrybits >= k) {
					n_digits[i] = carry & k_mask;
					carry = carry >> k;
					carrybits -= k;
				} else {
					var uintD next =
					  (n_LSDptr==n_MSDptr ? 0 : lsprefnext(n_LSDptr));
					n_digits[i] = (carry | (next << carrybits)) & k_mask;
					carry = next >> (k-carrybits);
					carrybits = intDsize - (k-carrybits);
				}
			}
		}
		// Compute maximum odd base-2^k digit.
		var uintL maxodd = 1;
		if (k <= 8) {
			for (var uintC i = 0; i < nnk; i++) {
				var uintL d = n_digits[i];
				if (d > 0) {
					d = odd_table[d];
					if (d > maxodd) maxodd = d;
				}
			}
		} else {
			for (var uintC i = 0; i < nnk; i++) {
				var uintL d = n_digits[i];
				if (d > 0) {
					var uintL d2; ord2_32(d,d2=);
					d = d>>d2;
					if (d > maxodd) maxodd = d;
				}
			}
		}
		maxodd = (maxodd+1)/2; // number of odd powers we need
		var _cl_MI* x_oddpow = cl_alloc_array(_cl_MI,maxodd);
		var _cl_MI x2 = (maxodd > 1 ? R->_square(x) : R->_zero());
		{
			init1(_cl_MI, x_oddpow[0]) (x);
			for (var uintL i = 1; i < maxodd; i++)
				init1(_cl_MI, x_oddpow[i]) (R->_mul(x_oddpow[i-1],x2));
		}
		var _cl_MI a;
		// Compute a = x^n_digits[nnk-1].
		{
			var uintL d = n_digits[nnk-1];
			if (d == 0) throw runtime_exception();
			var uintL d2;
			if (k <= 8)
				d2 = ord2_table[d];
			else
				ord2_32(d,d2=);
			d = d>>d2; // d := d/2^ord2(d)
			d = floor(d,2);
			a = x_oddpow[d]; // x^d
			// Square d2 times.
			if (d==0 && maxodd > 1 && d2>0) {
				a = x2; d2--;
			}
			if (!(d2 < k)) throw runtime_exception();
			for ( ; d2>0; d2--)
				a = R->_square(a);
		}
		for (var sintC i = nnk-2; i >= 0; i--) {
			// Compute a := a^(2^k) * x^n_digits[i].
			var uintL d = n_digits[i];
			var uintL d2;
			if (d > 0) {
				if (k <= 8)
					d2 = ord2_table[d];
				else
					ord2_32(d,d2=);
				d = d>>d2; // d/2^ord2(d)
				d = floor(d,2);
				for (var sintL j = k-d2; j>0; j--)
					a = R->_square(a);
				a = R->_mul(a,x_oddpow[d]);
			} else
				d2 = k;
			// Square d2 times.
			if (!(d2 <= k)) throw runtime_exception();
			for ( ; d2>0; d2--)
				a = R->_square(a);
		}
		{
			for (var uintL i = 0; i < maxodd; i++)
				x_oddpow[i].~_cl_MI();
		}
		return cl_MI(R,a);
	}
#endif
}

static const cl_MI_x std_expt (cl_heap_modint_ring* R, const _cl_MI& x, const cl_I& y)
{
	if (!minusp(y)) {
		if (zerop(y))
			return R->one();
		else
			return cl_MI(R,R->_expt_pos(x,y));
	} else
		return R->_recip(R->_expt_pos(x,-y));
}

static cl_modint_setops std_setops = {
	std_fprint,
	modint_equal,
	std_random
};
static cl_modint_addops std_addops = {
	std_zero,
	std_zerop,
	std_plus,
	std_minus,
	std_uminus
};
static cl_modint_mulops std_mulops = {
	std_one,
	std_canonhom,
	std_mul,
	std_square,
	std_expt_pos,
	std_recip,
	std_div,
	std_expt,
	std_reduce_modulo,
	std_retract
};

class cl_heap_modint_ring_std : public cl_heap_modint_ring {
	SUBCLASS_cl_heap_modint_ring()
public:
	// Constructor.
	cl_heap_modint_ring_std (const cl_I& m);
	// Virtual destructor.
	~cl_heap_modint_ring_std () {}
};

static void cl_heap_modint_ring_std_destructor (cl_heap* pointer)
{
	(*(cl_heap_modint_ring_std*)pointer).~cl_heap_modint_ring_std();
}

cl_class cl_class_modint_ring_std = {
	cl_heap_modint_ring_std_destructor,
	cl_class_flags_modint_ring
};

// Constructor.
inline cl_heap_modint_ring_std::cl_heap_modint_ring_std (const cl_I& m)
	: cl_heap_modint_ring (m, &std_setops, &std_addops, &std_mulops)
{
	type = &cl_class_modint_ring_std;
}

}  // namespace cln
