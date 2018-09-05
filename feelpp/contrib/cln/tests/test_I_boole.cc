#include "test_I.h"

int test_I_boole (int iterations)
{
	int error = 0;
	int i;
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(boole(boole_clr,a,b) == 0, a,b);
		ASSERT2(boole(boole_set,a,b) == -1, a,b);
		ASSERT2(boole(boole_1,a,b) == a, a,b);
		ASSERT2(boole(boole_2,a,b) == b, a,b);
		ASSERT2(boole(boole_c1,a,b) == lognot(a), a,b);
		ASSERT2(boole(boole_c2,a,b) == lognot(b), a,b);
		ASSERT2(boole(boole_and,a,b) == logand(a,b), a,b);
		ASSERT2(boole(boole_ior,a,b) == logior(a,b), a,b);
		ASSERT2(boole(boole_xor,a,b) == logxor(a,b), a,b);
		ASSERT2(boole(boole_eqv,a,b) == logeqv(a,b), a,b);
		ASSERT2(boole(boole_nand,a,b) == lognand(a,b), a,b);
		ASSERT2(boole(boole_nor,a,b) == lognor(a,b), a,b);
		ASSERT2(boole(boole_andc1,a,b) == logandc1(a,b), a,b);
		ASSERT2(boole(boole_andc2,a,b) == logandc2(a,b), a,b);
		ASSERT2(boole(boole_orc1,a,b) == logorc1(a,b), a,b);
		ASSERT2(boole(boole_orc2,a,b) == logorc2(a,b), a,b);
	}
	return error;
}
