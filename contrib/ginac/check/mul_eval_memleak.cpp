/** @file mul_eval_memleak.cpp
 *
 *  mul_eval_memleak.cpp Test for memory leak in {mul,power}::eval
 *
 *  The bug was introduced by
 *
 *  commit f418c6ee4d558c852e1fb95533af07a3ae43f409
 *  Author: Alexei Sheplyakov <varg@theor.jinr.ru>
 *  Date:   Wed Jul 11 14:34:42 2007 +0400
 *  (it was commited into the official branch as
 *  commit a602d34c225dceb3e53742a7b3e19a4b5e280485
 *  Author: Jens Vollinga <vollinga@thep.physik.uni-mainz.de>
 *  Date:   Wed Jul 11 21:07:40 2007 +0000)
 */

/*
 *  GiNaC Copyright (C) 1999-2011 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <ginac/ginac.h>
using namespace GiNaC;

#include <iostream>
#include <stdexcept>
#include <string>
using namespace std;

const unsigned check_mul_eval_memleak(const unsigned N)
{
	// Expanding this expression into a Laurent series triggers the bug.
	static const string e_str("\
1/605927415293858601*tgamma(3-eps)^(-1)*tgamma(2-eps)*(5013234896802\
*(-1+2*eps)*((-2539001/2)*eps^(-2)+(-7617003/2+1692800*log(920*scale\
^(-1))+3385202*log(1301*scale^(-1)))*eps^(-1)+F)+(846400+2539001*eps\
)*(2539001*(-1+2*eps)*((-2539001/2)*eps^(-2)+(-7617003/2+1692800*log\
(920*scale^(-1))+3385202*log(1301*scale^(-1)))*eps^(-1)+F)+286523497\
2800*(-1+eps)^(-1)*eps^(-2)*(920*scale^(-1))^(-2*eps)*(1301*scale^(-\
1))^(-2*eps)*tgamma(1+eps)^2)+6061411748045832000*(-1+eps)^(-1)*eps^\
(-2)*(920*scale^(-1))^(-2*eps)*(1301*scale^(-1))^(-2*eps)*tgamma(1+e\
ps)^2+716056132401*((-2539001/2)*eps^(-2)+(-7617003/2+1692800*log(92\
0*scale^(-1))+3385202*log(1301*scale^(-1)))*eps^(-1)+F))+71656139360\
0/716056132401*tgamma(1-eps)*tgamma(2*eps)^(-1)*tgamma(eps)^2*tgamma\
(3-eps)^(-1)*(920*scale^(-1))^(-4*eps)*tgamma(-1+2*eps)-2/6059274152\
93858601*tgamma(3-eps)^(-1)*(2149010446400*(-1+2*eps)*((-2539001/2)*\
eps^(-2)+(-7617003/2+1692800*log(920*scale^(-1))+3385202*log(1301*sc\
ale^(-1)))*eps^(-1)+F)+2425134880977920000*(-1+eps)^(-1)*eps^(-2)*(9\
20*scale^(-1))^(-2*eps)*(1301*scale^(-1))^(-2*eps)*tgamma(1+eps)^2-7\
16056132401*((-2539001/2)*eps^(-2)+(-7617003/2+1692800*log(920*scale\
^(-1))+3385202*log(1301*scale^(-1)))*eps^(-1)+F)+1692601*(-1+2*eps)*\
(2539001*(-1+2*eps)*((-2539001/2)*eps^(-2)+(-7617003/2+1692800*log(9\
20*scale^(-1))+3385202*log(1301*scale^(-1)))*eps^(-1)+F)+28652349728\
00*(-1+eps)^(-1)*eps^(-2)*(920*scale^(-1))^(-2*eps)*(1301*scale^(-1)\
)^(-2*eps)*tgamma(1+eps)^2))*tgamma(2-eps)+(1/716056132401*I)*tgamma\
(-1+eps)*mb^2*(mb*scale^(-1))^(-2*eps)*((716392960000*I)*mb^(-2)*(92\
0*scale^(-1))^(-2*eps)*tgamma(-2+eps)-(2864898145201*I)*mb^(-2)*(130\
1*scale^(-1))^(-2*eps)*tgamma(-2+eps)-(716224526400*I)*tgamma(-1+eps\
)*mb^(-2)*(920*scale^(-1))^(-2*eps))-3385202/605927415293858601*tgam\
ma(3-eps)^(-1)*tgamma(2-eps)*(2539001*(-1+2*eps)*((-2539001/2)*eps^(\
-2)+(-7617003/2+1692800*log(920*scale^(-1))+3385202*log(1301*scale^(\
-1)))*eps^(-1)+F)+2865234972800*(-1+eps)^(-1)*eps^(-2)*(920*scale^(-\
1))^(-2*eps)*(1301*scale^(-1))^(-2*eps)*tgamma(1+eps)^2+846201*((-25\
39001/2)*eps^(-2)+(-7617003/2+1692800*log(920*scale^(-1))+3385202*lo\
g(1301*scale^(-1)))*eps^(-1)+F))\
");
	const symbol eps("eps"), scale("scale"), mb("mb"), F("F");
	const lst syms(eps, scale, mb, F);
	const ex e0(e_str, syms);

	unsigned i = 0;
	unsigned n_failures = 0;

	ex e;
	try {
		for (; i < N; i++)
			e = e0.series(eps, 1).subs(Euler==0).expand();
	} catch (std::bad_alloc) {
		return i;
	}
	return 0;
}
			
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

static void set_VM_limit(const unsigned long MB) {
	const unsigned mem_lim_kb = MB*1024*1024;
	struct rlimit lim;
	lim.rlim_cur = mem_lim_kb;
	lim.rlim_max = mem_lim_kb;
	setrlimit(RLIMIT_AS, &lim);
}

int main(int argc, char** argv) {
	static const unsigned max_mem = 32; // megabytes
	// otherwise one need wait for a long[er] time.
	set_VM_limit(max_mem); 
	static const unsigned n_of_tests = 10000;
	const unsigned n_loops = check_mul_eval_memleak(n_of_tests);
	if (n_loops) {
		cerr << "memory exhausted after " << n_loops << " loops" << endl;
		return 1;
	}
	return 0;
}
