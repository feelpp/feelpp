/** @file exam_powerlaws.cpp
 *
 *  Tests for power laws.  You shouldn't try to draw much inspiration from
 *  this code, it is a sanity check rather deeply rooted in GiNaC's classes. */

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

#include "ginac.h"
using namespace GiNaC;

#include <iostream>
using namespace std;

static unsigned exam_powerlaws1()
{
	// (x^a)^b = x^(a*b)
	
	symbol x("x");
	symbol a("a");
	symbol b("b");
	
	ex e1 = power(power(x,a), b);
	if (!(is_exactly_a<power>(e1) &&
	      is_exactly_a<power>(e1.op(0)) &&
	      is_exactly_a<symbol>(e1.op(0).op(0)) &&
	      is_exactly_a<symbol>(e1.op(0).op(1)) &&
	      is_exactly_a<symbol>(e1.op(1)) &&
	      e1.is_equal(power(power(x,a),b)) )) {
		clog << "(x^a)^b, x,a,b symbolic wrong" << endl;
		clog << "returned: " << e1 << endl;
		return 1;
	}
	
	ex e2 = e1.subs(a==1);
	if (!(is_exactly_a<power>(e2) &&
	      is_exactly_a<symbol>(e2.op(0)) &&
	      is_exactly_a<symbol>(e2.op(1)) &&
	      e2.is_equal(power(x,b)) )) {
		clog << "(x^a)^b, x,b symbolic, a==1 wrong" << endl;
		clog << "returned: " << e2 << endl;
		return 1;
	}
	
	ex e3 = e1.subs(a==-1);
	if (!(is_exactly_a<power>(e3) &&
	      is_exactly_a<power>(e3.op(0)) &&
	      is_exactly_a<symbol>(e3.op(0).op(0)) &&
	      is_exactly_a<numeric>(e3.op(0).op(1)) &&
	      is_exactly_a<symbol>(e3.op(1)) &&
	      e3.is_equal(power(power(x,-1),b)) )) {
		clog << "(x^a)^b, x,b symbolic, a==-1 wrong" << endl;
		clog << "returned: " << e3 << endl;
		return 1;
	}
	
	ex e4 = e1.subs(lst(a==-1, b==2.5));
	if (!(is_exactly_a<power>(e4) &&
	      is_exactly_a<power>(e4.op(0)) &&
	      is_exactly_a<symbol>(e4.op(0).op(0)) &&
	      is_exactly_a<numeric>(e4.op(0).op(1)) &&
	      is_exactly_a<numeric>(e4.op(1)) &&
	      e4.is_equal(power(power(x,-1),2.5)) )) {
		clog << "(x^a)^b, x symbolic, a==-1, b==2.5 wrong" << endl;
		clog << "returned: " << e4 << endl;
		return 1;
	}
	
	ex e5 = e1.subs(lst(a==-0.9, b==2.5));
	if (!(is_exactly_a<power>(e5) &&
	      is_exactly_a<symbol>(e5.op(0)) &&
	      is_exactly_a<numeric>(e5.op(1)) &&
	      e5.is_equal(power(x,numeric(-0.9)*numeric(2.5))) )) {
		clog << "(x^a)^b, x symbolic, a==-0.9, b==2.5 wrong" << endl;
		clog << "returned: " << e5 << endl;
		return 1;
	}
	
	ex e6 = e1.subs(lst(a==numeric(3)+numeric(5.3)*I, b==-5));
	if (!(is_exactly_a<power>(e6) &&
	      is_exactly_a<symbol>(e6.op(0)) &&
	      is_exactly_a<numeric>(e6.op(1)) &&
	      e6.is_equal(power(x,numeric(-15)+numeric(5.3)*numeric(-5)*I)) )) {
		clog << "(x^a)^b, x symbolic, a==3+5.3*I, b==-5 wrong" << endl;
		clog << "returned: " << e6 << endl;
		return 1;
	}
	
	return 0;
}

static unsigned exam_powerlaws2()
{
	// (a*x)^b = a^b * x^b
	
	symbol x("x");
	symbol a("a");
	symbol b("b");
	
	ex e1 = power(a*x,b);
	if (!(is_exactly_a<power>(e1) &&
	      is_exactly_a<mul>(e1.op(0)) &&
	      (e1.op(0).nops()==2) &&
	      is_exactly_a<symbol>(e1.op(0).op(0)) &&
	      is_exactly_a<symbol>(e1.op(0).op(1)) &&
	      is_exactly_a<symbol>(e1.op(1)) &&
	      e1.is_equal(power(a*x,b)) )) {
		clog << "(a*x)^b, x,a,b symbolic wrong" << endl;
		clog << "returned: " << e1 << endl;
		return 1;
	}
	
	ex e2 = e1.subs(a==3);
	if (!(is_exactly_a<power>(e2) &&
	      is_exactly_a<mul>(e2.op(0)) &&
	      (e2.op(0).nops()==2) &&
	      is_exactly_a<symbol>(e2.op(0).op(0)) &&
	      is_exactly_a<numeric>(e2.op(0).op(1)) &&
	      is_exactly_a<symbol>(e2.op(1)) &&
	      e2.is_equal(power(3*x,b)) )) {
		clog << "(a*x)^b, x,b symbolic, a==3 wrong" << endl;
		clog << "returned: " << e2 << endl;
		return 1;
	}
	
	ex e3 = e1.subs(b==-3);
	if (!(is_exactly_a<mul>(e3) &&
	      (e3.nops()==2) &&
	      is_exactly_a<power>(e3.op(0)) &&
	      is_exactly_a<power>(e3.op(1)) &&
	      e3.is_equal(power(a,-3)*power(x,-3)) )) {
		clog << "(a*x)^b, x,a symbolic, b==-3 wrong" << endl;
		clog << "returned: " << e3 << endl;
		return 1;
	}
	
	ex e4 = e1.subs(b==4.5);
	if (!(is_exactly_a<power>(e4) &&
	      is_exactly_a<mul>(e4.op(0)) &&
	      (e4.op(0).nops()==2) &&
	      is_exactly_a<symbol>(e4.op(0).op(0)) &&
	      is_exactly_a<symbol>(e4.op(0).op(1)) &&
	      is_exactly_a<numeric>(e4.op(1)) &&
	      e4.is_equal(power(a*x,4.5)) )) {
		clog << "(a*x)^b, x,a symbolic, b==4.5 wrong" << endl;
		clog << "returned: " << e4 << endl;
		return 1;
	}
	
	ex e5 = e1.subs(lst(a==3.2, b==3+numeric(5)*I));
	if (!(is_exactly_a<mul>(e5) &&
	      (e5.nops()==2) &&
	      is_exactly_a<power>(e5.op(0)) &&
	      is_exactly_a<numeric>(e5.op(1)) &&
	      e5.is_equal(power(x,3+numeric(5)*I)*
					  power(numeric(3.2),3+numeric(5)*I)) )) {
		clog << "(a*x)^b, x symbolic, a==3.2, b==3+5*I wrong" << endl;
		clog << "returned: " << e5 << endl;
		return 1;
	}
	
	ex e6 = e1.subs(lst(a==-3.2, b==3+numeric(5)*I));
	if (!(is_exactly_a<mul>(e6) &&
	      (e6.nops()==2) &&
	      is_exactly_a<power>(e6.op(0)) &&
	      is_exactly_a<numeric>(e6.op(1)) &&
	      e6.is_equal(power(-x,3+numeric(5)*I)*
					  power(numeric(3.2),3+numeric(5)*I)) )) {
		clog << "(a*x)^b, x symbolic, a==-3.2, b==3+5*I wrong" << endl;
		clog << "returned: " << e6 << endl;
		return 1;
	}
	
	ex e7 = e1.subs(lst(a==3+numeric(5)*I, b==3.2));
	if (!(is_exactly_a<power>(e7) &&
	      is_exactly_a<mul>(e7.op(0)) &&
	      (e7.op(0).nops()==2) &&
	      is_exactly_a<symbol>(e7.op(0).op(0)) &&
	      is_exactly_a<numeric>(e7.op(0).op(1)) &&
	      is_exactly_a<numeric>(e7.op(1)) &&
	      e7.is_equal(power((3+numeric(5)*I)*x,3.2)) )) {
		clog << "(a*x)^b, x symbolic, a==3+5*I, b==3.2 wrong" << endl;
		clog << "returned: " << e7 << endl;
		return 1;
	}
	
	return 0;
}

static unsigned exam_powerlaws3()
{
	// numeric evaluation

	ex e1 = power(numeric(4),numeric(1,2));
	if (e1 != 2) {
		clog << "4^(1/2) wrongly returned " << e1 << endl;
		return 1;
	}
	
	ex e2 = power(numeric(27),numeric(2,3));
	if (e2 != 9) {
		clog << "27^(2/3) wrongly returned " << e2 << endl;
		return 1;
	}
	
	ex e3 = power(numeric(5),numeric(1,2));
	if (!(is_exactly_a<power>(e3) &&
	      e3.op(0).is_equal(numeric(5)) &&
	      e3.op(1).is_equal(numeric(1,2)))) {
		clog << "5^(1/2) wrongly returned " << e3 << endl;
		return 1;
	}
	
	ex e4 = power(numeric(5),evalf(numeric(1,2)));
	if (!(is_exactly_a<numeric>(e4))) {
		clog << "5^(0.5) wrongly returned " << e4 << endl;
		return 1;
	}
	
	ex e5 = power(evalf(numeric(5)),numeric(1,2));
	if (!(is_exactly_a<numeric>(e5))) {
		clog << "5.0^(1/2) wrongly returned " << e5 << endl;
		return 1;
	}
	
	return 0;
}

static unsigned exam_powerlaws4()
{
	// test for mul::eval()
	
	symbol a("a");
	symbol b("b");
	symbol c("c");
	
	ex f1 = power(a*b,ex(1)/ex(2));
	ex f2 = power(a*b,ex(3)/ex(2));
	ex f3 = c;
	
	exvector v;
	v.push_back(f1);
	v.push_back(f2);
	v.push_back(f3);
	ex e1 = mul(v);
	if (e1!=a*a*b*b*c) {
		clog << "(a*b)^(1/2)*(a*b)^(3/2)*c wrongly returned " << e1 << endl;
		return 1;
	}
	
	return 0;
}

static unsigned exam_powerlaws5()
{
	// cabinet of slightly pathological cases
	
	symbol a("a");
	
	ex e1 = pow(1,a);
	if (e1 != 1) {
		clog << "1^a wrongly returned " << e1 << endl;
		return 1;
	}
	
	ex e2 = pow(0,a);
	if (!(is_exactly_a<power>(e2))) {
		clog << "0^a was evaluated to " << e2
		     << " though nothing is known about a." << endl;
		return 1;
	}
	
	return 0;
}

unsigned exam_powerlaws()
{
	unsigned result = 0;
	
	cout << "examining power laws" << flush;
	
	result += exam_powerlaws1();  cout << '.' << flush;
	result += exam_powerlaws2();  cout << '.' << flush;
	result += exam_powerlaws3();  cout << '.' << flush;
	result += exam_powerlaws4();  cout << '.' << flush;
	result += exam_powerlaws5();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_powerlaws();
}
