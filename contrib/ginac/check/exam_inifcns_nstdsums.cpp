/** @file exam_inifcns_nstdsums.cpp
 *
 *  This test routine applies assorted tests on initially known higher level
 *  functions. */

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
#include <fstream>
using namespace std;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  S exam
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/*
 * The data in the following include file has been produced by the following
 * Mathematica (V4.1) script:
 *
 *
 *    x={2/10,1,14/10,30/10}
 *    y={0,3/10,-14/10}
 *    st = OpenAppend["exam_inifcns_nstdsums_data.raw"]
 *    $NumberMarks = False
 *    Do[
 *      Do[
 *        Do[Write[st, i]; Write[st,j]; Write[st,x[[k]]+I*y[[l]]];
 *          Write[st,Chop[N[PolyLog[i,j,x[[k]]+I*y[[l]]],25]]],{i,3},{j,3}], {k,4}],{l,3}]
 *    Do[
 *      Do[
 *        Do[Write[st, i]; Write[st,j]; Write[st,-x[[k]]+I*y[[l]]];
 *          Write[st,Chop[N[PolyLog[i,j,-x[[k]]+I*y[[l]]],25]]],{i,3},{j,3}], {k,4}], {l,3}]
 *    Close[st]
 *
 *    
 * and postprocessed by the following shell script
 *
 *
 *    #/bin/sh
 *    IFS=$'\n'
 *    cat exam_inifcns_nstdsums_data.raw | sed -e 's/\*\^/E/g' > exam_inifcns_nstdsums_data.raw2
 *    echo 'const char *data[] = {' > exam_inifcns_nstdsums_data.raw3
 *    for i in `cat exam_inifcns_nstdsums_data.raw2`; do echo \"$i\",; done >> exam_inifcns_nstdsums_data.raw3
 *    echo '"-999"};' >> exam_inifcns_nstdsums.h
 *
 *
 */
#include "exam_inifcns_nstdsums.h"


// signals end of data
const int ENDMARK = -999;


static unsigned inifcns_test_S()
{
	int digitsbuf = Digits;
	// precision of data
	Digits = 22;
	ex prec = 5 * pow(10, -(ex)Digits);
	
	unsigned result = 0;
	
	int i = 0;
	while (true) {
		ex n(data[i++],symbol());
		if (n == ENDMARK) {
			break;
		}
		ex p(data[i++],symbol());
		ex x(data[i++],symbol());
		ex res(data[i++],symbol());
		ex res2 = S(n, p, x).evalf();
		if (abs(res-res2) > prec) {
			clog << "S(" << n << "," << p << "," << x << ") seems to be wrong:" << endl;
			clog << "GiNaC           : " << res2 << endl;
			clog << "Reference       : " << res << endl;
			clog << "Abs. Difference : " << res2-res << endl;
			if (res2 != 0) {
				ex reldiff = abs((res2-res)/res2);
				clog << "Rel. Difference : " << reldiff << endl;
			}
			result++;
		}
		if (i % 80) {
			cout << "." << flush;
		}
	}

	Digits = digitsbuf;

	return result;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  H/Li exam
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


static unsigned inifcns_test_HLi()
{
	using GiNaC::log;
	int digitsbuf = Digits;
	Digits = 17;
	ex prec = 5 * pow(10, -(ex)Digits);
	numeric almostone("0.999999999999999999");
	unsigned result = 0;

	lst res;
	
	res.append(H(lst(2,1),numeric(1)/2).hold() - (zeta(3)/8 - pow(log(2),3)/6));
	res.append(H(lst(2,1,3),numeric(1)/3).hold() - Li(lst(2,1,3),lst(numeric(1)/3,1,1)).hold());
	res.append(H(lst(2,1,3),numeric(98)/100).hold() - Li(lst(2,1,3),lst(numeric(98)/100,1,1)).hold());
	res.append(H(lst(2,1,3),numeric(245)/100).hold() - Li(lst(2,1,3),lst(numeric(245)/100,1,1)).hold());
	res.append(H(lst(4,1,1,1),numeric(1)/3).hold() - S(3,4,numeric(1)/3).hold());
	res.append(H(lst(4,1,1,1),numeric(98)/100).hold() - S(3,4,numeric(98)/100).hold());
	res.append(H(lst(4,1,1,1),numeric(245)/100).hold() - S(3,4,numeric(245)/100).hold());
	res.append(H(lst(2,2,3),almostone).hold() - zeta(lst(2,2,3)));
	res.append(H(lst(-3,-1,2,1),almostone).hold() - zeta(lst(3,1,2,1),lst(-1,1,-1,1)));
	res.append(H(lst(-2,1,3),numeric(1)/3).hold() - -Li(lst(2,1,3),lst(-numeric(1)/3,-1,1)).hold());
	res.append(H(lst(-2,1,3),numeric(98)/100).hold() - -Li(lst(2,1,3),lst(-numeric(98)/100,-1,1)).hold());
	res.append(H(lst(-2,1,3),numeric(245)/100).hold() - -Li(lst(2,1,3),lst(-numeric(245)/100,-1,1)).hold());
	res.append(H(lst(-3,1,-2,0,0),numeric(3)/10).hold() - convert_H_to_Li(lst(-3,1,-2,0,0),numeric(3)/10).eval());
	
	for (lst::const_iterator it = res.begin(); it != res.end(); it++) {
		ex diff = abs((*it).evalf());
		if (diff > prec) {
			clog << *it << " seems to be wrong: " << diff << endl;
			result++;
		}
		cout << "." << flush;
	}

	Digits = digitsbuf;

	// conjugate test
	numeric cdif = ex_to<numeric>(H(lst(2,2,1),5.0-5.0*I) - H(lst(2,2,1),5.0+5.0*I));
	numeric cadd = ex_to<numeric>(H(lst(2,2,1),5.0-5.0*I) + H(lst(2,2,1),5.0+5.0*I));
	if ((cdif.real() > prec) || (cadd.imag() > prec)) {
		clog << "complex conjugation test of H({2,2,1},5.0-5.0*I) seems to be wrong: " << cdif << " " << cadd << endl;
		result++;
	}

	return result;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  zeta exam
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


static unsigned inifcns_test_zeta()
{
	int digitsbuf = Digits;
	
	unsigned result = 0;

	lst res;
	
	res.append(zeta(lst(2,1)) - zeta(3));
	res.append(zeta(lst(2,1,1,1,1)) - zeta(6));
	res.append(zeta(lst(6,3)) - (zeta(9)*83/2 - zeta(2)*zeta(7)*21 - zeta(2)*zeta(2)*zeta(5)*12/5));
	res.append(zeta(lst(4,2,3)) - (-zeta(9)*59 + zeta(2)*zeta(7)*28 + pow(zeta(2),2)*zeta(5)*4 -
	                               pow(zeta(3),3)/3 + pow(zeta(2),3)*zeta(3)*8/21));
	res.append(zeta(lst(3,1,3,1,3,1,3,1)) - (2*pow(Pi,16)/factorial(18)));
	res.append(zeta(lst(2),lst(-1)) - -zeta(2)/2);
	res.append(zeta(lst(1,2),lst(-1,1)) - (-zeta(3)/4 - zeta(lst(1),lst(-1))*zeta(2)/2));
	res.append(zeta(lst(2,1,1),lst(-1,-1,1)) - (-pow(zeta(2),2)*23/40 - pow(zeta(lst(1),lst(-1)),2)*zeta(2)*3/4
	                                            - zeta(lst(3,1),lst(-1,1))*3/2 - zeta(lst(1),lst(-1))*zeta(3)*21/8));
	
	for (lst::const_iterator it = res.begin(); it != res.end(); it++) {
		Digits = 17;
		ex prec = 5 * pow(10, -(ex)Digits);
		ex diff = abs((*it).evalf());
		if (diff > prec) {
			clog << *it << " seems to be wrong: " << diff << endl;
			clog << "Digits: " << Digits << endl;
			result++;
		}
		cout << "." << flush;
		Digits = 40;
		prec = 5 * pow(10, -(ex)Digits);
		diff = abs((*it).evalf());
		if (diff > prec) {
			clog << *it << " seems to be wrong: " << diff << endl;
			clog << "Digits: " << Digits << endl;
			result++;
		}
		cout << "." << flush;
	}

	Digits = digitsbuf;

	return result;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  H/Li exam
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


static unsigned inifcns_test_LiG()
{
	int digitsbuf = Digits;
	Digits = 17;
	ex prec = 5 * pow(10, -(ex)Digits);
	numeric almostone("0.99999999999999999999");
	unsigned result = 0;

	lst res;
	
	res.append(Li(lst(4), lst(6)).hold() - Li(4, 6.0));
	res.append(G(lst(0,0,5.0,0,2.0,0,0,0,3.0),0.5).hold()
	           + Li(lst(3,2,4), lst(numeric(1,10), numeric(5,2), numeric(2,3))));
	res.append(Li(lst(2,1,1), lst(almostone, almostone, almostone)) - zeta(lst(2,1,1)));

	// check Li_{1,1} against known expression
	symbol x("x"), y("y");
	ex eps = 1e-30*I;
	ex s1 = Li(lst(1,1),lst(x,y));
	ex s2 = log(1-1/x/y-eps)*log((1-1/x-eps)/(1/x/y-1/x)) + Li(2,(1-1/x/y-eps)/(1/x-1/x/y))
			- log(-1/x/y-eps)*log((-1/x-eps)/(1/x/y-1/x)) - Li(2,(-1/x/y-eps)/(1/x-1/x/y))
			- log(-1/x/y-eps)*log(1-1/x-eps) + log(-1/x/y-eps)*log(-1/x-eps);
	res.append(s1.subs(lst(x==numeric(1)/2, y==3)) - s2.subs(lst(x==numeric(1)/2, y==3)));
	res.append(s1.subs(lst(x==numeric(3)/2, y==numeric(1)/2)) - s2.subs(lst(x==numeric(3)/2, y==numeric(1)/2)));
	res.append(s1.subs(lst(x==2, y==numeric(4)/5)) - s2.subs(lst(x==2, y==numeric(4)/5)));

	// shuffle and quasi-shuffle identities
	res.append(G(lst(0,0.2),1).hold() * G(lst(0.5),1).hold() - G(lst(0.5,0,0.2),1).hold()
			- G(lst(0,0.5,0.2),1).hold() - G(lst(0,0.2,0.5),1).hold());
	res.append(G(lst(0,0.5),1).hold() * G(lst(0.6),1).hold() - G(lst(0,0.5,0.5*0.6),1).hold()
			- G(lst(0.6,0,0.5*0.6),1).hold() + G(lst(0,0,0.5*0.6),1).hold());
	res.append(Li(lst(2),lst(numeric(1,5))).hold() * Li(lst(3),lst(7)).hold() - Li(lst(2,3),lst(numeric(1,5),7)).hold()
			- Li(lst(3,2),lst(7,numeric(1,5))).hold() - Li(lst(5),lst(numeric(7,5))).hold());
	symbol a1, a2, a3, a4;
	res.append((G(lst(a1,a2),1) * G(lst(a3,a4),1) - G(lst(a1,a2,a3,a4),1)
			- G(lst(a1,a3,a2,a4),1) - G(lst(a3,a1,a2,a4),1)
			- G(lst(a1,a3,a4,a2),1) - G(lst(a3,a1,a4,a2),1) - G(lst(a3,a4,a1,a2),1))
				.subs(lst(a1==numeric(1)/10, a2==numeric(3)/10, a3==numeric(7)/10, a4==5)));
	res.append(G(lst(-0.009),1).hold() * G(lst(-8,1.4999),1).hold() - G(lst(-0.009,-8,1.4999),1).hold()
			- G(lst(-8,-0.009,1.4999),1).hold() - G(lst(-8,1.4999,-0.009),1).hold());
	res.append(G(lst(sqrt(numeric(1)/2)+I*sqrt(numeric(1)/2)),1).hold() * G(lst(1.51,-0.999),1).hold()
			- G(lst(sqrt(numeric(1)/2)+I*sqrt(numeric(1)/2),1.51,-0.999),1).hold()
			- G(lst(1.51,sqrt(numeric(1)/2)+I*sqrt(numeric(1)/2),-0.999),1).hold()
			- G(lst(1.51,-0.999,sqrt(numeric(1)/2)+I*sqrt(numeric(1)/2)),1).hold());
	// checks for hoelder convolution which is used if one argument has a distance to one smaller than 0.01 
	res.append(G(lst(0, 1.2, 1, 1.01), 1).hold() - G(lst(0, 1.2, 1, numeric("1.009999999999999999")), 1).hold());

	for (lst::const_iterator it = res.begin(); it != res.end(); it++) {
		ex diff = abs((*it).evalf());
		if (diff > prec) {
			clog << *it << " seems to be wrong: " << diff << endl;
			result++;
		}
		cout << "." << flush;
	}

	return result;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  legacy exam - checking for historical bugs
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


static unsigned inifcns_test_legacy()
{
	Digits = 17;
	ex prec = 5 * pow(10, -(ex)Digits);

	unsigned result = 0;

	ex r1 = zeta(lst(1,1,1,1,1,1),lst(-1,-1,-1,1,1,1));
	if ((r1.evalf() - numeric("-0.0012588769028204890704")) > prec) {
		clog << "zeta({1,1,1,1,1,1},{-1,-1,-1,1,1,1}) seems to be wrong." << endl;
		result++;
	}

	ex x1 = exp(2*Pi*I/13).evalf();
	ex x2 = exp(24*Pi*I/13).evalf();
  	ex r2 = Li(lst(2),lst(x1)).hold().evalf();
	ex r3 = Li(lst(2),lst(x2)).hold().evalf();
	if ( abs(r2-conjugate(r3)) > prec ) {
		clog << "Legacy test 2 seems to be wrong." << endl;
		result++;
	}

  	ex x3 = exp(5*Pi*I/3).evalf();
	ex r4 = Li(lst(3),lst(x3)).hold().evalf();
	if ( abs(r4 - numeric("0.40068563438653142847-0.95698384815740185713*I")) > prec ) {
		clog << "Legacy test 3 seems to be wrong." << endl;
		result++;
	}

	Digits = 100;
	prec = 5 * pow(10, -(ex)Digits);
	ex x0 = 1.;
	   x1 = exp(Pi*I/3).evalf();
	   x2 = exp(2*Pi*I/3).evalf();
	   x3 = -1.;
	ex x4 = exp(4*Pi*I/3).evalf();
	ex x5 = exp(5*Pi*I/3).evalf();

	ex r5 = Li(lst(1,1,1,1),lst(x2,x4,x3,x0)).hold().evalf();
	ex r6 = Li(lst(1,1,1,1),lst(x4,x2,x3,x0)).hold().evalf();
	if ( abs(r5-conjugate(r6)) > prec ) {
		clog << "Legacy test 4 seems to be wrong." << endl;
		result++;
	}

	ex r7 = Li(lst(1,2,1),lst(x3,x2,x4)).hold().evalf()
		+Li(lst(1,1,2),lst(x3,x2,x4)).hold().evalf()
		+Li(lst(1,1,1,1),lst(x3,x0,x2,x4)).hold().evalf()
		+Li(lst(1,1,1,1),lst(x3,x2,x0,x4)).hold().evalf()
		+Li(lst(1,1,1,1),lst(x3,x2,x4,x0)).hold().evalf()
		+Li(lst(1,2,1),lst(x2,x1,x0)).hold().evalf()
		+Li(lst(1,1,2),lst(x2,x3,x4)).hold().evalf()
		+Li(lst(1,1,1,1),lst(x2,x4,x3,x0)).hold().evalf()
		+Li(lst(1,1,1,1),lst(x2,x3,x4,x0)).hold().evalf()
		+Li(lst(1,1,1,1),lst(x2,x3,x0,x4)).hold().evalf()
		+Li(lst(2,2),lst(x5,x4)).hold().evalf()
		+Li(lst(2,1,1),lst(x5,x0,x4)).hold().evalf()
		+Li(lst(2,1,1),lst(x5,x4,x0)).hold().evalf()
		-Li(lst(1,1),lst(x3,x0)).hold().evalf()*Li(lst(1,1),lst(x2,x4)).hold().evalf();
	if ( abs(r7) > prec ) {
		clog << "Legacy test 5 seems to be wrong." << endl;
		result++;
	}

	return result;
}


unsigned exam_inifcns_nstdsums(void)
{
	unsigned result = 0;
	
	cout << "examining consistency of nestedsums functions" << flush;
	
	result += inifcns_test_zeta();
	result += inifcns_test_S();
	result += inifcns_test_HLi();
	result += inifcns_test_LiG();
	result += inifcns_test_legacy();
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_inifcns_nstdsums();
}
