/** @file time_lw_C.cpp
 *
 *  Test C from the paper "Comparison of Polynomial-Oriented CAS" by Robert H.
 *  Lewis and Michael Wester. */

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
#include "timer.h"
using namespace GiNaC;

#include <iostream>
#include <vector>
using namespace std;

static unsigned test()
{
	numeric x(13*17*31);
	numeric y(13*19*29);
	
	for (int i=1; i<200; ++i)
		gcd(pow(x,300+(i%181)),pow(y,200+(i%183)));
	
	ex lastgcd = gcd(pow(x,300+(200%181)),pow(y,200+(200%183)));
	if (lastgcd != numeric("53174994123961114423610399251974962981084780166115806651505844915220196792416194060680805428433601792982500430324916963290494659936522782673704312949880308677990050199363768068005367578752699785180694630122629259539608472261461289805919741933")) {
		clog << "gcd(" << x << "^" << 300+(200%181) << ","
		     << y << "^" << 200+(200%183) << ") erroneously returned "
		     << lastgcd << endl;
		return 1;
	}
	return 0;
}

unsigned time_lw_C()
{
	unsigned result = 0;
	unsigned count = 0;
	timer rolex;
	double time = .0;
	
	cout << "timing Lewis-Wester test C (gcd of big integers)" << flush;
	
	rolex.start();
	// correct for very small times:
	do {
		result = test();
		++count;
	} while ((time=rolex.read())<0.1 && !result);
	cout << '.' << flush;
	cout << time/count << 's' << endl;
	
	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_lw_C();
}
