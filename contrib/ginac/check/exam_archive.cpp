/** @file exam_archive.cpp
 *
 *  Here we test GiNaC's archiving system. */

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

#include <fstream>
#include <iostream>
using namespace std;

unsigned exam_archive()
{
	unsigned result = 0;
	
	cout << "examining archiving system" << flush;

	symbol x("x"), y("y"), mu("mu"), dim("dim", "\\Delta");
	ex e, f;

	// This expression is complete nonsense but it contains every type of
	// GiNaC object
	e = -42 * x * pow(y, sin(y*Catalan)) * dirac_ONE()
	    * epsilon_tensor(idx(fail(), 3), idx(0, 3), idx(y/2, 3))
	  + lorentz_g(
	      varidx(lst(x, -11*y, acos(2*x).series(x==3-5*I, 3)) * color_ONE()
	        * metric_tensor(varidx(log(cos(128.0/(x*y))), 5), varidx(2, 5)), zeta(3)),
	      varidx(diag_matrix(lst(-1, Euler, atan(x/y==-15*I/17)))
	        * delta_tensor(idx(x, 2), idx(wild(7), 3)), zeta(3), true),
	      true
	    )
	  + dirac_gamma(varidx(mu, dim)) * dirac_gamma(varidx(mu, 4-dim, true))
	    * color_T(idx(x, 8), 1) * color_h(idx(x, 8), idx(y, 8), idx(2, 8))
	    * indexed(x, sy_anti(), idx(2*y+1, x), varidx(-mu, 5))
	  - 2.4275 * spinor_metric(spinidx(0, 2, false, true), spinidx(y))
	  + abs(x).series(x == y, 4);

	archive ar;
	ar.archive_ex(e, "expr 1");
	{
		std::ofstream fout("exam.gar", std::ios_base::binary);
		fout << ar;
	}
	ar.clear();
	{
		std::ifstream fin("exam.gar", std::ios_base::binary);
		fin >> ar;
	}
	f = ar.unarchive_ex(lst(x, y, mu, dim), "expr 1");

	ex difference = (f - e).expand();
	if (!difference.is_zero()) {
		clog << "archiving/unarchiving " << e << endl
		     << "erroneously returned " << f << endl;
		++result;
	}

	return result;
}

int main(int argc, char** argv)
{
	return exam_archive();
}
