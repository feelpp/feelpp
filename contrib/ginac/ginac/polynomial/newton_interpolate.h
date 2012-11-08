/** @file newton_interpolate.h
 *
 *  Functions for Newton interpolation. */

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

#ifndef GINAC_PGCD_NEWTON_INTERPOLATE_H
#define GINAC_PGCD_NEWTON_INTERPOLATE_H

#include "ex.h"
#include "numeric.h"
#include "smod_helpers.h"

namespace GiNaC {

/**
 * Newton interpolation -- incremental form.
 *
 * Given a polynomials e1 \in Z_p[Y], prev \in Z_p[x][Y], and evaluation
 * points pt1, prevpts \in Z_p compute a polynomial P \in Z_[x][Y] such
 * that P(x = pt1) = e1, P(x = pt \in prevpts) = prev(x = pt).
 *
 * @var{prevpts} enumerates previous evaluation points, should have a form
 * (x - pt_2) (x - pt_3) \ldots (x - pt_n).
 * @var{prev} encodes the result of previous interpolations.
 */
ex newton_interp(const ex& e1, const long pt1,
		 const ex& prev, const ex& prevpts,
		 const ex& x, const long p)
{
	const numeric pnum(p);
	const numeric nc = ex_to<numeric>(prevpts.subs(x == numeric(pt1)).smod(pnum));
	const numeric nc_1 = recip(nc, p);
	// result = prev + prevpts (e1 - prev|_{x = pt_1})/prevpts|_{x = pt_1)
	ex tmp = prev.subs(x == numeric(pt1)).smod(pnum);
	tmp = (((e1 - tmp).expand().smod(pnum))*nc_1).smod(pnum);
	tmp = (prevpts*tmp).expand().smod(pnum);
	tmp = (prev + tmp).expand().smod(pnum);
	return tmp;
}

} // namespace GiNaC

#endif // ndef GINAC_PGCD_NEWTON_INTERPOLATE_H
