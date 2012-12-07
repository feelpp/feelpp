/** @file gcd_uvar.cpp
 *
 *  Several GCD algorithms. */

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

#include "upoly.h"
#include "sr_gcd_uvar.h"
#include "heur_gcd_uvar.h"

#include <stdexcept>

namespace GiNaC {

upoly sr_gcd(const upoly& a, const upoly& b)
{
	upoly g;
	bool found = sr_gcd_priv(g, a, b);
	if (found)
		return g;

	throw std::runtime_error("failed to compute gcd");
}

bool heur_gcd_z(upoly& g, const upoly& a, const upoly& b)
{
	return heur_gcd_z_priv(g, a, b);
}

upoly pseudoremainder(const upoly& a, const upoly& b)
{
	upoly r;
	pseudoremainder(r, a, b);
	return r;

}

} // namespace GiNaC
