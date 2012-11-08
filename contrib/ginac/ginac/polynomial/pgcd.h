/** @file pgcd.h
 *
 *  Interface to GCD functions for polynomials over prime fields. */

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

#ifndef GINAC_CHINREM_GCD_PGCD_H
#define GINAC_CHINREM_GCD_PGCD_H

#include "ex.h"

namespace GiNaC {

/// Exception to be thrown when modular GCD algorithm fails
struct pgcd_failed
{
	virtual ~pgcd_failed() { }
};

/**
 * @brief Compute the GCD of two polynomials over a prime field Z_p
 * 
 * @param vars variables
 * @param p designates the coefficient field Z_p
 * @param A polynomial \in Z_p[vars]
 * @param B second polynomial \in Z_p[vars]
 */
extern ex
pgcd(const ex& A, const ex& B, const exvector& vars, const long p);

} // namespace GiNaC

#endif // ndef GINAC_CHINREM_GCD_PGCD_H
