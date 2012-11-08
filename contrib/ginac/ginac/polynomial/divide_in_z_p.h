/** @file divide_in_z_p.h
 *
 *  Interface to polynomial division in Z/Zp. */

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

#ifndef GINAC_CHINREM_GCD_DIVIDE_IN_Z_P_H
#define GINAC_CHINREM_GCD_DIVIDE_IN_Z_P_H

#include "ex.h"

namespace GiNaC {

/** 
 * Exact polynomial division of a, b \in Z_p[x_0, \ldots, x_n]
 * It doesn't check whether the inputs are proper polynomials, so be careful
 * of what you pass in.
 *  
 * @param a  first multivariate polynomial (dividend)
 * @param b  second multivariate polynomial (divisor)
 * @param q  quotient (returned)
 * @param var variables X iterator to first element of vector of symbols
 *
 * @return "true" when exact division succeeds (the quotient is returned in
 *          q), "false" otherwise.
 */
extern bool
divide_in_z_p(const ex &a, const ex &b, ex &q, const exvector& vars, const long p);

} // namespace GiNaC

#endif // ndef GINAC_CHINREM_GCD_DIVIDE_IN_Z_P_H
