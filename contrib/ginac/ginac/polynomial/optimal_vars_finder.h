/** @file optimal_vars_finder.h
 *
 *  Interface to a function that optimizes the choice of variable for GCD
 *  computations. */

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

#ifndef GINAC_CHINREM_GCD_OPTIMAL_SYMBOL_FINDER_H
#define GINAC_CHINREM_GCD_OPTIMAL_SYMBOL_FINDER_H

#include "ex.h"

#include <vector>

namespace GiNaC {

/**
 * @brief Find the order of variables which is optimal for GCD computation.
 *
 * Collects statistical information about the highest and lowest degrees
 * of all variables that appear in two polynomials. Sorts the variables
 * by minimum degree (lowest to highest). The information gathered by
 * this function is used by GCD routines to find out the main variable
 * for GCD computation.
 */
extern exvector gcd_optimal_variables_order(const ex& A, const ex& B);

} // namespace GiNaC

#endif // GINAC_CHINREM_GCD_OPTIMAL_SYMBOL_FINDER_H
