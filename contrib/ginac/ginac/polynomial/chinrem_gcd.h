/** @file chinrem_gcd.h
 *
 *  Interface to GCD functions using Chinese remainder algorithm. */

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

#ifndef GINAC_CHINREM_GCD_H
#define GINAC_CHINREM_GCD_H

#include "ex.h"

namespace GiNaC {

extern ex chinrem_gcd(const ex& A_, const ex& B_, const exvector& vars);
extern ex chinrem_gcd(const ex& A, const ex& B);

struct chinrem_gcd_failed
{
	virtual ~chinrem_gcd_failed() { }
};

} // namspace GiNaC

#endif // ndef GINAC_CHINREM_GCD_H
