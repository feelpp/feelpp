/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 12 Sep 2014

   Copyright (C) 2014 Feel++ Consortium

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/*
 * Initially Written by Soumyajit De (2013) for Shogun and
 * KRYLSTAT Copyright 2011 by Erlend Aune <erlenda@math.ntnu.no> under GPL2+
 *
 * few parts rewritten and adjusted for Feel++
 */

#define BOOST_TEST_MODULE jacobi elliptic functions testsuite

#include <testsuite/testsuite.hpp>
#include <feel/feelcore/feelcomplex.hpp>
#include <feel/feelmath/jacobiellipticfunctions.hpp>
using namespace Feel;

typedef float64_t Real;
typedef complex128_t Complex;

BOOST_AUTO_TEST_SUITE( Jacobiellipticfunctions )

BOOST_AUTO_TEST_CASE( test_ellipkkp )
{
	Real K, Kp;
	math::ellipkkp(0.9999999999, K, Kp);
	BOOST_CHECK_CLOSE(K, 1.57153044117197637775, 5E-13);
	BOOST_CHECK_CLOSE(Kp, 4.52953569660335286784, 5E-13);

	math::ellipkkp(10.1, K, Kp);
	BOOST_CHECK_CLOSE(K, 1.57079632679489655800, 5E-13);
	BOOST_CHECK_CLOSE(Kp, 33.11638016237679948972, 5E-13);

	math::ellipkkp(10.0, K, Kp);
	BOOST_CHECK_CLOSE(K, 1.57079632679489655800, 5E-13);
	BOOST_CHECK_EQUAL(Kp, std::numeric_limits<Real>::max());

	math::ellipkkp(0.0, K, Kp);
	BOOST_CHECK_EQUAL(K, std::numeric_limits<Real>::max());
	BOOST_CHECK_CLOSE(Kp, 1.5707963267948965580, 5E-13);
}

BOOST_AUTO_TEST_CASE( test_ellipjc )
{
	Complex sn, cn, dn;
	math::ellipjc(Complex(0.5,0.0), 0.0, sn, cn, dn);
	BOOST_CHECK_CLOSE(sn.real(), 0.47942553860420300538, 5E-13);
	BOOST_CHECK_CLOSE(cn.real(), 0.87758256189037275874, 5E-13);
	BOOST_CHECK_CLOSE(dn.real(), 1.0, 5E-13);
	BOOST_CHECK_CLOSE(sn.imag(), 0.0, 5E-13);
	BOOST_CHECK_CLOSE(cn.imag(), 0.0, 5E-13);
	BOOST_CHECK_CLOSE(dn.imag(), 0.0, 5E-13);

	math::ellipjc(Complex(0.5,0.5), 0.0, sn, cn, dn);
	BOOST_CHECK_CLOSE(sn.real(), 0.54061268571315335141, 5E-13);
	BOOST_CHECK_CLOSE(cn.real(), 0.98958488339991990124, 5E-13);
	BOOST_CHECK_CLOSE(dn.real(), 1.0, 5E-13);
	BOOST_CHECK_CLOSE(sn.imag(), 0.45730415318424921800, 5E-13);
	BOOST_CHECK_CLOSE(cn.imag(), -0.24982639750046153893, 5E-13);
	BOOST_CHECK_CLOSE(dn.imag(), 0.0, 5E-13);

	math::ellipjc(Complex(0.5,0.5), 0.5, sn, cn, dn);
	BOOST_CHECK_CLOSE(sn.real(), 0.55284325098354925032, 5E-13);
	BOOST_CHECK_CLOSE(cn.real(), 0.96941944802337476350, 5E-8);//low precision
	BOOST_CHECK_CLOSE(dn.real(), 0.97703467549246159063, 5E-13);
	BOOST_CHECK_CLOSE(sn.imag(), 0.43032986483620772056, 5E-13);
	BOOST_CHECK_CLOSE(cn.imag(), -0.24540972636400421036, 5E-8);//low precision
	BOOST_CHECK_CLOSE(dn.imag(), -0.12174847394819815483, 5E-13);

	math::ellipjc(Complex(0.2,0.0), 0.99, sn, cn, dn);
	BOOST_CHECK_CLOSE(sn.real(), 0.19738823726686019477, 5E-13);
	BOOST_CHECK_CLOSE(cn.real(), 0.98032539689058428856, 5E-7);//low precision
	BOOST_CHECK_CLOSE(dn.real(), 0.98052409707808552142, 5E-13);
	BOOST_CHECK_CLOSE(sn.imag(), 0.0, 5E-13);
	BOOST_CHECK_CLOSE(cn.imag(), 0.0, 5E-13);
	BOOST_CHECK_CLOSE(dn.imag(), 0.0, 5E-13);

	math::ellipjc(Complex(0.7,0.4),
		0.9999999999999999, sn, cn, dn);
	BOOST_CHECK_CLOSE(sn.real(), 0.66873789889962576005, 5E-13);
	BOOST_CHECK_CLOSE(cn.real(), 0.81197156278765458826, 5E-13);
	BOOST_CHECK_CLOSE(dn.real(), 0.81197156278765469928, 5E-13);
	BOOST_CHECK_CLOSE(sn.imag(), 0.25191557357137611683, 5E-13);
	BOOST_CHECK_CLOSE(cn.imag(), -0.20747708305429035658, 5E-13);
	BOOST_CHECK_CLOSE(dn.imag(), -0.20747708305429032882, 5E-13);

	math::ellipjc(Complex(0.5,0.5), 1.0-1E-16,
		sn, cn, dn);
	BOOST_CHECK_CLOSE(sn.real(), 0.56408314126749847794, 5E-13);
	BOOST_CHECK_CLOSE(cn.real(), 0.94997886761549465984, 5E-13);
	BOOST_CHECK_CLOSE(dn.real(), 0.94997886761549465984, 5E-13);
	BOOST_CHECK_CLOSE(sn.imag(), 0.40389645531602574868, 5E-13);
	BOOST_CHECK_CLOSE(cn.imag(), -0.23982763093808803778, 5E-13);
	BOOST_CHECK_CLOSE(dn.imag(), -0.23982763093808801003, 5E-13);
}

BOOST_AUTO_TEST_SUITE_END()
