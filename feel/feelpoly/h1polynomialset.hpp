/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 24 août 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_H1POLYNOMIALSET_HPP
#define FEELPP_H1POLYNOMIALSET_HPP 1

#include <boost/type_traits/is_base_of.hpp>

namespace Feel {

/**
 * H1 conforming polynomialset base class
 */
class H1PolynomialSet {};

/**
 * type traits for H1 conforming polynomialset
 * @return true_type if H1 conforming polynomialset, false_type otherwise
 */
template<typename P>
class is_h1_conforming : public boost::is_base_of<H1PolynomialSet,P>
{
};

}


#endif
