/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 27 sept. 2015

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

#ifndef FEELPP_PTREETOOLS_HPP
#define FEELPP_PTREETOOLS_HPP

#include <string>
#include <feel/feelcore/feelmacros.hpp>
#include <boost/property_tree/ptree.hpp>

namespace Feel
{
namespace pt =  boost::property_tree;
/**
 * remove c and c++ comments from string \p str_with_comments if any
 */
FEELPP_EXPORT std::string removeComments( std::string str_with_comments );

FEELPP_EXPORT void editPtreeFromOptions( pt::ptree& p, std::string const& prefix="" );

FEELPP_EXPORT void mergePtree( pt::ptree& pa, const pt::ptree& pb );

}
#endif
