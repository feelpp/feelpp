/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Mar 2015

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
#ifndef FEELPP_MODELMARKERS_HPP
#define FEELPP_MODELMARKERS_HPP 1

#include <feel/feelmodels/modelindexes.hpp>

namespace Feel {

class FEELPP_EXPORT ModelMarkers : public std::set<std::string>
{
    using super_type = std::set<std::string>;
  public:
    ModelMarkers() = default;
    ModelMarkers( ModelMarkers const& ) = default;
    ModelMarkers( ModelMarkers&& ) = default;
    explicit ModelMarkers( super_type const& set)
        : super_type(set) {}
    explicit ModelMarkers( std::string const& marker );
    ModelMarkers& operator=( ModelMarkers const& ) = default;
    ModelMarkers& operator=( ModelMarkers&& ) = default;
    void setPTree( pt::ptree const& p, ModelIndexes const& indexes = ModelIndexes() );
    void setup( nl::json const& jarg, ModelIndexes const& indexes = ModelIndexes() );
    bool empty() const;
};

}

#endif
